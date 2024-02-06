#rm(list=ls())
library(tmvtnorm)
library("stats")
library(invgamma)#install.packages("invgamma")



########### 
burnin= 500
thin = 100
nmc = 5000*thin # the number of mcmc iterations

############### Data generation #################################
n <- 100 # sample size
set.seed(123123)
x1<-matrix(c(rep(1, n)), ncol=1)
x2_x3 <- matrix(rnorm(n*2), nrow = n, ncol = 2)
x <- cbind(x1, x2_x3)
eps<-matrix(rnorm(n),nrow =n, ncol=1)
beta <- matrix(c(1, 2, 0.5), ncol = 1)
y <- x%*%beta + eps



############## posterior sample will be stored ##################
beta_post_sample <-  matrix(0, nrow = nmc/thin, ncol = 3)
logsigma2_post_sample <- rep(0, l = nmc/thin) #log-transformed sigma2


# #create list for parameters gamma (vector) and omega(matrtix): parameters needed for target distribution of beta
# gamma= matrix(0, nrow = n/thin, ncol = 3)
# omega= replicate(n/thin ,matrix(0, nrow =3, ncol = 3))

# #create vector list for parameters a and b: parameters need for target distribution for sigma2
# a=rep(0, l=n+burnin)
# b=rep(0, l=n+burnin)

###################### Hyper-parameters
shape <- 0.01 # parameters for prior sigma2's inverse_gamma
rate <- 0.01 # parameters for prior sigma2's inverse_gamma
set.seed(1)


############## initial values ####################################
beta_old <-  as.vector( beta )
logsigma2_old <- 0
counter <- 1

#XtX <- t(x)%*%x
XtX <- crossprod(x) # crossprod is much faster than t(x)%*%x
XtY <- crossprod(x,y)
p <- nrow(XtX)
# XtX == t(x)%*%x
# XtY == t(x)%*%y

#MH algorithm for beta and logsigma2
##################### This is MH within Gibbs sampler #########################
####################### we run Gibbs steps, where each Gibbs step runs MH #####

for (j in 1:(burnin + nmc) ){
  # Sampling beta
  # (We don't need truncated distirbution as their support is (-inf,inf) )

  # This way, we need to calculate X'X and X'y for each iteration
  # That is inefficient because they are constant over ierations.
  # gamma[j-1,] <- solve(diag(3)+(1/sigma2_post_sample_1[j-1])*t(x)%*%x)%*%t(x)%*%y
  # omega[,,j-1] <- solve((diag(3)+(1/sigma2_post_sample_1[j-1])*t(x)%*%x))
  
  beta_proposed <- rtmvt(1, mean = beta_old, sigma = diag(3), df = 3,
                         lower = rep(-Inf, length = p), upper = rep(Inf, length = p) )
  beta_proposed <- as.vector(beta_proposed)
  
  omega <- solve( diag(3)+(1/sigma2_old)*XtX )
  gamma <- as.vector( omega %*% XtY/sigma2_old )
  
  # density values: target density is the full conditional distribution of beta which is "multivaraite normal"
  target_den_new <- dmvnorm(beta_proposed, mean=gamma,  sigma=omega )
  target_den_old <- dmvnorm(beta_old,      mean=gamma,  sigma=omega )
  
  # density values of proposal distribution
  # beta_new | beta_old (new given old)
  prop_den_new <- dtmvt(beta_proposed, mean = beta_old, sigma = diag(3), df = 3,
                        lower = rep(-Inf, length = p), upper = rep(Inf, length = p) )
  # beta_old | beta_new (old given new)
  prop_den_old <- dtmvt(beta_old, mean = beta_proposed, sigma = diag(3), df = 3,
                        lower = rep(-Inf, length = p), upper = rep(Inf, length = p) )
  
  # prop_den_new == prop_den_old
  # but I am just adding them for your reference
  numerator <- (target_den_new * prop_den_old)
  denominator <- (target_den_old * prop_den_new)
  # calculate acceptance probability
  alpha <- min( numerator/denominator, 1)
  
  if(runif(1)<=alpha){
    beta_new <- beta_proposed
  }else{
    beta_new <- beta_old
  }
  
  ###################### sigma2 (log-transformed) part
  a <- n/2 + shape
  e <- as.vector( y - x%*%beta_new ) # residual
  b <- sum(e^2)/2 + rate
  # note sum(e^2) = t(y - x%*% beta_new)%*%(y - x%*% beta_new)
  logsigma2_old <- log(sigma2_old)
  
  #proposal from univariate truncated student t distribution for log-transformed sigma2
  #Y_2 <- rtt(1,location=logsigma2_old, scale=1)
  
  logsigma2_proposed <- rtmvt(1, mean = logsigma2_old, sigma = diag(1), df = 3,
                           lower = rep(-Inf, length = 1), upper = rep(Inf, length = 1) )
  logsigma2_proposed <- as.vector( logsigma2_proposed )
                  
  # density values: target density is the full conditional distribution of sigma2 which is "inverse_gamma", since we did log transformed, so need to add jacobian adjustment
  logsig2_target_den_new <-  dinvgamma(exp(logsigma2_proposed), shape = a, rate = b)*exp(logsigma2_proposed) 
  logsig2_target_den_old <-  dinvgamma(exp(logsigma2_old), shape = a, rate = b)*exp(logsigma2_old) 
  # calculate acceptance probability
  
  logsig2_prop_den_new <- dmvt(logsigma2_proposed - logsigma2_old, sigma = diag(1), df = 3, log=FALSE )
  logsig2_prop_den_old <- dmvt(logsigma2_old - logsigma2_proposed, sigma = diag(1), df = 3, log=FALSE )
  
  
  # Again, logsig2_prop_den_new == logsig2_prop_den_old
  # but I am just adding them for your reference
  numerator <- (logsig2_target_den_new * logsig2_prop_den_old)
  denominator <- (logsig2_target_den_old * logsig2_prop_den_new)
  # calculate acceptance probability
  alpha <- min( numerator/denominator, 1)
  
  if(runif(1)<=alpha){
    logsigma2_new <- logsigma2_proposed
  }else{
    logsigma2_new <- logsigma2_old
  }
  
  # burning and thin steps for adjust the convergency and autocorrelation in MCMC samplings
  if( j > burnin & (j%%thin) == 0){
    # save new sample
    beta_post_sample[counter,] <- beta_new
    logsigma2_post_sample[counter] <- logsigma2_new
    counter = counter + 1
    if( counter%%10 == 0){ print(counter) }
  }

  beta_old <-  beta_new
  sigma2_old <- exp(logsigma2_new)
  
}


par(mfrow=c(2,2))
for( jj in 1:p){
  hist( beta_post_sample[,jj], prob=TRUE, main = paste0("Posterior beta",jj))
}
hist( exp(logsigma2_post_sample), prob=TRUE, main = "Posterior sigma2")


colMeans( beta_post_sample )
as.vector( beta )

par(mfrow=c(2,2))
for( jj in 1:p){
  acf(beta_post_sample[,3], main = paste0("Posterior beta",jj))
}
acf( logsigma2_post_sample, main = "Posterior sigma2"  )




