#sampling from cauchy distribution, with proposal density function as "normal"
n=1000000
x=rep(0, l=n)
x[1]=rnorm(1)
for (j in 2:n){
  Y<-rnorm(1)
  alpha<-min(dcauchy(Y)*dnorm(x[j-1])/(dcauchy(x[j-1])*dnorm(Y)),1)
  if(runif(1)<=alpha){
       x[j] <- Y
  }else{
       x[j] <- x[j-1]}
}
x
hist(x,50,freq=F,col="red",main="", xlim=c(-10,10)) 

xs = seq(-10,10,l=5000)
lines(xs, dcauchy(xs), type="l")

curve(dcauchy(x),from=-3,to=3,add=T,col="black")

acf(x)

thin_ind = seq(1,length(x), by = 50)

acf( x[thin_ind] )

hist( x[thin_ind], prob = TRUE, 50)
lines(xs, dcauchy(xs), type="l")

xx = rcauchy(100000)
hist( xx[abs(xx)<30], prob = TRUE, 1000, xlim=c(-10,10))
lines(xs, dcauchy(xs), type="l", lwd=2, col = 2)

#Below "random walk" 1/5/2024
#-----------------------------------------------------------------
n=1000000
x=rep(0, l=n)
x[1]=rnorm(1)
for (j in 2:n){
  Y<-rnorm(1, mean = x[j-1], sd = 1)
  alpha <- min( dcauchy(Y)/dcauchy(x[j-1]) , 1 )
  if(runif(1)<=alpha){
    x[j] <- Y
  }else{
    x[j] <- x[j-1]}
}
x
acf(x)
hist(x,50,freq=F,col="red",main="", xlim=c(-10,10)) 
xs = seq(-10,10,l=5000)
lines(xs, dcauchy(xs), type="l")


hist(x[abs(x)<30],100,freq=F,col="red",main="",xlim=c(-10,10)) 
xs = seq(-10,10,l=5000)
lines(xs, dcauchy(xs), type="l", lwd=2, col=3)

#1/5/2024-------------------------------------------------------------------
#MH for one poesterior parameter beta where sigma2=1 with proposal Multivariate normal distribution
rm(list=ls())
source("./dataGeneration.R")
#1) directly generate 10,000 beta value from posterior disbribution
sigma2 <- 1
gamma<-(solve(diag(3)+(1/sigma2)*t(x)%*%x)%*%t(x)%*%y)/sigma2
omega<-solve(diag(3)+(1/sigma2)*t(x)%*%x)
library(MASS)
sample_beta_1<-mvrnorm(n = 10000,mu = gamma,Sigma = omega) 
par(mfrow=c(1 ,3))
beta_1<-hist(sample_beta_1[,1],freq=F)
mean(sample_beta_1[,1])
acf(sample_beta_1[,1])
var(sample_beta_1[,1])
beta_2<-hist(sample_beta_1[,2],freq=F)
beta_3<-hist(sample_beta_1[,3],freq=F)

#2)MH sampler for beta values-----------1/11/2024
rm(list=ls())
source("./dataGeneration.R")
install.packages("mvtnorm")
install.packages("ggplot2")
library("mvtnorm")
set.seed(100)
burnin = 500
thin = 50
n = 5000*thin
beta_post_sample <- matrix(0, nrow = n/thin, ncol = 3)
beta_old = as.vector( beta )
counter = 1

sigma2 <- 1
gamma<-(solve(diag(3)+(1/sigma2)*t(x)%*%x)%*%t(x)%*%y)/sigma2
omega<-solve(diag(3)+(1/sigma2)*t(x)%*%x)


for (j in 1:(burnin + n) ){
  # proposal
  Y <- as.vector( mvrnorm(1, mu = beta_old, Sigma=diag(3)) )
  # density values
  den_new = dmvnorm(Y, mean=gamma, sigma=omega)
  den_old = dmvnorm(beta_old,mean=gamma,sigma=omega)
  # calculate acceptance probability
  alpha <- min( den_new/den_old, 1)
    
  if(runif(1)<=alpha){
    beta_new <- Y
  }else{
    beta_new <- beta_old
  }
  
  if( j > burnin & (j%%thin) == 0){
    # save new sample
    beta_post_sample[counter,] = beta_new
    counter = counter + 1
  }
  
  beta_old = beta_new
  
}

# Set figure size
options(repr.plot.width = 20, repr.plot.height = 15)

# Create a new plot area with 2 rows and 3 columns
par(mfrow = c(2, 3))

# Loop to create histograms for the first row
for (jj in 1:3) {
  # Set graphical parameters for smaller caption
  par(cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
  
  # Create histogram plot
  hist(beta_post_sample[, jj], prob = TRUE, main = paste0("Posterior beta", jj))
  
  # Add vertical lines at x=3, x=4, x=5
  abline(v = c(1, 2, 0.5), col = "red", lty = 1, lwd = c(2, 2, 2))
}

# Loop to create ACF plots for the second row
for (jj in 1:3) {
  # Set graphical parameters for smaller caption
  par(cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
  
  # Create ACF plot
  acf(beta_post_sample[, jj], main = paste0("Posterior beta", jj))
}

#1/11/2024-------------------------------------------------------------------------------
#MH for two posterior parameters beta and sigma2(in log transformed) with proposal truncated t distribution 
rm(list=ls())
source("./dataGeneration.R")

#create vector list for beta samplings and sigma2 samplings
burnin= 500
thin = 50
nmc = 5000*thin
beta_post_sample_1= matrix(0, nrow =nmc/thin, ncol = 3)
logsigma2_post_sample=rep(0, l=nmc/thin) #log-transformed sigma2

#computation to prepare for prior beta and sigma2
n=100 # sample sizes from data
alpha=0.1 # parameters for prior sigma2's inverse_gamma
beta_prime=0.1 # parameters for prior sigma2's inverse_gamma

#assign initial values for beta and sigma2 before MH_algorithm
beta_old_MH = as.vector( beta)
sigma2_old_MH= 1 #initial values for sigma2
logsigma2_old_MH<-log(sigma2_old_MH)#log transform for sigma2
counter=1
p=nrow(t(x)%*%x)

#MH algorithm for beta and logsigma2
for (j in 2:(burnin + nmc) ){
   # beta part
   # proposal from Multivariate student t distribution
  sigma2_old_MH<-exp(logsigma2_old_MH)
  omega <- solve( diag(3) + t(x)%*%x/sigma2_old_MH )
  gamma  <- omega %*% t(x) %*% y/sigma2_old_MH
  #set up the truncated range for the each sampling parameter
  Y_1 <- rtmvt(1, mean = beta_old_MH, sigma = diag(3), df = 3,
                         lower = rep(-Inf, length = p), upper = rep(Inf, length = p))
  Y_1 <- as.vector(Y_1)
  
  # density values: target density is the full conditional distribution of beta which is "Multivariate normal"
  den_new_1 = dmvnorm(Y_1, mean=gamma, sigma=omega)
  den_old_1 = dmvnorm(beta_old_MH,mean=gamma,sigma=omega)
  # calculate acceptance probability
  alpha <- min(den_new_1/den_old_1, 1) 
  
  if(runif(1)<=alpha){
    beta_new_MH <- Y_1 
  }else{
    beta_new_MH <- beta_old_MH
  }
  
  # sigma2 (log-transformed) part
  a = n/2 + alpha
  e = as.vector(y - x%*%as.matrix(beta_new_MH))
  b = as.numeric(crossprod(e)/2 + beta_prime)
  
  #proposal from univariate truncated student t distribution for log-transformed sigma2
  Y_2 <- rtmvt(1, mean = logsigma2_old_MH, sigma = diag(1), df = 3,
                              lower = rep(-Inf, length = 1), upper = rep(Inf, length = 1) ) 
  Y_2 <- as.vector(Y_2)
     
  den_new_2 = dinvgamma(exp(Y_2), shape=a, rate = b)*exp(Y_2) # add Jacobian adjustment term for the target density
  den_old_2 = dinvgamma(exp(logsigma2_old_MH),shape=a,rate = b)*exp(logsigma2_old_MH) # add Jacobian adjustment term for the target density
  # calculate acceptance probability
  alpha <- min(den_new_2/den_old_2, 1)
  
  if(runif(1)<=alpha){
    logsigma2_new_MH <- Y_2
  }else{
    logsigma2_new_MH <- logsigma2_old_MH
  }
  
  if( j > burnin & (j%%thin) == 0){
    # save new sample
    beta_post_sample_1[counter,] = beta_new_MH
    logsigma2_post_sample[counter]=logsigma2_new_MH
    counter = counter + 1
  }
  
  beta_old_MH = beta_new_MH
  logsigma2_old_MH = logsigma2_new_MH
  
}


par(mfrow=c(2,2))
for( jj in 1:3){
  hist( beta_post_sample_1[,jj], prob=TRUE, main = paste0("Posterior beta",jj))
}
hist((logsigma2_post_sample), prob=TRUE,main="logsigma2")

par(mfrow=c(2,2))
for( jj in 1:3){
  acf(beta_post_sample_1[,jj], main = paste0("Posterior beta",jj))
}
acf((logsigma2_post_sample), main = "logsigma2")



#--------------------------------------------------------------------------
dnorm(20)/dnorm(50)
#log scale for calculating the accpetance ratio: more stable
log( runif(1) ) <= ( dnorm(50, log=TRUE) - dnorm(-20, log=TRUE) )

#2/9/2024----------------------------------------------------------------------------
#MH for jointly beta and logsigma2
rm(list=ls())
source("./dataGeneration.R")

#computation to prepare for prior beta and sigma2
burnin= 500
thin = 50
nmc = 5000*thin
n=100 # sample sizes from data
alpha=0.1 # parameters for prior sigma2's inverse_gamma
beta_prime=0.1 # parameters for prior sigma2's inverse_gamma
beta_logsigma2_post_sample_1= matrix(0, nrow =nmc/thin, ncol = 4)
beta_logsigma2_post_1 <- matrix(c(1, 2, 0.5, 1), ncol = 1)

#assign initial values for beta and sigma2 before MH_algorithm
beta_logsigma2_old_MH = as.vector( beta_logsigma2_post_1 )
counter=1

#MH algorithm for beta and logsigma2 
for (j in 2:(burnin + nmc) ){
  # proposal from multivariate student t distribution 
  #set up the parameters for proposal multivariate_student t distribution
  Y_2 <- rtmvt(1, mean = beta_logsigma2_old_MH, sigma = diag(4), df = 3,
               lower = rep(-Inf, length = 4), upper = rep(Inf, length = 4))
  Y_2 <- as.vector(Y_2)
  
  # density values: target density--------------------------------------------
  e_new = y - x%*%as.matrix(Y_2[1:3])
  e_old = y - x%*%as.matrix(beta_logsigma2_old_MH[1:3])
  
  b_new = as.numeric( crossprod(e_new)/2*exp(Y_2[4]))
  b_old = as.numeric( crossprod(e_old)/2*exp(beta_logsigma2_old_MH[4]))
  
  c_new = as.numeric(crossprod(Y_2[1:3]))
  c_old = as.numeric(crossprod(beta_logsigma2_old_MH[4]))
  
  log_den_new = (-n/2-alpha+1)*Y_2[4] - b_new - c_new - (beta_prime/exp(Y_2[4]))
  log_den_old = (-n/2-alpha+1)*beta_logsigma2_old_MH[4] - b_old - c_old - (beta_prime/exp(beta_logsigma2_old_MH[4]))
  
  # calculate acceptance probability
  log_alpha <- (log_den_new - log_den_old) 
  
  if(log(runif(1))<=log_alpha){
    beta_logsigma2_new_MH <- Y_2 
  }else{
    beta_logsigma2_new_MH <- beta_logsigma2_old_MH
  }
  
  # burning and thin steps for adjust the convergence and autocorrelation in MCMC samplings
  if( j > burnin & (j%%thin) == 0){
    # save new sample
    beta_logsigma2_post_sample_1[counter,] = beta_logsigma2_new_MH
    counter = counter + 1
  }
  
  beta_logsigma2_old_MH = beta_logsigma2_new_MH
  
}

# sigma2 (log-transformed) part
a = n/2 + alpha
e = as.vector(y - x%*%as.matrix(beta_new_MH))
b = as.numeric( crossprod(e)/2 + beta_prime )

#proposal from univariate truncated student t distribution for log-transformed sigma2
#Y_2 <- rtt(1,location=logsigma2_old_MH, scale=1, df=1)#???#
Y_2 <- rtmvt(1, mean = logsigma2_old_MH, sigma = diag(1), df = 3,
             lower = rep(-Inf, length = 1), upper = rep(Inf, length = 1) ) #??? why not univarate truncated t disbribution?
Y_2 <- as.vector(Y_2)

# density values: target density is the full conditional distribution of sigma2 which is "inverse_gamma", since we did log transformed, so need to add jacobian adjustment
den_new_2 = dinvgamma(exp(Y_2), shape=a, rate = b)*exp(Y_2) # add jacobian adjustment term for the target density
den_old_2 = dinvgamma(exp(logsigma2_old_MH),shape=a,rate = b)*exp(logsigma2_old_MH) # add jacobian adjustment term for the target density

# calculate acceptance probability
alpha <- min(den_new_2/den_old_2, 1)

if(runif(1)<=alpha){
  logsigma2_new_MH <- Y_2
}else{
  logsigma2_new_MH <- logsigma2_old_MH
}



