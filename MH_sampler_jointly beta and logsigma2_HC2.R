#MH for jointly beta and logsigma2
rm(list=ls())
source("./dataGeneration.R")

#computation to prepare for prior beta and sigma2
#p=4
burnin= 500
thin = 50
nmc = 5000*thin
n=100 # sample sizes from data
alpha=0.1 # parameters for prior sigma2's inverse_gamma
beta_prime=0.1 # parameters for prior sigma2's inverse_gamma
beta_logsigma2_post_sample_1= matrix(0, nrow =nmc/thin, ncol = 4) # sampling vectors
beta_logsigma2_post_1 <- matrix(c(1, 2, 0.5, 1), ncol = 1) #initial values

#assign initial values for beta and sigma2 before MH_algorithm
beta_logsigma2_old_MH = as.vector(beta_logsigma2_post_1)
counter=1

# param <- beta_logsigma2_old_MH
log_post_den <- function(param){
  
  e = y - x%*%as.matrix(param[1:3])
  b = as.numeric( crossprod(e)/( 2*exp(param[4])) )
  c = as.numeric(crossprod(param[1:3]))/2
  
  log_den = (-n/2-alpha)*( param[4] ) - b - c - (beta_prime/exp(param[4]))
  
  return(log_den)
  
}
param <- c(beta, log(2) )
optimum <- optim(param, log_post_den, control=list(fnscale=-1))#???? delta methods??

log_post_den(optimum$par)
#param <- optimum$par

sigma_prop <- -solve( numDeriv::hessian( log_post_den, optimum$par) )
# sigma_prop <- diag(4)


counter=1
#MH algorithm for beta and logsigma2 
for (j in 2:(burnin + nmc) ){
  # proposal from multivariate student t distribution 
  #set up the parameters for proposal multivariate_student t distribution
  Y_2 <- rtmvt(1, mean = beta_logsigma2_old_MH, sigma = sigma_prop, df = 5,
               lower = rep(-Inf, length = 4), upper = rep(Inf, length = 4))
  Y_2 <- as.vector(Y_2)
  
  # target density--------------------------------------------
  e_new = y - x%*%as.matrix(Y_2[1:3])
  e_old = y - x%*%as.matrix(beta_logsigma2_old_MH[1:3])
  
  b_new = as.numeric( crossprod(e_new)/( 2*exp(Y_2[4])) )
  b_old = as.numeric( crossprod(e_old)/( 2*exp(beta_logsigma2_old_MH[4])) )
  
  c_new = as.numeric(crossprod(Y_2[1:3]))/2
  c_old = as.numeric(crossprod(beta_logsigma2_old_MH[1:3]))/2
  
  log_den_new = (-n/2-alpha)*( Y_2[4] ) - b_new - c_new - (beta_prime/exp(Y_2[4]))
  log_den_old = (-n/2-alpha)*( beta_logsigma2_old_MH[4] ) - b_old - c_old - (beta_prime/exp(beta_logsigma2_old_MH[4]))
  
  # calculate acceptance probability in log scale
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

dim(beta_logsigma2_post_sample_1)
colMeans( beta_logsigma2_post_sample_1[,1:3] )
as.vector( beta )

mean( beta_logsigma2_post_sample_1[,4] )
plot( beta_logsigma2_post_sample_1[,4] )
hist( beta_logsigma2_post_sample_1[,4] )
acf( beta_logsigma2_post_sample_1[,4] )


par(mfrow=c(2,2))
for( jj in 1:p){
  hist( beta_logsigma2_post_sample_1[,jj], prob=TRUE, main = paste0("Posterior beta and logsigma2",jj))
}

par(mfrow=c(2,2))
for( jj in 1:p){
  acf(beta_logsigma2_post_sample_1[,jj], main = paste0("Posterior beta",jj))
}



plot( beta_logsigma2_post_sample_1[,1] )
plot( beta_logsigma2_post_sample_1[,2] )





