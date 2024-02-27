#MH for jointly beta and logsigma2
rm(list=ls())
source("./dataGeneration.R")

#computation to prepare for prior beta and sigma2
#p=4
burnin= 1000
thin = 50
nmc = 5000*thin
n=100 # sample sizes from data
alpha=0.1 # parameters for prior sigma2's inverse_gamma
beta_prime=0.1 # parameters for prior sigma2's inverse_gamma
beta_logsigma2_post_sample_1= matrix(0, nrow =nmc/thin, ncol = 4) # sampling vectors for four dimensions
beta_logsigma2_post_1 <- matrix(c(1, 2, 0.5, 1), ncol = 1) #initial values

#assign initial values for beta and sigma2 before MH_algorithm
beta_logsigma2_old_MH = as.vector(beta_logsigma2_post_1)
counter=1

#define log_target functions to make calculation accurate and efficiency
log_target = function(beta, tau2){
  
  e = y - x%*%as.matrix(beta)
  
  b = crossprod(e)/(2*exp(tau2))
  
  c = crossprod(beta)/2
  
  log_den = (-(n/2)-alpha)*tau2 - b - c - (beta_prime/exp(tau2))
  
  return(log_den)
  
}

#MH algorithm for beta and logsigma2 
for (j in 2:(burnin + nmc) ){
  # proposal from multivariate student t distribution 
  #set up the parameters for proposal multivariate_student t distribution
  Y_2 <- rtmvt(1, mean = beta_logsigma2_old_MH, sigma = diag(4), df = 3,
               lower = rep(-Inf, length = 4), upper = rep(Inf, length = 4))
  Y_2 <- as.vector(Y_2)
  
  # target density--------------------------------------------
  log_den_new = log_target(Y_2[1:3],Y_2[4])
  log_den_old = log_target(beta_logsigma2_old_MH[1:3],beta_logsigma2_old_MH[4])
  
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

plot(beta_logsigma2_post_sample_1)
par(mfrow=c(2,2))
for( jj in 1:3){
  hist( beta_logsigma2_post_sample_1[,jj], prob=TRUE, main = paste0("Posterior beta",jj))
}
hist(exp(beta_logsigma2_post_sample_1[,4]), prob=TRUE,main="sigma2")

par(mfrow=c(2,2))
for( jj in 1:3){
  acf(beta_logsigma2_post_sample_1[,jj], main = paste0("Posterior beta",jj))
}
acf(exp(beta_logsigma2_post_sample_1[,4]), main = "sigma2")

dev.off()
#due to using identical proposal density function lead to kindly high rejection from the sampling. Use Hessian matrix to replace the identity matrix
#in the proposal density function, to solve this high rejection issues.
plot(exp(beta_logsigma2_post_sample_1[,1]), cex=0.5)
plot(exp(beta_logsigma2_post_sample_1[,2]), cex=0.5) 
plot(exp(beta_logsigma2_post_sample_1[,3]), cex=0.5) 
plot(exp(beta_logsigma2_post_sample_1[,4]), cex=0.5) 
