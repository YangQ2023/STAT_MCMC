#2)MH sampler for beta values-----------1/11/2024
rm(list=ls())
source("./dataGeneration.R")
burnin = 500
thin = 50
n = 5000*thin
beta_post_sample <- matrix(0, nrow = n/thin, ncol = 3)
beta_old = as.vector( gamma )
counter = 1

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

par(mfrow=c(1,3))
for( jj in 1:3){
  hist( beta_post_sample [,jj], prob=TRUE, main = paste0("Posterior beta",jj))
}

par(mfrow=c(1,3))
for( jj in 1:3){
  acf( beta_post_sample [,jj], prob=TRUE, main = paste0("Posterior beta",jj))
}
