#coped codes from the lecture_southhampton unviersity
gibbs_sampler=function(x0 , p0 , N0 , n, J){
  # Single path method for generating n samples of (X,p,N) value# sampler. In total J+n samples are generated but the first J # (x0,p0,N0) is the initial value.
  x.seq <- p.seq <- N.seq <- rep(NA , J+n+1)
  x.seq[1] <- x0; p.seq[1] <- p0; N.seq[1] <- N0
  for(j in 2:(J+n+1)) {
    x.seq[j] <- rbinom(1, N.seq[j-1], p.seq[j -1])
    p.seq[j] <- rbeta(1, (x.seq[j] + 2),
                      (N.seq[j -1] - x.seq[j] + 4))
    N.seq[j] <- rpois(1, 16 * (1 - p.seq[j])) + x.seq[j]
  }
  result <- list(X=x.seq[(J+2):(J+n+1)] ,p=p.seq[(J+2):(J+n+1)] ,
                 N=N.seq[(J+2):(J+n+1)])
  return(result)}

#run above function
set.seed(2024)# needed?
n=1000
mysamples =gibbs_sampler(8, 0.5, 16, n=n, J=100)
par(mfrow=c(1 ,3))
hist(mysamples$X)
hist(mysamples$p)
hist(mysamples$N)

#1/11/2024_________________________________________________________________________________________________________
#Gibbs sampling for two posterior parameters beta and sigma2......
rm(list=ls())
source("./dataGeneration.R")
#create vector list for Gibbs sampling
burnin= 1000
thin = 50
n = 5000*thin
beta_Gibbs_sample= matrix(0, nrow =n/thin, ncol = 3)
sigma2_Gibbs_sample=rep(0, l=n/thin)

#assign values for parameters 
sample_size=100 # sample sizes from data
alpha=0.1 # shape parameters for prior sigma2's inverse_gamma
beta_prime=0.1 # rate parameters for prior sigma2's inverse_gamma

#assign inital sigma2 values and assign counter=1
sigma2_old_gibbs=1
counter=1

#Gibbs algorithms 
for(j in 1:(burnin + n) ) {
  omega <- solve( diag(3) + t(x)%*%x/sigma2_old_gibbs )
  gamma  <- omega %*% t(x) %*% y/sigma2_old_gibbs
  
  #beta's full conditional density 
  beta_new_gibbs <- as.vector(mvrnorm(n = 1,mu = gamma, Sigma = omega))
  a = sample_size/2 + alpha
  e = y - x%*%as.matrix(beta_new_gibbs)
  b = as.numeric( crossprod(e)/2 + beta_prime )
  
  #sigma2's full conditional density
  sigma2_new_gibbs <- 1/rgamma(1, shape = a, rate = b)
  
  #adding burnin and thin steps in Gibbs sampling
  if( j > burnin & (j%%thin) == 0){
    # save new sample
    beta_Gibbs_sample[counter,] = beta_new_gibbs
    sigma2_Gibbs_sample[counter] = sigma2_new_gibbs
    counter = counter + 1
  }
  
  sigma2_old_gibbs=sigma2_new_gibbs
}
#-----------------------------------------------------------------------------------------
dev.off()
par(mfrow=c(2,2))
hist(beta_Gibbs_sample[, 1], prob = TRUE, 
     main = paste0("Posterior beta1"))
abline(v = 1, col = "red", lty = 1, lwd = 2)


hist(beta_Gibbs_sample[, 2], prob = TRUE, 
     main = paste0("Posterior beta2"))
abline(v = 2, col = "red", lty = 1, lwd = 2)


hist(beta_Gibbs_sample[, 3], prob = TRUE, 
     main = paste0("Posterior beta3"))
abline(v = 0.5, col = "red", lty = 1, lwd = 2)


hist(sigma2_Gibbs_sample, prob=TRUE,main="sigma2")
abline(v = 1, col = "red", lty = 1, lwd = 2)

##################################################################################
par(mfrow=c(2,2))
for( jj in 1:3){
  acf(beta_Gibbs_sample[,jj], main = paste0("Posterior beta",jj))
}
acf((sigma2_Gibbs_sample), main = "sigma2")



