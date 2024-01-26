#coped codes from the lecture_southhampton unvi
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
# data prepare steps......................
install.packages("dplyr")
set.seed(100)
x1<-matrix(c(rep(1, 100)), ncol=1)
x2_x3 <- matrix(rnorm(200), nrow = 100, ncol = 2)
x <- cbind(x1, x2_x3)
eps < -matrix(rnorm(100),nrow =100, ncol=1)
beta <- matrix(c(1, 2, 0.5), ncol = 1)
y <- x%*%beta + eps

#create vector list for Gibbs sampling
burnin= 1000
thin = 50
n = 5000*thin
beta_Gibbs_sample= matrix(0, nrow =n/thin, ncol = 3)
sigma2_Gibbs_sample=rep(0, l=n/thin)


#computation to prepare for prior beta and sigma2
library("stats")
sample_size=100 # sample sizes from data
alpha=0.1 # parameters for prior sigma2's inverse_gamma
beta_prime=0.1 # parameters for prior sigma2's inverse_gamma
set.seed(1)

#gamma[1,]<-solve(diag(3)+1/sigma2_Gibbs_sample[1]*t(x)%*%x)%*%t(x)%*%y
#omega[,,1]<-matrix(solve((diag(3)+1/sigma2_Gibbs_sample[1]*t(x)%*%x)))
#a[1]=sample_size/2+alpha
#b[1]=(t(y-x%*%as.matrix(beta_Gibbs_sample[1,])))%*%(y-x%*%as.matrix(beta_Gibbs_sample[1,]))/2 + beta_prime
#library(mvtnorm)
sigma2_new_gibbs=1
counter=1

#Gibbs algorithms 
for(j in 2:(burnin + n) ) {
  omega <- solve( diag(3) + t(x)%*%x/sigma2_new_gibbs )
  gamma  <- omega %*% t(x) %*% y/sigma2_new_gibbs
  
  #beta's full conditional density 
  beta_new_gibbs <- as.vector(mvrnorm(n = 1,mu = gamma, Sigma = omega))
  a = sample_size/2 + alpha
  e = y - x%*%as.matrix(beta_new_gibbs)
  b = as.numeric( crossprod(e)/2 + beta_prime )
  #sigma2's full conditional density
  sigma2_new_gibbs <- 1/rgamma(1, shape = a, rate = b)
  
  #1/rgamma(1, shape = a, rate = b)
  #rinvgamma(1, shape = a, rate = b)
  
  
  #adding burnin and thin steps in gibbs sampling
  if( j > burnin & (j%%thin) == 0){
    # save new sample
    beta_Gibbs_sample[counter,] = beta_new_gibbs
    sigma2_Gibbs_sample[counter] = sigma2_new_gibbs
    counter = counter + 1
  }
}
#--------------------------------------------------------------------------------
par(mfrow=c(1,3))
hist(beta_Gibbs_sample[,1])
hist(beta_Gibbs_sample[,2])
hist(beta_Gibbs_sample[,3])

acf(sigma2_Gibbs_sample )
mean(sigma2_Gibbs_sample)
hist(sigma2_Gibbs_sample,freq=F)

?rinvgamma


plot(sigma2_Gibbs_sample)
hist(sigma2_Gibbs_sample,freq=F)

plot(gamma[,1] )
plot(gamma[,2] )

sigma2_Gibbs_sample
acf(sigma2_Gibbs_sample)
hist(sigma2_Gibbs_sample,freq=F)
hist(beta_Gibbs_sample[,3],freq=F)
acf(beta_Gibbs_sample[,3])
beta_Gibbs_sample


