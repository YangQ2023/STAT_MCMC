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
# data preprocessing steps......
install.packages("dplyr")
set.seed(1)
x1<-matrix(c(rep(1, 100)), ncol=1)
x2_x3 <- matrix(rnorm(200), nrow = 100, ncol = 2)
x <- cbind(x1, x2_x3)
eps<-matrix(rnorm(100),nrow =100, ncol=1)
beta <- matrix(c(1, 2, 0.5), ncol = 1)
y<-x%*%beta+eps


#Gibbs algorithms 
n=5000
sigma2_Gibbs_sample=rep(0, l=n)
beta_Gibbs_sample<- matrix(0, nrow =n, ncol = 3)


#computation steps prepare for prior beta and sigma2
sample_size=100
alpha=0.1
beta_prime=0.1
set.seed(2)
sigma2[1] <- rinvgamma(1,0.1,0.1)
dim(beta)
beta_Gibbs_sample[1,]<-beta
dim(beta_Gibbs_sample[1,])
gamma[1]<-solve((diag(3)+(1/sigma2[1])*t(x)%*%x))%*%t(x)%*%y
omega[1]<-solve((diag(3)+(1/sigma2[1])*t(x)%*%x))
a[1]=sample_size/2+alpha
b[1]=(t(y-x%*%beta_Gibbs_sample[1,])%%(y-x%*%beta_Gibbs_sample[1,]))/2+beta_prime

#Loop structure for gibbs sampling interations  
for(j in 2:5000) {
  x.seq[j] <- rbinom(1, N.seq[j-1], p.seq[j -1])
  p.seq[j] <- rbeta(1, (x.seq[j] + 2),
                    (N.seq[j -1] - x.seq[j] + 4))
  N.seq[j] <- rpois(1, 16 * (1 - p.seq[j])) + x.seq[j]
}


