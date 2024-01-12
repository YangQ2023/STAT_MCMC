#exercise
i=1:10^4
x[1]=rnorm(1)
?rnorm
r=0.9
for (i in 2:10^4){
  x[i]=r*x[i-1]+rnorm(1) }
hist(x,freq=F,col="wheat2",main="")
curve(dnorm(x,sd=1/sqrt(1-r^2)),add=T,col="tomato")
nBoot=200
B=array(0,dim=c(nBoot, 2))

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

#--------------------------------------------------------------------------------
#1/5/2024
# data preprocessing steps......
install.packages("dplyr")
set.seed(1)
x1<-matrix(c(rep(1, 100)), ncol=1)
x2_x3 <- matrix(rnorm(200), nrow = 100, ncol = 2)
x <- cbind(x1, x2_x3)
eps<-matrix(rnorm(100),nrow =100, ncol=1)
beta <- matrix(c(1, 2, 0.5), ncol = 1)
y<-x%*%beta+eps

#1) directly generate 10,000 beta value from posterior disbribution
sigma2 <- 1
gamma<-solve((diag(3)+(1/sigma2)*t(x)%*%x))%*%t(x)%*%y
omega<-solve((diag(3)+(1/sigma2)*t(x)%*%x))
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
install.packages("mvtnorm")
install.packages("ggplot2")
library("mvtnorm")
mean<-rep(0,3)
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


beta = beta_post_sample
beta_1_MH<-hist(beta[,1],freq=F,xlim=c(0,2))
acf(beta[,1])
mean(beta[,1])
var(beta[,1])
beta_2_MH<-hist(beta[,2],freq=F)
acf(beta[,2])
beta_3_MH<-hist(beta[,3],freq=F,xlim=c(0,1.5))


plot( beta_1, col=rgb(0,0,1,1/4), xlim=c(-4,4))  # first histogram
plot( beta_1_MH, col=rgb(1,0,0,1/4), xlim=c(-4,4), add=TRUE)  


par(mfrow=c(1,3))
hist(sample_beta_1[,1],freq=F, col=rgb(0,0,1,1/4), xlim=c(0,2), ylim=c(0,5) )
hist(beta[,1], freq=F, col=rgb(1,0,0,1/4), xlim=c(0,2), add=TRUE)

hist(sample_beta_1[,2],freq=F, col=rgb(0,0,1,1/4), xlim=c(1.2,2.8), ylim=c(0,5) )
hist(beta[,2], freq=F, col=rgb(1,0,0,1/4), xlim=c(0,2), add=TRUE)

hist(sample_beta_1[,3],freq=F, col=rgb(0,0,1,1/4), xlim=c(0,1), ylim=c(0,5) )
hist(beta[,3], freq=F, col=rgb(1,0,0,1/4), xlim=c(0,2), add=TRUE)



library(invgamma)#install.packages("invgamma")
dev.off()
?dev.off()

xs= seq(0,20,l=1000)

plot( xs, dinvgamma(xs, 0.1, 0.1), type="l" )

#1/11/2024------------------------------
print("hello world")





