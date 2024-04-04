rm(list=ls())
library(invgamma)
library(mvtnorm)

a=0.1
b=0.1
sigmasq <- 2
n <- 100
p <- 3
I <- diag(n)
y <- rnorm(n)
X <- matrix(rnorm(n*p),n,p)
dim(X)
# y = e where e ~ N(0,1)
# precalculatyion
XXt <- X%*%t(X)
dim(XXt)
V <- sigmasq*I + XXt


## log density of sigmasq | y
-(a+1)*log(sigmasq) - b/sigmasq + -0.5*determinant(V)$mod[1] - sum( y*solve(V,y) )/(2) 

#check density values for equation of log density of sigmasq |y##################
dinvgamma(sigmasq, shape=a, rate=b, log=TRUE)
-(a+1)*log(sigmasq) - b/sigmasq + (a*log(b) - log(gamma(a)) )

-0.5*determinant(V)$mod[1] - sum( y*solve(V,y) )/(2) + (-n/2)*log(2*pi)
dmvnorm(y, mean=rep(0,n), V, log=TRUE)


#################################################################################
Ip <- diag(p)
XtX <- crossprod(X)
#Matrix determinant lemma: using low dimension matrix for fast calculation########
V <- sigmasq*I + XXt
M <- Ip + XtX/sigmasq
# n x n determinant
determinant(V)$mod[1]

# p x p determinant
n*log(sigmasq) + determinant(M)$mod[1] 

# p x p determinant is much faster
microbenchmark::microbenchmark(determinant(V)$mod[1], n*log(sigmasq) + determinant(M)$mod[1])

#replaced the lower dimension matrix to the high dimension matrix for fast calculation#########
-(a+1)*log(sigmasq) - b/sigmasq + -0.5*determinant(V)$mod[1] - sum( y*solve(V,y) )/(2)
-(a+1)*log(sigmasq) - b/sigmasq + -0.5*( n*log(sigmasq) + determinant(M)$mod[1] ) - sum( y*solve(V,y) )/(2)



# variable transformation #######################################################################
sigmasq <- exp(tausq)
tausq <- log(sigmasq)
Mtau <- Ip + XtX/exp(tausq)
Vtau <- sigmasq*I + XXt 
-(a+1)*tausq - b/exp(tausq) + -0.5*determinant(Vtau)$mod[1] - sum( y*solve(Vtau,y) )/(2) + tausq
-a*tausq - b/exp(tausq) + -0.5*determinant(Vtau)$mod[1] - sum( y*solve(Vtau,y) )/(2)
-a*tausq - b/exp(tausq) + -0.5*( n*tausq + determinant(Mtau)$mod[1] ) - sum( y*solve(Vtau,y) )/(2)

################################################################################################
logden_tausq <- function(tausq){

  sigmasq <- exp(tausq)
  Mtau <- Ip + XtX/exp(tausq)

  Xty <- crossprod(X,y)
  quad <-  0.5*( sum(y^2) - sum( Xty * solve( sigmasq*Ip + XtX, Xty) ) )/sigmasq

  logden <- -(a+1)*log(sigmasq) - b/sigmasq + -0.5*( n*log(sigmasq) + determinant(Mtau)$mod[1] ) - quad

  return(logden)

}

tausqs <- seq(-2,4,l=1000)
par(mfrow=c(1,2))
plot( tausqs, Vectorize(logden_tausq)(tausqs), type="l")
plot( tausqs, exp( Vectorize(logden_tausq)(tausqs) ), type="l")
dev.off()

## find the posterior mode: define the potential mean and varaicne of the posterior distribution########
pmode <- optimize(logden_tausq, c(-0.3,0.3), maximum = TRUE )$maximum
pmode <- optimize(logden_tausq, c(-0.5,0.5), maximum = TRUE )$maximum #Q1; why change the area from 0.3 to 0.5 produced different pmode?
pvar <- as.numeric( -1/numDeriv::hessian(logden_tausq, pmode))

dev.off()

################## normal is not a good proposal 
tausqs <- seq(-1.5,1.5,l=1000)
plot(tausqs, exp( Vectorize(logden_tausq)(tausqs) ), type="l")
par(new=TRUE)
plot(tausqs, dnorm(tausqs, pmode, sqrt(pvar) ), type="l", col = 4)


plot(tausqs, exp(Vectorize(logden_tausq)(tausqs) ), type="l")
par(new=TRUE)
plot(tausqs, dnorm(tausqs, pmode, sqrt(pvar) ), type="l", col = 4)


logden_tausq(pmode)  - dnorm(pmode, pmode, sqrt(pvar), log=TRUE) 

logC <-  58 

(logC+logden_tausq(pmode))  < (dnorm(pmode, pmode, sqrt(pvar), log=TRUE))


(logC+logden_tausq(pmode)) / dnorm(pmode, pmode, sqrt(pvar), log=TRUE)  



plot( tausqs,  exp(logC + Vectorize(logden_tausq)(tausqs)), type="l")
par(new=TRUE)
lines(tausqs, dnorm(tausqs, pmode, sqrt(pvar) ), type="l", col = 4)


plot(tausqs,  exp(logC + Vectorize(logden_tausq)(tausqs)), type="l", ylim=c(0,1e-10), xlim=c(1,1.5))
lines(tausqs, dnorm(tausqs, pmode, sqrt(pvar) ), type="l", col = 4, ylim=c(0,1e-10), xlim=c(1,1.5))

plot(tausqs,  exp(logC + Vectorize(logden_tausq)(tausqs)), type="l", ylim=c(0,1e-10), xlim=c(1,1.5))
lines(tausqs, 10*dnorm(tausqs, pmode, sqrt(pvar)), type="l", col = 4, ylim=c(0,1e-10), xlim=c(1,1.5))

dev.off()
#################change the proposal from normal to student t###############################################################
library(tmvtnorm)

tausqs <- seq(-1.5,1.5,l=1000)
plot(tausqs, exp( Vectorize(logden_tausq)(tausqs) ), type="l")
par(new=TRUE)
plot( tausqs, dt( (tausqs - pmode)/sqrt(pvar), df=5  ), type="l", col = 4)

logC <- 57
  
plot( tausqs, exp( logC  + Vectorize(logden_tausq)(tausqs) ), type="l")
lines( tausqs, 1.5*dt( (tausqs - pmode)/sqrt(pvar), df=5  ), type="l", col = 4) #Q2 Is "1.5" the constant "M" in rejection sampling? Why we selected 1.5 specifically??


plot( tausqs,  exp( logC + Vectorize(logden_tausq)(tausqs) ), type="l", ylim=c(0,1e-3), xlim=c(1,1.5) )
lines( tausqs, 1.5*dt( (tausqs - pmode)/sqrt(pvar), df=5  ) , type="l", col = 4, ylim=c(0,1e-3), xlim=c(1,1.5))

1/1.5 #Q3 why we need this line ???

#######################################################
log_prop_den <- function(x){
  
  val <- dt((x - pmode)/sqrt(pvar), df=5, log=TRUE  )
  
  return(val)
}


sample_prop <- function(){
  
  dt((x - pmode)/sqrt(2*pvar), df=5, log=TRUE  )
  
}


pmode + rt(1, df=5 )*sqrt(pvar) #Q4 confused on this line??

# E( pmode + rt( 1, df=5 )*sqrt(pvar) ) = pmode
# Var( pmode + rt( 1, df=5 )*sqrt(pvar) ) =  Var( rt( 1, df=5 )*sqrt(pvar) )
# = pvar * Var( rt( 1, df=5 ) )

###############################################################################
tausq_sample <- rep(0,1e5)
beta_sample<- matrix(0, nrow = 1e5, ncol = 3)

# rs for sigmasq part#####################################################################
counter <- 0
while(counter < 1e5){
  
  proposal <- pmode + rt(1, df=5 )*sqrt(pvar)
  logratio <- logden_tausq(proposal) + logC - (log_prop_den(proposal) + log(1.5)) #Q5 why is "+" log(1.5) instead of "-" ??
  logu <- log(runif(1))

  if( logu < logratio){
    counter <- counter + 1
    tausq_sample[counter] <- proposal
    
    # Computation for beta part based on proposal
    gamma_rs <- (solve(diag(3) + (1 / exp(proposal)) * t(X)%*%X) %*% t(X)%*%y) / (exp(proposal))
    omega_rs <- solve(diag(3) + (1 / exp(proposal)) * t(X)%*%X)
    beta_sample[counter] <- mvrnorm(1, mu = gamma_rs, Sigma = omega_rs)
  }
}

##################################################################################################
plot( tausq_sample, cex=0.01 )
hist( tausq_sample, 100, prob=TRUE)

# un-normalize density of tau2
den_tausq <- function(x){
  val <- exp( logden_tausq(x) )
}

inte <- integrate(Vectorize(den_tausq), lower=-1, upper=1)
# normalizing constant
nc <- 1/inte$value

hist( tausq_sample, 100, prob=TRUE)
lines( tausqs,  nc*exp( Vectorize(logden_tausq)(tausqs) ), type="l", col=2, lwd=2 )


###############################################################################

hist(exp(tausq_sample), 100, prob=TRUE)






