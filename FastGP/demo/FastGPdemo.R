#A demonstration of functions from the GPPack
library(FastGP)
library(mvtnorm)
library(MASS)
library(rbenchmark)

N <- 200
sigma <- 1 #variance parameter in the covariance function
phi <- 1 #decay parameter for the exponential kernel
Sig <- as.matrix(sigma*exp(-as.matrix(dist(seq(1,N)))^2/phi)) #test covariance function

#RcppEigen based inverse vs. R inverse
benchmark(solve(Sig),rcppeigen_invert_matrix(Sig))
#Rcpp based Toeplitz inverse vs. R inverse
benchmark(solve(Sig),tinv(Sig))
#RcppEigen based determinant vs. R determinant
benchmark(det(Sig),rcppeigen_get_det(Sig))
#Rcpp based dmvnorm vs. mvtnorm package (with Toeplitz flag)
benchmark(rcpp_log_dmvnorm(S=Sig, mu=rep(0,N), x=rep(1,N),istoep=TRUE),dmvnorm(rep(1,N),mean=rep(0,N),sigma = Sig,log=T))
#Rcpp based dmvnorm vs. mvtnorm package (without Toeplitz flag)
benchmark(rcpp_log_dmvnorm(S=Sig, mu=rep(0,N), x=rep(1,N),istoep=FALSE),dmvnorm(rep(1,N),mean=rep(0,N),sigma = Sig,log=T))
#Rcpp based distance versus R distance
benchmark(dist(as.matrix(seq(1,N))),rcpp_distance(matrix(seq(1,N),nrow=N),N,1))
#Rcpp based rmvnorm versus mvtnorm rmvnorm
benchmark(rcpp_rmvnorm(10,Sig,rep(0,N)),rmvnorm(10, mean = rep(0, N), sigma = Sig))
#Rcpp based rmvnorm versus MASS mvtnorm 
benchmark(rcpp_rmvnorm(10,Sig,rep(0,N)),mvrnorm(10, mu = rep(0, N), Sigma = Sig))

#Now a demo of elliptical slice sampling

#Relevant Parameters:
A <- 1 #amplitude of the sin function used for our signal
T <- 5 #period of the sin function used as our signal
sigma <- 10 #variance parameter in the covariance function
phi <- 100 #decay parameter for the exponential kernel
N <- 100 #number of time points 

#A function that creates a signal, which can be toggled for multiple shapes. In the present iteration we use a sin function with variable period.
signal <- function(t) {
  return(A*sin(t/T))
}

#Define our covariance matrix using the exponential kernel
S <- as.matrix(sigma*exp(-as.matrix(dist(seq(1,N)))^2/phi))
#Generate a copy of the signal on the time scale of 1...N
t <- mvrnorm(n = 1, mu = rep(0,100), Sigma = S , tol = 1e-6)
s <- signal(seq(1,N)+t)+rnorm(N,sd=sqrt(.001))

#log-lik function
log_lik <- function(w,s){
  return (dmvnorm(s, mean = signal(seq(1:N)+w),sigma=.001*diag(1,N),log=T))
}

#rcpp based log-lik function
log_lik_rcpp <- function(w,s){
  return (rcpp_log_dmvnorm(S = .001*diag(1,N),mu=signal(seq(1:N)+w),s,istoep=TRUE))
}
X <- matrix(seq(1,N),nrow=N)
Sig <- sigma*exp(-rcpp_distance(X,dim(X)[1],dim(X)[2])^2/phi)
mcmc_samples <- ess(log_lik_rcpp,s,Sig,1000,250,100,TRUE)
par(mfrow=c(1,1))
plot(colMeans(mcmc_samples), type="l")
points(t, type="l", lwd=3, col="red")
points(colMeans(mcmc_samples) + 2* apply(mcmc_samples, 2, sd), type="l", lty="dashed")
points(colMeans(mcmc_samples) - 2* apply(mcmc_samples, 2, sd), type="l", lty="dashed")

#test ess with rcpp likelihood versus non rcpp likelihood
benchmark(ess(log_lik_rcpp,s,Sig,1000,250,100,FALSE),ess(log_lik,s,Sig,1000,250,100,TRUE),replications=2)