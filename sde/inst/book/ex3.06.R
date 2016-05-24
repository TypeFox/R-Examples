# ex3.06.R
require(sde)

d <- function(t,x,theta) theta[1]-theta[2] * x
dx <- function(t,x,theta) -theta[2]
dxx <- function(t,x,theta) 0
dt <- function(t,x,theta) 0
s <- function(t,x,theta)  theta[3]*sqrt(x)
sx <- function(t,x,theta) theta[3]/(2*sqrt(x))
sxx <- function(t,x,theta) -theta[3]/(4*x^1.5)

d2 <- function(t,x,theta){
 (4*theta[1]-theta[3]^2)/(2*x*theta[3]^2) -theta[2]*x/2 }
d2x <- function(t,x,theta){
 -theta[2]/2 - (4*theta[1]-theta[3]^2)/(2*x^2*theta[3]^2)}
d2xx <- function(t,x,theta) (4*theta[1]-theta[3]^2)/(x^3*theta[3]^2) 
d2t <- function(t,x,theta) 0
s2 <- function(t,x,theta)  1
s2x <- function(t,x,theta) 0
s2xx <- function(t,x,theta) 0

Euler.LIK <- function(theta) {
   sum(dcEuler(X[2:n], t[2:n], X[1:(n-1)], t[1:(n-1)],
     c(0.5, theta, sqrt(0.05)), d,s, TRUE),na.rm=TRUE)
}

Elerian.LIK <- function(theta) {
   sum(dcElerian(W[2:n], t[2:n], W[1:(n-1)], t[1:(n-1)], 
    c(0.5, theta, sqrt(0.05)), d2, s2, s2x, TRUE) 
	 -0.5*log(X[2:n]*0.05),na.rm=TRUE)
}

Ozaki.LIK <- function(theta) {
   sum(dcOzaki(W[2:n], t[2:n], W[1:(n-1)], t[1:(n-1)], 
     c(0.5, theta, sqrt(0.05)), d2, d2x, s2, TRUE)
	  -0.5*log(X[2:n]*0.05),na.rm=TRUE)
}

Shoji.LIK <- function(theta) {
   sum(dcShoji(W[2:n], t[2:n], W[1:(n-1)], t[1:(n-1)], 
    c(0.5, theta, sqrt(0.05)), d2, d2x,d2xx, d2t, s2, TRUE)
	 -0.5*log(X[2:n]*0.05),na.rm=TRUE)
}

True.LIK <- function(theta) {
   sum(dcCIR(X[2:n], deltat(X), X[1:(n-1)], 
    c(0.5, theta, sqrt(0.05)), TRUE),na.rm=TRUE)
}

pTrue <- function(x) True.LIK(x)
pEuler <- function(x) Euler.LIK(x)
pElerian <- function(x) Elerian.LIK(x)
pOzaki <- function(x) Ozaki.LIK(x)
pShoji <- function(x) Shoji.LIK(x)


# ex3.06.R (cont)
set.seed(123)
X1 <- sde.sim(model="CIR", theta=c(0.5, 0.2, sqrt(0.05)), 
 X0=2,delta=.001, N=500000)
xx <- seq(0.001,0.4, length=50)

par(mfrow=c(2,2))
est <- NULL
for(dt in c(4, 2,1,.5)){
 X <- window(X1, deltat=dt)
 W <- 2*sqrt(X)/sqrt(0.05)
 t <- as.numeric(time(X))
 n <- length(X)
 cat(sprintf("number of observations: %d, Delta=%3.2f\n",n,dt))
 dEuler <- sapply(xx, pEuler)
 dTrue <- sapply(xx, pTrue)
 dElerian <- sapply(xx, pElerian)
 dOzaki <- sapply(xx, pOzaki)
 dShoji <- sapply(xx, pShoji)

 mx <- max(c(dTrue,dEuler,dShoji,dOzaki,dElerian,na.rm=TRUE))
 mn <- min(c(dTrue,dEuler,dShoji,dOzaki,dElerian,na.rm=TRUE))

 matplot(xx,cbind(dTrue,dEuler,dOzaki,dShoji,dElerian),type="l",
  ylim=c(mn,mx),xlab="",ylab="approx", 
  main=sprintf("N=%d, Delta=%3.2f",n,dt),lty=1:5,col=1:5)
 legend(.1,0.7*(mx+mn), lty=1:5,col=1:5,
  legend=c("True", "Euler", "Ozaki", "Shoji","Elerian"))

 tmp <- c(n, dt, optimize(pTrue, c(0.01,.4),maximum=T)$max,
  optimize(pEuler, c(0,.4),maximum=T)$max,
  optimize(pElerian, c(0,1),maximum=T)$max,
  optimize(pOzaki, c(0,.4),maximum=T)$max,
  optimize(pShoji, c(0,.4),maximum=T)$max)
 est <- rbind(est, tmp)
}
dimnames(est)[[2]] <- c("N","Delta","True","Euler", 
 "Elerian", "Ozaki", "Shoji")
dimnames(est)[[1]] <- 1:4
print(est)
par(mfrow=c(1,1))

