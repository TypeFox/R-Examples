rm(list=ls())
if(FALSE) {
    files <- list.files("../R/",full.names=TRUE)
    files <- files[-grep("\\~",files)]
    for(ff in files) source(ff)
} else {
    library(regress)
}

library(MASS) ## needed for mvrnorm
n <- 100
mu <- c(1,2)
Sigma <- matrix(c(10,5,5,10),2,2)
Y <- mvrnorm(n,mu,Sigma)

## this simulates multivariate normal rvs
y <- as.vector(t(Y))
X <- kronecker(rep(1,n),diag(1,2))
V1 <- matrix(c(1,0,0,0),2,2)
V2 <- matrix(c(0,0,0,1),2,2)
V3 <- matrix(c(0,1,1,0),2,2)
sig1 <- kronecker(diag(1,n),V1)
sig2 <- kronecker(diag(1,n),V2)
gam <- kronecker(diag(1,n),V3)

reg.obj <- regress(y~X-1,~sig1+sig2+gam,identity=FALSE,start=c(1,1,0.5),verbose=2)
reg.obj

Sig <- var(Y)
round(Sig, 3)

if(sum((reg.obj$sigma - Sig[c(1,4,2)])^2) > 0.001) stop()

n <- 101
x <- (0:100)/100
prb <- (1 + cos(2*pi*x))/2

## indicator which x-values are observed
obs <- (runif(n) < prb)
mu <- 1 + x + sin(3*log(x+0.1))/(2+x)

## full data
y <- mu + rnorm(n,0,0.1)

## 2. Fit the cubic spline to observed data
one <- rep(1,n)
d <- abs(x %*% t(one) - one %*% t(x))
V <- d^3
X <- model.matrix(y~1+x)
fit <- regress(y[obs]~X[obs,],~V[obs,obs],
identity=TRUE)

## 3. BLP at all x given the observed data
wlsfit <- X%*%fit$beta
blp <- wlsfit + fit$sigma[1]*V[,obs] %*% fit$W%*%(y-wlsfit)[obs]

## 4. Plot of results
plot(x[obs],y[obs], xlab="x",ylab="y")
lines(x,mu,lty=2)
lines(x,blp)
