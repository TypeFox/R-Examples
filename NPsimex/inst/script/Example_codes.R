##------------------------
## Example codes for the package "NPsimex"
library(NPsimex)

plot.simex.density <- function(X.simex,X,...){
	plot(X.simex$x, X.simex$y, type="l", xlab="x", ylab="density", lwd=3, lty=2, col=2,...)
	lines(density(X, bw="SJ"), lwd=3)
	}

 
############### Homoscadestic error
N <- 1000
set.seed(123); X <- c(rnorm(N/2, mean=-2), rnorm(N/2,mean=2)); U <- rnorm(N,sd=1)
msigma <- 0.5
W <- X + msigma*U

par(mfrow=c(2,2))
X.simex1 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=1)
plot.simex.density(X.simex1, X, ylim=c(0,0.25))

X.simex2 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=3)
plot.simex.density(X.simex2, X, ylim=c(0,0.25))

X.simex3 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=8)
plot.simex.density(X.simex3, X, ylim=c(0,0.25))

X.simex4 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=35)
plot.simex.density(X.simex4, X, ylim=c(0,0.25))

#---- Select lambda span
par(mfrow=c(1,2))
spans <- span.select(W, msigma)
plot(spans$span, spans$ISE, type="o", xlab="span", ylab="ISE") 

X.simex <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=spans$span[order(spans$ISE)[1]])
plot.simex.density(X.simex, X, ylim=c(0,0.25))



############### Heteroscadestic error
N <- 1000
set.seed(123); X <- c(rnorm(N/2, mean=-2), rnorm(N/2,mean=2)); U <- rnorm(N,sd=1); msigma <- runif(N,min=0.3,max=0.5)
W <- X + msigma*U

par(mfrow=c(2,2))
X.simex1 <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=1)
plot.simex.density(X.simex1, X, ylim=c(0,0.25))

X.simex2 <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=3)
plot.simex.density(X.simex2, X, ylim=c(0,0.25))

X.simex3 <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=8)
plot.simex.density(X.simex3, X, ylim=c(0,0.25))

X.simex4 <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=35)
plot.simex.density(X.simex4, X, ylim=c(0,0.25))

#---- Select the optimal lambda span
par(mfrow=c(1,2))
spans <- span.H.select(W, msigma, span=c(2,4,6,8,10,12,16,25), approx=TRUE)
plot(spans$span, spans$ISE, type="o", xlab="span", ylab="ISE") 

X.simex <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=spans$span[order(spans$ISE)[1]])
plot.simex.density(X.simex, X,ylim=c(0,0.25))

## more computational time needed if approx=FALSE.

# spans <- span.H.select(W, msigma, span=c(2,4,6,8,10,12,16,25), approx=FALSE)
# plot(spans$span, spans$ISE, type="o", xlab="span", ylab="ISE") 
# 
# X.simex <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=spans$span[order(spans$ISE)[1]])
# plot.simex.density(X.simex, X,ylim=c(0,0.25))





