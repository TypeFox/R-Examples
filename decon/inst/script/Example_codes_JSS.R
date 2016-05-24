## Example codes for the paper to Journal of Statisitcal Software
## Xiao-Feng Wang and Bin Wang
## the package "decon" can be installed from terminal
## using the command "R CMD INSTALL decon_1.0-5.tar.gz &"
## All codes have been tested at R version 2.10.1 under MAC OS X version 10.6.2!!!

library(decon)

#####################
## the R function to estimate the smooth distribution function
SDF <- function (x, bw = bw.nrd0(x), n = 512, lim=1){
		dx <- lim*sd(x)/20 
		xgrid <- seq(min(x)-dx, max(x)+dx, length = n)
		Fhat <- sapply(x, function(x) pnorm((xgrid-x)/bw)) 
		return(list(x = xgrid, y = rowMeans(Fhat)))
	}


#####################
## Simultion studies
## Deconvolution: the case of homoscedastic errors
## Case 1.1: homoscedastic Laplacian errors
n1 <- 500
x1 <- rnorm(n1, sd=1)
sig1 <- .5
u1 <- ifelse(runif(n1) > 0.5, 1, -1) * rexp(n1,rate=1/sig1)
w1 <- x1+u1
## the rule-of-thumb method may not be accurate, you may try the bootstrap method
bw1 <- bw.dnrd(w1,sig=sig1, error="laplacian")
(f1 <-  DeconPdf(w1,sig1,error='laplacian',bw=bw1, fft=TRUE))

## Case 1.2: homoscedastic normal errors
n2 <- 1000
x2 <- c(rnorm(n2/2,-3,1),rnorm(n2/2,3,1))
sig2 <- .8
u2 <- rnorm(n2, sd=sig2)
w2 <- x2+u2
# estimate the bandwidth with the bootstrap method with resampling
bw2 <- bw.dboot2(w2,sig=sig2, error="normal")
# estimate the unknown density with measurement error
(f2 <-  DeconPdf(w2,sig2,error='normal',bw=bw2, fft=TRUE))
# estimate the distribution function with measurement error
F2 <-  DeconCdf(w2,sig2,error='normal',bw=bw2)

# plot the results
#postscript("pdf_1.eps", pointsize = 16)
par(mfrow=c(2,2))
plot(f1,  col="red", lwd=3, lty=2, xlab="x", ylab="f(x)", main="")
lines(density(x1, from=min(w1), to=max(w1)), lwd=3, lty=1)
lines(density(w1), col="blue", lwd=3, lty=3)
plot(f2,  col="red", lwd=3, lty=2, xlab="x", ylab="f(x)", main="")
lines(density(x2, from=min(w2), to=max(w2)), lwd=3, lty=1)
lines(density(w2), col="blue", lwd=3, lty=3)

plot(F2,  col="red", lwd=3, lty=2, xlab="x", ylab="F(x)", main="")
lines(SDF(x2), lwd=3, lty=1)
lines(SDF(w2), col="blue", lwd=3, lty=3)
# dev.off()


## Deconvolution: the case of heteroscedastic errors
## Case 2: heteroscedastic normal errors
n3 <- 2000
x3 <- rchisq(n3, df=1.5, ncp=0)
sig3 <- 0.7+ x3/max(x3)
u3 <- sapply(sig3, function(x) rnorm(1, sd=x))
w3 <- x3+u3
# estimate the bandwidth using the bootstrap method withou resampling
bw3 <- bw.dboot1(w3,sig=sig3, error="normal")
# estimate the unknown density with measurement error
(f3 <-  DeconPdf(w3,sig3,error="normal", bw=bw3, fft=TRUE))
# estimate the distribution function with measurement error
(F3 <-  DeconCdf(w3,sig3,error="normal", bw=bw3))

# plot the results
#postscript("pdf_cdf_2.eps", pointsize = 16)
par(mfrow=c(1,2))
plot(f3,  col="red", lwd=3, lty=2, ylim=c(0,0.4), xlab="x", ylab="f(x)", main="")
lines(density(x3, adjust=2), lwd=3, lty=1)
lines(density(w3, adjust=2), col="blue", lwd=3, lty=3)

plot(F3,  col="red", lwd=3, lty=2, xlab="x",ylab="F(x)", main="")
lines(SDF(x3), lwd=3, lty=1)
lines(SDF(w3), col="blue", lwd=3, lty=3)
#dev.off()


## Nonparametric regression with errors-in-variables
## Case 3: homoscedastic normal errors
n <- 2000
x <- c(rnorm(n/2,2,1), rnorm(n/2, -2,1))
sig <- .8
u <- sig*rnorm(n)
w <- x+u
e <- rnorm(n, sd=0.2)
y <- x^2-2*x+e
bw1 <- bw.dboot1(w, sig)
# estimate the unknown density with measurement error
(m1 <-  DeconNpr(w, sig, y ,error="normal", from=0.9*min(x), to=0.9*max(x)))
# plot the results
#postscript("npr_1.eps", pointsize = 16)
plot(m1,  col="red", lwd=3, lty=2, xlab="x", ylab="m(x)", main="", zero.line=FALSE)
lines(ksmooth(x,y, kernel = "normal", 2, range.x=c(0.9*min(x),0.9*max(x))), lwd=3, lty=1)
lines(ksmooth(w,y, kernel = "normal", 2, range.x=c(0.9*min(x),0.9*max(x))), col="blue", lwd=3, lty=3)
#dev.off()


#####################
## Real data anlaysis
## Case 1: the study of framingham data
data(framingham)

## take the average of the two measurements at each examination
SBP1 <- (framingham$SBP11+framingham$SBP12)/2
SBP2 <- (framingham$SBP21+framingham$SBP22)/2
mean(SBP2)
var(SBP2)
# measurement errors look quite normal
# So, we can estimate the measurement error variance assuming a normal error model
par(mfrow=c(1,1))
hist(SBP1-SBP2, probability=TRUE, ylim=c(0,0.035))
lines(density(SBP1-SBP2))
## Note for normal variables: if X1|X, X2|X ~ N(X, s^2), X1-X2 ~ N(0, 2*s^2).
## Estimate the measurement error variance
(sig2 <- 0.5*var(SBP1-SBP2))
sig <- sqrt(0.5*var(SBP1-SBP2))


## We consider an iterative bootstrap method with resampling to estimate the bandwidth
(bw0= bw.dnrd(SBP2, sig=sig, error="normal"))
temp=bw0
ibw=rep(0,20)
for (i in 1:20){
	temp= bw.dboot2(SBP2, sig=sig, h0=temp, error="normal", B=1000)
	ibw[i]=temp
}
(ibw1 <- mean(ibw))
## Result is 3.297 (It may be slightly different in each estimation)

SBP2.dec <- DeconPdf(SBP2, sig=sig, error="normal", bw=ibw1,  fft=TRUE)

#postscript("framingham.eps", pointsize = 16)
par(mfrow=c(1,2))
hist(SBP1-SBP2, main="")
lines(density(SBP1-SBP2, adjust=1.2), lwd=3)
plot(SBP2.dec, lwd=3, main="", xlab="SBP2")
lines(density(SBP2, adjust=1.6), lty=2, lwd=3, col="blue")
#dev.off()
## When esitmating the desnity of SBP2 ignoring measurement errors, the bandwidth was selected by eye to be as small as possible while retaining smoothness. So, adjust=1.6 here.


## Case 2: the study of galaxy data
data(galaxy)
## The two bandwidths were chosen by eye to be as small as possible while retaining smoothness.
(m1 <-  DeconNpr(galaxy$V, galaxy$Err, galaxy$Rkpc ,error="normal", bw=9.3))
#postscript("galaxy.eps", pointsize = 16)
plot(m1, xlim=c(0,250), ylim=c(0,15),lwd=3, main="", xlab="x")
lines(ksmooth(galaxy$V, galaxy$Rkpc, kernel="n", 61.8), lwd=3, lty=2, col=4)
#dev.off()




