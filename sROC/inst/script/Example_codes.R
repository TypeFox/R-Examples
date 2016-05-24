##------------------------
## Example codes for the package "sROC"
library(sROC)

 
############### 
set.seed(100)
n <- 200
x <- c(rnorm(n/2, mean=-2, sd=1), rnorm(n/2, mean=3, sd=0.8))
x.CDF <- kCDF(x)
x.CDF
summary(CI.CDF(x.CDF))
plot(x.CDF, alpha=0.05, main="Kernel estimate of distribution function")
curve(pnorm(x, mean=-2, sd=1)/2 + pnorm(x, mean=3, sd=0.8)/2, from =-6, to=6, add=TRUE, lty=2, col="blue")

#####################
set.seed(100)
n <- 200
x <- rgamma(n,2,1)
y <- rnorm(n)

xy.ROC <- kROC(x,y, bw.x="pi_sj",bw.y="pi_sj")
plot(xy.ROC)

#####################

set.seed(100)
n <- 200
x <- rlnorm(n, mean=2, sd=1)
y <- rnorm(n,mean=2,sd=2)

xy.ROC <- kROC(c(x,NA,NA),c(y,1.2, NA), na.rm=TRUE)
plot(xy.ROC)
AUC(xy.ROC)