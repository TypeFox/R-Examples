\dontrun{
## Spheroids with lognormal distributed length of axes
## set number of cpu cores (optional)
# options(par.unfoldr=8)

## Intensity: mean number of spheroids per unit volume
lam <- 2500

## simulation parameters
theta <- list("size"=list("meanlog"=-2.5,"sdlog"=0.5),
		      "shape"=list("alpha"=0.5),"orientation"=list("kappa"=1.5))
## simualtion
set.seed(1234)
S <- simSpheroidSystem(theta,lam,size="rlnorm",
			orientation="rbetaiso",box=list(c(0,5)),stype="prolate",pl=101)
## unfolding
sp <- verticalSection(S,2.5)
ret <- unfold(sp,c(15,12,10),kap=1.25)
cat("Intensities: ", sum(ret$N_V)/25, "vs.",lam,"\n")

## uncomment for plot of joint 3d histogram
# require("rgl")
# cols <- c("#0000FF","#00FF00","#FF0000","#FF00FF","#FFFF00","#00FFFF")
# trivarHist(ret$N_V,col=cols,scale=0.9)

param3d <- parameters3d(S)
paramEst <- parameterEstimates(ret$N_V,ret$breaks)

## Marginal histograms
op <- par(mfrow = c(3, 2))
hist(param3d$a[param3d$a<max(ret$breaks$size)],
 main=expression(paste("3D Histogram ", a)),
 breaks=ret$breaks$size,col="gray",right=FALSE,freq=FALSE,xlab="a",ylim=c(0,25))
hist(paramEst$a,
 main=expression(paste("Estimated histogram ",hat(a))),
 breaks=ret$breaks$size,
 right=FALSE,freq=FALSE,col="gray",
 xlab=expression(hat(a)),ylim=c(0,25))
hist(param3d$Theta[param3d$Theta<max(ret$breaks$angle)],
 main=expression(paste("3D Histogram ", theta)),
 breaks=ret$breaks$angle,col="gray",right=FALSE,freq=FALSE,
 xlab=expression(theta),ylim=c(0,2))
hist(paramEst$Theta,
 main=expression(paste("Estimated Histogram ", hat(theta))),
 breaks=ret$breaks$angle,
 right=FALSE,freq=FALSE,col="gray",
 xlab=expression(hat(theta)),ylim=c(0,2))
hist(param3d$s,main=expression(paste("3D Histogram ", s)),
 col="gray",breaks=ret$breaks$shape,
 right=FALSE,freq=FALSE,xlab="s")
hist(paramEst$s,main=expression(paste("Estimated Histogram ", hat(s))),
 breaks=ret$breaks$shape,
right=FALSE,freq=FALSE,col="gray",xlab=expression(hat(s)))
par(op)
}
