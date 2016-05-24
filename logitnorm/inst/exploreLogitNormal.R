# TODO: Add comment
# 
# Author: twutz
###############################################################################

library(ggplot2)

xGrid = c(0+.Machine$double.eps,seq(0,1, length.out=81)[-c(1,81)],1-.Machine$double.eps)

#explore chanign sigma at mu=0
theta0 <- expand.grid(mu=seq(0,2,length.out=9), sigma=10^seq(-0.5,0.5,length.out=5))
n <- nrow(theta0)


.calcDensityGrid <- function(
	### Calculate lognormal desnity for given combinations
	theta0	##<< matrix with columns mu and sigma
	,xGrid = seq(0,1, length.out=81)[-c(1,81)]
){ 
	dx <- apply( theta0, 1, function(theta0i){
			dx <- dlogitnorm(xGrid, mu=theta0i[1], sigma=theta0i[2])
		})
	dimnames(dx) <- list(iX=NULL,iTheta=NULL)
	ds <- melt(dx)
	ds[1:10,]
	ds$mu <- rep(as.factor(round(theta0[,1],2)), each=length(xGrid))
	ds$sigma <- rep(as.factor(round(theta0[,2],2)), each=length(xGrid))
	ds$x <- rep(xGrid, nrow(theta0))
	ds
	### data frame with columns value,x,mu and sigma 
}

.calcCdfGrid <- function(
	### Calculate lognormal desnity for given combinations
	theta0	##<< matrix with columns mu and sigma
	,xGrid = seq(0,1, length.out=81)[-c(1,81)]
){ 
	dx <- apply( theta0, 1, function(theta0i){
			dx <- plogitnorm(xGrid, mu=theta0i[1], sigma=theta0i[2])
		})
	dimnames(dx) <- list(iX=NULL,iTheta=NULL)
	ds <- melt(dx)
	ds[1:10,]
	ds$mu <- rep(as.factor(round(theta0[,1],2)), each=length(xGrid))
	ds$sigma <- rep(as.factor(round(theta0[,2],2)), each=length(xGrid))
	ds$x <- rep(xGrid, nrow(theta0))
	ds
	### data frame with columns value,x,mu and sigma 
}




#qplot(xGrid,value,data=ds[ds$mu==0,], geom="line", color=sigma )
windows(width=7,height=4.5)

ds <- .calcDensityGrid(theta0,xGrid=xGrid)
qplot(xGrid,value,data=ds, geom="line", color=sigma, ylab="density")+ facet_wrap(~mu,scales="free")+opts(axis.title.x = theme_blank())
#savePlot("logitnormDensityGrid",type="pdf")

windows(width=4,height=4)
ds <- .calcDensityGrid(theta0,xGrid=xGrid)
qplot(xGrid,value,data=ds[ds$mu %in% c(0,1),], geom="line", color=sigma, ylab="density")+ facet_grid(mu~.,scales="free")+opts(axis.title.x = theme_blank())
#savePlot("logitnormDensityGrid2",type="pdf")

ds <- .calcCdfGrid(theta0,xGrid=xGrid)
qplot(xGrid,value,data=ds[ds$mu %in% c(0,1),], geom="line", color=sigma, ylab="cumulative density")+ facet_grid(mu~.,scales="free")+opts(axis.title.x = theme_blank())
#savePlot("logitnormCdfGrid2",type="pdf")



mle=0.8
mu <- seq(0,logit(mle),length.out=10)[-10]
tmp <- twSigmaLogitnorm(mle,mu)
sigma2 <- (logit(mle)-mu)/(2*mle-1)   
theta0 <- cbind(mu,sigma=as.numeric(sqrt(sigma2)))

ds <- .calcDensityGrid(theta0,xGrid=xGrid)
qplot(xGrid,value,data=ds, geom="line", color=mu )
qlogitnorm(p=0.975,mu=theta0[,1],sigma=theta0[,2])
qlogitnorm(p=0.99,mu=theta0[,1],sigma=theta0[,2])

mle=0.8
#mtrace(.ofLogitnormMLE)
plot( .ofLogitnormMLE(mu,mle=mle,quant=0.95,perc=0.999) ~ mu) 
#mtrace(twCoefLogitnormMLE)
theta <- twCoefLogitnormMLE( mle=mle,quant=0.95,perc=0.99)
plot( dlogitnorm(xGrid,mu=theta[1],sigma=theta[2])~xGrid, type="l")
c( (logit(mle)-theta[1])/(2*mle-1), theta[2]^2 )
abline(v=c(mle,0.99),col="gray")

mle=0.5
#mtrace(twCoefLogitnormMLE)
theta <- twCoefLogitnormMLE( mle=mle,quant=0.99,perc=0.999)
q2 <- qlogitnorm(0.999,mu=theta[1],sigma=theta[2])
plot( dlogitnorm(xGrid,mu=theta[1],sigma=theta[2])~xGrid, type="l")
abline(v=c(mle,q2),col="gray")

mle=0.1
#mtrace(twCoefLogitnormMLE)
theta <- twCoefLogitnormMLE( mle=mle,quant=0.5,perc=0.99)
quant2 <- qlogitnorm(0.99,mu=theta[1],sigma=theta[2])
plot( dlogitnorm(xGrid,mu=theta[1],sigma=theta[2])~xGrid, type="l")
abline(v=c(mle,quant2),col="gray")

mle=0.9
#mtrace(twCoefLogitnormMLE)
theta <- twCoefLogitnormMLE( mle=mle,quant=0.98,perc=0.99)
#theta <- twCoefLogitnormMLE( mle=mle,quant=0.99,perc=0.999)
quant2 <- qlogitnorm(0.99,mu=theta[1],sigma=theta[2])
plot( dlogitnorm(xGrid,mu=theta[1],sigma=theta[2])~xGrid, type="l")
abline(v=c(mle,quant2),col="gray")



#------- sigma for the flattest case
mle=0.9
mu=seq(0,if(mle<0.5) mle else 1-mle,length.out=10)[-10]
mu=seq(0,0.25,length.out=10)[-10]
tmp <- twSigmaLogitnorm(mle,mu)
sigma2 <- (logit(mle)-mu)/(2*mle-1)   
theta0 <- cbind(mu,sigma=as.numeric(sqrt(sigma2)))
ds <- .calcDensityGrid(theta0,xGrid=xGrid)
qplot(xGrid,value,data=ds, geom="line", color=mu, ylim=c(0.9,1.2) )

