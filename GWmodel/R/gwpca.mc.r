# Monte-Carlo test of variability - randomises locations and computes standard deviation of the FIRST EIGENVALUE.
# Note - could make this more general to all other eigenvalues...
# Arguments as before, but nsims controls number of simulations

montecarlo.gwpca.1<- function(data, bw, vars, k = 2, nsims=99,robust = FALSE, kernel = "bisquare",
                  adaptive = FALSE,  p = 2, theta = 0, longlat = F,
                  dMat)
{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
  }
  else if (is(data, "data.frame")&&(!missing(dMat)))
     data<-data
  else
     stop("Given data must be a Spatial*DataFrame or data.frame object")
  if (missing(dMat))
  {
    dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
  }
  if(missing(bw)||bw<=0)
    stop("Bandwidth is not specified incorrectly")
  actual <- sd(gwpca(data=data, vars=vars, k = k, robust = robust, kernel = kernel,
                  adaptive = adaptive, bw=bw, dMat=dMat)$var[,1])
  n<-nrow(dp.locat)
	res <- numeric(nsims)
	for (i in 1:nsims) {
    mcs <- sample(n)
		dMat[mcs,]<-dMat[1:n,]
		dMat[,mcs]<-dMat[,1:n]
		res[i] <- sd(gwpca(data=data, vars=vars, k = k, robust = robust, kernel = kernel,
                  adaptive = adaptive, bw=bw, dMat=dMat)$var[,1])}
	simres <- list(actual=actual,sims=res)
	class(simres) <- "mcsims"
	simres}


# As above but instead of simulating using a given bandwidth this one runs the automatic choice each simulation
# More computationally demanding, but more realistic, as it allows for the fact that using cross-validation
# 'mines' for the bandwidth for a given sample data set.
montecarlo.gwpca.2 <- function(data, vars, k = 2, nsims=99,robust = FALSE, kernel = "bisquare",
                  adaptive = FALSE,  p = 2, theta = 0, longlat = F,
                  dMat)
{
  if (is(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
  }
  else if (is(data, "data.frame")&&(!missing(dMat)))
     data<-data
  else
     stop("Given data must be a Spatial*DataFrame or data.frame object")
  if (missing(dMat))
  {
    dMat <- gw.dist(dp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
  }
  bw <- bw.gwpca(data=data, vars=vars, k = k, robust = robust, kernel = kernel,
                  adaptive = adaptive, dMat=dMat)
  actual <- sd(gwpca(data=data, vars=vars, k = k, robust = robust, kernel = kernel,
                  adaptive = adaptive, bw=bw, dMat=dMat)$var[,1])
  n<-nrow(dp.locat)
	res <- numeric(nsims)
	for (i in 1:nsims) {
		mcs <- sample(n)
		dMat[mcs,]<-dMat[1:n,]
		dMat[,mcs]<-dMat[,1:n]
		bw <- bw.gwpca(data=data, vars=vars, k = k, robust = robust, kernel = kernel,
                  adaptive = adaptive, dMat=dMat)
		res[i] <- sd(gwpca(data=data, vars=vars, k = k, robust = robust, kernel = kernel,
                  adaptive = adaptive, bw=bw, dMat=dMat)$var[,1])}
	simres <- list(actual=actual,sims=res)
	class(simres) <- "mcsims"
	simres}


#  The above two methods produce 'mcsims' objects - the following two are print and plot methods for these
#  'sname' is the name of the test statistic
print.mcsims <- function(x,...) {
	cat("Result of Monte-Carlo ",length(x$sims),"simulations:\n")
	cat("Computed p=value is ",sum(x$actual < x$sims)/length(x$sims),"\n")}

plot.mcsims <- function(x,sname="SD of local eigenvalues from randomisations",...) {
	dist.info <- hist(x$sims,plot=FALSE,...)
	plot(c(dist.info$mids,x$actual),c(dist.info$density,0),type='n',xlab=sname,ylab="Density",
  main="Test statistic for eigenvalue nonstationarity")
	lines(dist.info$mids,dist.info$density,type='b')
	abline(v=x$actual)
	text(x$actual,max(dist.info$density)/2,
		paste("Observed SD of local eigenvalues: Estimated p=",
		sprintf("%6.3f\n",sum(x$actual < x$sims)/length(x$sims)),sep=''),srt=90,cex=0.7)}