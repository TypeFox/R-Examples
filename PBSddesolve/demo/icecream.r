require(PBSddesolve)
if (!require(PBSmodelling)) stop("The package PBSmodelling must be installed for this demo")

local(env=.PBSddeEnv, expr={

#close any existing windows
closeWin("window")

#store working directory
if (!exists("oldwd") || getwd() != system.file(package = "PBSddesolve")) {
	oldwd <- getwd()
	setwd(system.file(package = "PBSddesolve")) }

# Sine waves used for swithing
sinWave <- function(t,aa,pp) { sin(2*pi*(aa + (t/pp)) ) }

# the switch function
mySwitch <- function(t,y,parms) {
	a <- parms$a; p <- parms$p;
	c( sinWave(t,a[1],p[1]), sinWave(t,a[2],p[2]) ) }

# the map function
myMap <- function(t,y,swID,parms) {
	Y <- parms$Y; f <- parms$f;
	if (swID==1) y <- y + Y else y <- (1-f)*y }
myMapCollect <- function(t,y,swID,parms) {
	Y <- parms$Y; f <- parms$f;
	swTimes = .PBSddeEnv$swTimes
	swTimes[[swID]] <- c( swTimes[[swID]], t ) #collect switch times
	assign("swTimes", swTimes, envir=.PBSddeEnv)
	if (swID==1) y <- y + Y else y <- (1-f)*y }

# the gradient function
myGrad <- function(t,y,parms) {
	-parms$r*y }

# modified by RH (2014-05-14)
#  changed yout1$t to yout1$time, yout2$t to yout2$time, to eliminate partial matching warnings.
icePlot <- function() {
	assign("swTimes", list(c(), c()),envir=.PBSddeEnv)   # store times when map is called
	getWinVal(scope="L");

	# calculate actual time offset for deliveries
	if (rgive) agive <- runif(1,min=0.0001,max=max(0.0001,ogive)) else agive <- ogive
	if (rtake) atake <- runif(1,min=0.0001,max=max(0.0001,otake)) else atake <- otake
	setWinVal(list(agive=round(agive,2),atake=round(atake,2)));
	
	a <- c(agive,atake); p <- c(tgive,ttake); 
	gelati <- list(a=a,p=p,r=r,Y=Y,f=f);
	
	tt1 <- seq(from,to,by); tt2 <- seq(from,to,see);
	# solve the dde system
	yout1 <- dde(y=y0, times=tt1, func=myGrad, switchfunc=mySwitch, mapfunc=myMapCollect, parms=gelati)
	yout2 <- dde(y=y0, times=tt2, func=myGrad, switchfunc=mySwitch, mapfunc=myMap, parms=gelati)
	xadd <- swTimes[[1]][-1]; yadd <- sinWave(xadd,a[1],p[1]);
	xsub <- swTimes[[2]][-1]; ysub <- sinWave(xsub,a[2],p[2]);
	zadd <- is.element(yout1$time,xadd); zsub <- is.element(yout1$time,xsub);

	#plot results
	frame(); resetGraph();
	expandGraph(mfrow=c(2,1),mar=c(3,3,2,1),mgp=c(1.75,.5,0));
	plot( yout1$time, yout1$ice, type="l", lwd=2, xlab="days", ylab="ice cream", 
		main="DDE with two switches", col="cornflowerblue" );
	points(yout2$time,yout2$ice,col="blue",cex=1.5,pch=16);
	points(yout1$time[zadd], yout1$ice[zadd], col="green3",cex=1.5,pch=16)
	points(yout1$time[zsub], yout1$ice[zsub], col="red",cex=1.5,pch=16)

	#display switch functions (only for demo illustration)
	plot(tt1,sinWave(tt1,a[1],p[1]), type="l", col="green3", 
		main="Switch function and corresponding map function call",
		xlab="days", ylab="y")
	lines(tt1,sinWave(tt1,a[2],p[2]), type="l", col="red" )
	points(xadd, yadd, col="green3",cex=1.5,pch=16)
	points(xsub, ysub, col="red",cex=1.5,pch=16)
}
#restore working directory once demo is done
#onClose <- function() { setwd(oldwd); }
# Now handled by `.onClosePBSddeExamples`

createWin("demo_files/icecream_win.txt")
})
