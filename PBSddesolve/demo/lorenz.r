require(PBSddesolve)
if (!require(PBSmodelling)) stop("The package PBSmodelling must be installed for this demo")

local(env=.PBSddeEnv, expr={

#close any existing windows
closeWin("window")

#store working directory
if (!exists("oldwd") || getwd() != system.file(package = "PBSddesolve")) {
	oldwd <- getwd()
	setwd(system.file(package = "PBSddesolve"))
}

runPlot <- function()
{
	#extract variables from GUI
	getWinVal(scope="L")
	
	#function to calculate gradient at a given time
	myGrad <- function(t, y, parms=NULL) {
		y1 <- sigma*(y[2]-y[1]) 
		y2 <- y[1]*(tau-y[3]) - y[2]
		y3 <- y[1]*y[2] - rho*y[3]
		if (derivative=="yes") {
			deriv <- c(dy1=y1,dy2=y2,dy3=y3) #R will merge previous name found in y
			names( deriv ) <- c( "dy1", "dy2", "dy3" ) #so they must be reset here (LAME)
			return(list(c(y1,y2,y3), deriv))
		}
		else
			return(list(c(y1,y2,y3), NULL))
	}
	#initial values
	yinit <- c(y1=y1,y2=y2,y3=y3)
	#solve ODE
	#if (solver=="PBSddesolve") {
		x <- dde(y=yinit, func=myGrad, times=seq(t0,t1,timestep), hbsize=0)
	#}
	#else if (solver=="deSolve") {
	#	require(deSolve) || stop("deSolve is required")
	#	x <- lsoda(y=yinit, times=seq(from=t0, to=t1, by=timestep), func=myGrad, parms=NULL, rtol=1e-6, atol=1e-4)
	#	x <- as.data.frame(x)
	#	if (derivative=="yes") #something weird happens to the labels with odesolve
	#		colnames(x)<-c("time", "y1", "y2", "y3", "dy1", "dy2", "dy3")
	#}
	frame(); resetGraph();
	pairs(x, pch=15, gap=0, cex=0.2) # pch=183 does not work on UNIX, pch="." is too small
}

#function to restore working directory once demo is done
#onClose <- function() { setwd(oldwd); }
# Now handled by `.onClosePBSddeExamples`

createWin("demo_files/lorenz_win.txt")

})
