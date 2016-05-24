# Newton's Law of Cooling
# Models the cooling of a cup of coffee left in a room (ACB)
# see http://en.wikipedia.org/wiki/Heat_conduction#Newton.27s_law_of_cooling
require(PBSddesolve)
if (!require(PBSmodelling)) stop("The package PBSmodelling must be installed for this demo")

local(env=.PBSddeEnv, expr={

closeWin("window")

#store working directory
if (!exists("oldwd") || getwd() != system.file(package = "PBSddesolve")) {
	oldwd <- getwd()
	setwd(system.file(package = "PBSddesolve"))
}

runPlot <- function() {
	#extract variables from window
	getWinVal(scope="L")
	
	# `y' is the estimated value of the variable at time `t'
	# y[1] represents the temperature of a cup of coffee at a given time `t'.
	# dy/dt ~ y - room tempurature. 
	# so we can write dy/dt = -rho(y - roomTemp)
	myGrad <- function(t, y) {
		y1 <- -rho*(y[1]-Tenv)
		return(y1) 
	}
	#solve the ODE from t0..t1 - (ignore hbsize, setting to zero may crash)
	x <- dde(y=Tcup, func=myGrad, times=seq(t0,t1,length=100), hbsize=0) 
	frame(); resetGraph();
	plot(x, type="l", main="Cooling of a cup of coffee", ylab="Tempurature", xlab="Time")
}

#restore working directory once demo is done
#onClose <- function() { setwd(oldwd); }
# Now handled by `.onClosePBSddeExamples`

createWin("demo_files/cooling_win.txt")
})
