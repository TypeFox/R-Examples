plot.YPmodel.IntervalBands <-
function(x=c(), Internal=c(), ...)
#fun.plot2 <-
#function(b,bt,ru,gama,p,pl,deni,sm,kall,Data,GroupData)
########################################################
#fun.plot2(b,bt,ru,gama,p,pl,deni,sm,kall,Data,GroupData)
#######################################################
# version 0.1
# May 19, 2012
# Junlong Sun
# no-return
#######################################################
# May 19, 2012 - v0.1 Create
#######################################################
{

	IntervalBands <- x
	Estimate <- x$Estimate
	Data <- x$Data

	if(is.null(Internal)){
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}

	hr <- IntervalBands$hr
	ld2 <- IntervalBands$ld2
	ud2 <- IntervalBands$ud2
	upp22 <- IntervalBands$upp22
	low22 <- IntervalBands$low22
	upp3 <- IntervalBands$upp3
	low3 <- IntervalBands$low3
	upp90 <- IntervalBands$upp90
	low90 <- IntervalBands$low90
	GroupData <- Data$GroupData

#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
n <- Data$length
Z <- Data$Z
X <- Data$X
Delta <- Data$Delta

#-----------------------------------------------------------------#
## main function
#-----------------------------------------------------------------#

ylimMax=max(upp22[ld2:ud2],upp3[ld2:ud2],hr[ld2:ud2],upp90[ld2:ud2])

#dev.new()

#plot(X[ld2:ud2]*365,upp22[ld2:ud2],"l",col="blue",xlab="Days", ylab="Hazard ratio",ylim=c(0, ylimMax))
#lines(X[ld2:ud2]*365,low22[ld2:ud2],"l",col="blue")
#lines(X[ld2:ud2]*365,upp3[ld2:ud2],"l",col="magenta",lty=3)
#lines(X[ld2:ud2]*365,low3[ld2:ud2],"l",col="magenta",lty=3)
#lines(X[ld2:ud2]*365,hr[ld2:ud2],"l",col="green",lty=3)
#lines(X[ld2:ud2]*365,upp90[ld2:ud2],"l",col="cyan")
#lines(X[ld2:ud2]*365,low90[ld2:ud2],"l",col="cyan")
#legend("topright", c("Estimated hazard ratio function","Pointwise confidence limits","95% EQ confidence bands","90% EQ confidence bands"),lty=c(3,3,1,1),col=c("green","magenta","blue","cyan"))
#title(main="Confidence intervals and bands of the hazard ration")

plot(X[ld2:ud2]*365,upp22[ld2:ud2],"l",col="blue",xlab="Days", ylab="Hazard ratio",ylim=c(0, ylimMax), lwd=1.5)
lines(X[ld2:ud2]*365,low22[ld2:ud2],"l",col="blue", lwd=1.5)
lines(X[ld2:ud2]*365,upp3[ld2:ud2],"l",col="magenta",lty=2, lwd=1.5)
lines(X[ld2:ud2]*365,low3[ld2:ud2],"l",col="magenta",lty=2, lwd=1.5)
lines(X[ld2:ud2]*365,hr[ld2:ud2],"l",col="green",lty=2, lwd=1.5)
#lines(X[ld2:ud2]*365,upp90[ld2:ud2],"l",col="cyan")
#lines(X[ld2:ud2]*365,low90[ld2:ud2],"l",col="cyan")
legend("topright", c("Estimated hazard ratio function","Pointwise confidence limits","95% EQ confidence bands"),lty=c(2,2,1),col=c("green","magenta","blue"))
title(main="Confidence intervals and bands of the hazard ration")


}
