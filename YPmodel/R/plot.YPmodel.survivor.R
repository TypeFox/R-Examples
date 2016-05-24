plot.YPmodel.survivor <-
function(x=c(), Internal=c(), ...)
########################################################
#fun.plot1(bt,ru,fb,Data)
#######################################################
# version 0.1
# May 19, 2012
# Junlong Sun
# no-return
#######################################################
# May 19, 2012 - v0.1 Create
#necessary Data, InternalParameters
#######################################################
{
	Data <- x$Data
	Parameters  <- x$Parameters

	#-----------------------------------------------------------------#

	if(is.null(Data)){
		stop(paste(fun.errorMessage('DataSet')))
	}
	#-----------------------------------------------------------------#
		
	if(is.null(Internal)){
		if(is.null(x)){
			if(is.null(Parameters)){
				warning(paste(fun.errorMessage('DefaultParameter')))
				Parameters <- YPmodel.setParameter()
			}
			x <- YPmodel.estimate(Data=Data)
		}
		Estimate <- x
		if(is.null(Parameters)){
			warning(paste(fun.errorMessage('DefaultParameter')))
			Parameters <- YPmodel.setParameter()
		}		
		Internal <- fun.internalParameters(Data=Data, Estimate=Estimate)
	}
	bt <- Internal$bt
	ru <- Internal$ru
	fb <- Internal$fb

#-----------------------------------------------------------------#
## loading data
#-----------------------------------------------------------------#
	fb1 <- fb$Num1
	fb2 <- fb$Num2
	X <- Data$X

#-----------------------------------------------------------------#
## main function
#-----------------------------------------------------------------#

m <- 1

rr <- 1/bt
r <- ru

lam2 <- log(1+bt[2]/bt[1]*r)/bt[2]

fbh2 <- exp(-lam2)
fbh1 <- 1/(1+r)

plot(X*365,fb1,"s",col="black",xlab="Days", ylab="Survival",)
lines(X*365,fb2,"s",col="red" , lwd=1.5)
lines(X*365,fbh1,"s",col="blue",lty=4 , lwd=1.5)
lines(X*365,fbh2,"s",col="magenta",lty=4 , lwd=1.5)
legend("topright", c("KM, Group1","KM, Group2","Fitted, Group1","Fitted, Group2"),lty=c(1,1,4,4),col=c("black","red","blue","magenta")) 
title(main="Survivor function estimates and model-based estimates")

##return [p1,q1]
#data8 <- fun.stairs(oy,fb1)
#p1 <- data8$p
#q1 <- data8$q

##return [p2,q2]
#data9 <- fun.stairs(oy,fb2)
#p2 <- data9$p
#q2 <- data9$q

##return [p3,q3]
#data10 <- fun.stairs(oy,fbh1)
#p3 <- data10$p
#q3 <- data10$q

##return [p4,q4]
#data11 <- fun.stairs(oy,fbh2)
#p4 <- data11$p
#q4 <- data11$q

#plot(p1,q1,'k-', p2,q2,'r-', p3,q3,'b-.', p4,q4,'m:','LineWidth',1.5);
#plot(p1,q1,"l",col="black")
#points(p2,q2,col="red")
#points(p3,q3,col="blue",lty=2)
#points(p4,q4,col="magenta",lty=2)

#dev.new()
#hr=(1+r)/(1/rr[1]+r/rr[2])
#plot(X,hr,"l",col="yellow",xlab="Days", ylab="Ration function")

}
