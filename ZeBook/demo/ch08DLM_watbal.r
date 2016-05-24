################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-04-17
############################### MAIN PROGRAM ###################################
library(dlm)
library(ZeBook)

TAB=watbal.simobsdata

#x : simulation for each day
x<-matrix(TAB$WATp,ncol=1)

#variance empirique

Vemp<-mean(TAB$WATp_SF.var[!is.na(TAB$WATp_SF.var)])


buildFun<-function(theta) {
	modWAT<-dlmModReg(x,dV=Vemp,dW=c(exp(theta[1]),exp(theta[2])))
	return(modWAT)
	}

# Y : observation
Y<-TAB$WATp_SF.mean

# Parameter estimation
fit<-dlmMLE(Y,parm=c(0,0),build=buildFun)
fit

#
fitted.modWAT<-buildFun(fit$par)

#
modFilter<-dlmFilter(Y,mod=fitted.modWAT)
modFilter$m

# graphical representation
plot(1:length(x),100*Y,ylim=c(0,50),pch=19, xlab="Time (days)", ylab="Soil water content (%)")
lines(1:length(x),100*x,lty=2)
lines(1:length(x),100*(modFilter$m[-1,1]+modFilter$m[-1,2]*x))
