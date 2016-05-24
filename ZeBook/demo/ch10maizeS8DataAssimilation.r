################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
library(dlm)

# choose a year and a site
year=2010
idsite=92

# build a table with ALL simulation and observation
simobs = merge(maize.data_EuropeEU[(maize.data_EuropeEU$idsite==idsite)&(maize.data_EuropeEU$year==year),],maize.model2(maize.define.param()["nominal",],maize.weather(working.year=year,working.site=idsite,weather_all=weather_EuropeEU), 100, 250),  by=c("day"), all.y = TRUE)[,c("day","LAIobs","Bobs","TT","LAI","B")]
simobs

# function of DLM
x<-matrix(simobs$LAI,ncol=1)
buildFun<-function(theta) {
	modMAIZE<-dlmModReg(x,dV=0.75,dW=c(exp(theta[1]),exp(theta[2])))
	return(modMAIZE)
	}
Y<-simobs$LAIobs
fit<-dlmMLE(Y,parm=c(0,0),build=buildFun)
fit
fitted.modMAIZE<-buildFun(fit$par)
modFilter<-dlmFilter(Y,mod=fitted.modMAIZE)
modSmooth<-dlmSmooth(Y,mod=fitted.modMAIZE)
SmoothedLAI<-modSmooth$s[,1][-1]+modSmooth$s[,2][-1]*x
SmoothedLAI[SmoothedLAI<0]<-0
plot(1:length(x),Y,ylim=c(0,10),pch=19, xlab="Time (days)", ylab="LAI")
lines(1:length(x),x,lty=1, lwd=2)
lines(1:length(x),modFilter$m[-1,1]+modFilter$m[-1,2]*x, lty=2, lwd=3)
lines(1:length(x),SmoothedLAI, lty=3, lwd=1)

# end of file
