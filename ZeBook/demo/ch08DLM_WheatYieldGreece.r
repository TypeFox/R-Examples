################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : Makowski. 2013-03-04
############################### MAIN PROGRAM ###################################
library(tseries)
library(dlm)
library(ZeBook)
################################################################################                               
#Data for Wheat Yield in Greece
TAB_Yield<-WheatYieldGreece

Year<-TAB_Yield[,1]
Yield<-TAB_Yield[,2]/10000

par(mfrow=c(1,1),oma=c(5,1,5,1))
plot(Year,Yield,ylab="Yield (t ha-1)", type="l",lwd=1)

################################################################################
#Definition of model 1
MyModel<-function(x) {
	return(dlmModPoly(1, dV=exp(x[1]), dW=exp(x[2])))
	}

#Estimation of parameters of the model
fitMyModel<-dlmMLE(Yield,parm=c(0,0), build=MyModel)
print(fitMyModel)

#aVar<-solve(fitMyModel$hessian)
#sqrt(diag(aVar))

#Filtrage, Lissage, Prediction

FittedModel<-MyModel(fitMyModel$par)
#FittedModel<-MyModel(c(0,-5))

YieldFilter<-dlmFilter(Yield, FittedModel)
YieldSmooth<-dlmSmooth(Yield, FittedModel)
YieldFilter_1<-YieldFilter

par(mfrow=c(1,1),oma=c(5,1,5,1))
plot(Year,Yield,ylab="Yield (t ha-1)", type="l",lwd=1)
lines(Year,YieldFilter$m[-1],lwd=2, lty=3)
lines(Year,YieldSmooth$s[-1],lwd=2)

################################################################################
#Definition of model 2
MyModel<-function(x) {
	return(dlmModPoly(2, dV=exp(x[1]), dW=c(exp(x[2]), exp(x[3]))))
	}

#Estimation of parameters of the model

fitMyModel<-dlmMLE(Yield,parm=c(0,0,0), build=MyModel)
print(fitMyModel)

#aVar<-solve(fitMyModel$hessian)
#sqrt(diag(aVar))

#Filtrage, Lissage, Prediction

FittedModel<-MyModel(fitMyModel$par)
#FittedModel<-MyModel(c(0,-5,-5))

YieldFilter<-dlmFilter(Yield, FittedModel)
YieldSmooth<-dlmSmooth(Yield, FittedModel)
YieldFilter_2<-YieldFilter


par(mfrow=c(1,1),oma=c(5,1,5,1))
plot(Year,Yield,ylab="Yield (t ha-1)", type="l",lwd=1)
lines(Year,YieldFilter$m[,1][-1]+YieldFilter$m[,2][-1],lwd=2, lty=3)
lines(Year,YieldSmooth$s[,1][-1]+YieldSmooth$s[,2][-1],lwd=2)


################################################################################
#Extraction of variances of slops smoothed
Var<-dlmSvd2var(YieldSmooth$U.S,YieldSmooth$D.S)
Var.slope.sm<-rep(NA,40)
for (i in 1:51) {
    Var.slope.sm[i]<-Var[[i]][2,2]
	}
                                  

par(mfrow=c(1,1),oma=c(5,1,5,1))
plot(Year,YieldSmooth$s[,2][-1],ylab="Yield increase rate (t ha-1 year-1)", type="l",lwd=2,ylim=c(-0.025,0.15))
lines(Year,YieldSmooth$s[,2][-1]+qnorm(0.75)*sqrt(Var.slope.sm[-1]),lty=2)
lines(Year,YieldSmooth$s[,2][-1]-qnorm(0.75)*sqrt(Var.slope.sm[-1]),lty=2)

################################################################################
#Predictions for next 20 years
foreYield_1<-dlmForecast(YieldFilter_1,nAhead=20)
foreYield_2<-dlmForecast(YieldFilter_2,nAhead=20)

FuturYear<-seq(Year[50]+1,Year[50]+20)

par(mfrow=c(1,2),oma=c(5,1,5,1))
plot(FuturYear,foreYield_1$f,xlab="Year", ylab="Forecasted Yield (t ha-1)", ylim=c(2,5), type="l")
#Confidence intervals
lines(FuturYear,foreYield_1$f+qnorm(0.75)*sqrt(unlist(foreYield_1$Q)),lty=2)
lines(FuturYear,foreYield_1$f-qnorm(0.75)*sqrt(unlist(foreYield_1$Q)),lty=2)
title("A                              ")

plot(FuturYear,foreYield_2$f,xlab="Year", ylab="Forecasted Yield (t ha-1)", ylim=c(2,5), type="l")

#Confidence intervals
lines(FuturYear,foreYield_2$f+qnorm(0.75)*sqrt(unlist(foreYield_2$Q)),lty=2)
lines(FuturYear,foreYield_2$f-qnorm(0.75)*sqrt(unlist(foreYield_2$Q)),lty=2)
title("B                              ")


TABslopeFrance<-data.frame(Year,YieldSmooth$s[,2][-1],Var.slope.sm[-1])
names(TABslopeFrance)<-c("Year","Slope","VarSlope")

TABslopeFrance

# end of file