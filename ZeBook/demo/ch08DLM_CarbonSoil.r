################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : David Makowski, 2013-03-04
############################### MAIN PROGRAM ###################################
library(dlm)
set.seed(3)

################################################################################
##Data generation without trend
N<-40
Yc<-rep(NA,40)
Yc[1]<-rnorm(1, 16000,sqrt(100000))
for (i in 2:N) {
	Yc[i]<-Yc[i-1]-0.01*Yc[i-1]+0.2*2000+rnorm(1,0,sqrt(100000))
	}
Ym<-Yc+rnorm(N,0,sqrt(500000))

################################################################################
# Model C 1
ModelC<-function(x) {
	mGG<-matrix(0,2,2)
	mGG[1,1]<-x[1]
	mGG[2,2]<-1
	mGG[1,2]<-0.2
	mW<-matrix(0,2,2)
	mW[1,1]<-exp(x[2])
	mW[2,2]<-0
	mC0<-matrix(0,2,2)
	mC0[1,1]<-500000
	mC0[2,2]<-0
	return(dlm(FF=matrix(c(1,0),1,2), V=500000, GG=mGG, W=mW, m0=c(16000,2000),C0=mC0))
	}

FittedModel<-ModelC(c(0.999, 13))

CarbFilter<-dlmFilter(Ym, FittedModel)
CarbSmooth<-dlmSmooth(Ym, FittedModel)
Smooth<-CarbSmooth$s

#Prediction
foreCarb<-dlmForecast(CarbFilter,nAhead=10)
Pred<-foreCarb$f

#Extraction of variances
Var<-dlmSvd2var(CarbSmooth$U.S,CarbSmooth$D.S)

VarC<-1:length(Ym)
for (i in 1:length(VarC)) {
	VarC[i]<-Var[[i]][1,1]
	}

#Graphics
par(mfrow=c(1,1),oma=c(1,1,1,1))

plot(1:40,Ym,xlab="Year",ylab="Soil C kg ha-1",xlim=c(1,50),ylim=c(15000,27000))
lines(1:N,Smooth[,1][-1],lwd=2)
lines(1:N,Smooth[,1][-1]+qnorm(0.95)*sqrt(VarC),lty=2)
lines(1:N,Smooth[,1][-1]-qnorm(0.95)*sqrt(VarC),lty=2)

lines(41:50,Pred)
lines(41:50,Pred+qnorm(0.95)*sqrt(unlist(foreCarb$Q)),lty=2)
lines(41:50,Pred-qnorm(0.95)*sqrt(unlist(foreCarb$Q)),lty=2)

################################################################################
## Model C 2 
ModelC<-function(x) {
	mGG<-matrix(0,2,2)
	mGG[1,1]<-x[1]
	mGG[2,2]<-1
	mGG[1,2]<-0.2
	mW<-matrix(0,2,2)
	mW[1,1]<-exp(x[2])
	mW[2,2]<-0
	mC0<-matrix(0,2,2)
	mC0[1,1]<-500000
	mC0[2,2]<-0
	return(dlm(FF=matrix(c(1,0),1,2), V=500000, GG=mGG, W=mW, m0=c(16000,2000),C0=mC0))
	}
	
fitTemp<-dlmMLE(Ym,parm=c(0.5,30), build=ModelC, hessian=T, lower=c(0, 1), upper=c(1, 100))
print(fitTemp)	

FittedModel<-ModelC(fitTemp$par)

CarbFilter<-dlmFilter(Ym, FittedModel)
CarbSmooth<-dlmSmooth(Ym, FittedModel)
Smooth<-CarbSmooth$s

#Prediction
foreCarb<-dlmForecast(CarbFilter,nAhead=10)
Pred<-foreCarb$f

#Extraction of variances
Var<-dlmSvd2var(CarbSmooth$U.S,CarbSmooth$D.S)

VarC<-1:length(Ym)

for (i in 1:length(VarC)) {
	VarC[i]<-Var[[i]][1,1]
	}

#Graphiques
	
par(mfrow=c(1,1),oma=c(1,1,1,1))

plot(1:40,Ym,xlab="Year",ylab="Soil C kg ha-1",xlim=c(1,50),ylim=c(15000,27000))
lines(1:N,Smooth[,1][-1],lwd=2)
lines(1:N,Smooth[,1][-1]+qnorm(0.95)*sqrt(VarC),lty=2)
lines(1:N,Smooth[,1][-1]-qnorm(0.95)*sqrt(VarC),lty=2)

lines(41:50,Pred)
lines(41:50,Pred+qnorm(0.95)*sqrt(unlist(foreCarb$Q)),lty=2)
lines(41:50,Pred-qnorm(0.95)*sqrt(unlist(foreCarb$Q)),lty=2)

# end of file