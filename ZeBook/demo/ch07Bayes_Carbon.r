
set.seed(3)
library(triangle)

##Data generation without trend

N<-40
Yc<-rep(NA,40)
Yc[1]<-rnorm(1, 16000,sqrt(100000))

for (i in 2:N) {
	
	Yc[i]<-Yc[i-1]-0.01*Yc[i-1]+0.2*2000+rnorm(1,0,sqrt(100000))
	
	}
	
Ym<-Yc+rnorm(N,0,sqrt(500000))


##Model function

CarbonMod<-function(Year,R,Sc0) {

Sc<-rep(NA,length(Year))

for (t in 1:length(Year)) {
Sc.t<-Sc0[t]
		
if (Year[t]>1) {	
for (y in 2:Year[t]) {
Sc.t<-Sc.t-R*Sc.t+0.2*2000}	
}

Sc[t]<-Sc.t
}

return(Sc)		
}

##############################################
##Simulation with initial value (prior mean)##
##############################################

Ysim<-Ym
for (i in 2:length(Yc)) {
Ysim[i]<-CarbonMod(i,0.02,Ym[1])
}

par(mfrow=c(1,1),oma=c(1,1,1,1))
plot(1:40,Ym,xlab="Year",ylab="Soil C kg ha-1",xlim=c(1,40),ylim=c(15000,25000))
lines(1:40,Ysim)

########################
##Estimation using OLS##
########################

Data<-data.frame(Time=1:40,Ym,Sc0=rep(Ym[1],40))
Data

Fit<-nls(Ym~CarbonMod(Time,R,Sc0), start=list(R=0.02), trace=T, data=Data)


Ysim<-Ym
for (i in 2:length(Yc)) {
Ysim[i]<-CarbonMod(i,coef(Fit),Ym[1])
}

lines(1:40,Ysim, lty=2)

##################
##Posterior mode##
##################

CarbonModPost<-function(Year,R,Sc0,SigmaYm,SigmaR) {

Sc<-rep(NA,length(Year))

for (t in 1:length(Year)) {
Sc.t<-Sc0[t]
		
if (Year[t]>1) {	
for (y in 2:Year[t]) {
Sc.t<-Sc.t-R*Sc.t+0.2*2000}	
}

Sc[t]<-Sc.t
}

return(c(Sc/SigmaYm,R/SigmaR))		
}


SigmaRvec<-sort(c(0.01,0.005, 0.001, 0.002,0.003,0.004, 0.0005,0.0004,0.0003,0.0002, 0.0001, 0.00001))
ParamVal<-SigmaRvec

for (k in 1:length(SigmaRvec)) {

SigmaR<-SigmaRvec[k]
YmPost<-c(Ym/sqrt(500000),0.02/SigmaR)
Time<-1:40
Sc0=rep(Ym[1],length(Time))
SigmaYm<-rep(sqrt(500000),length(Time))

Fit.post<-nls(YmPost~CarbonModPost(Time,R,Sc0,SigmaYm,SigmaR), start=list(R=0.02), trace=T)

ParamVal[k]<-coef(Fit.post)

}


dev.new()
par(mfrow=c(1,1), oma=c(5,5,5,5))

plot(SigmaRvec, ParamVal, xlab="Prior parameter standard deviation", ylab="Posterior mode", ylim=c(0.011,0.021), pch=19, type="b")
abline(h=coef(Fit), lty=2)
abline(h=0.02, lty=3)

#######################################
##Importance sampling. Gaussian prior##
#######################################
set.seed(1)

par(mfrow=c(2,2),oma=c(1,1,1,1))

Npara<-10

Rsample<-rnorm(Npara,0.02,0.01)
LogLike_vec<-rep(NA,Npara)


for (i in 1:Npara) {

Simul.i<-CarbonMod(Data$Time,Rsample[i],Data$Sc0)
LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,SigmaYm)))
LogLike_vec[i]<-LogLikelihood.i+1000
}

Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
Weight
plot(Rsample,Weight, xlim=c(0.005,0.025), pch=20, xlab="Value of the parameter R")
abline(v=ParamVal[SigmaRvec==0.01])
title("A. Q=10                             ")

Npara<-100

Rsample<-rnorm(Npara,0.02,0.01)
LogLike_vec<-rep(NA,Npara)


for (i in 1:Npara) {

Simul.i<-CarbonMod(Data$Time,Rsample[i],Data$Sc0)
LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,SigmaYm)))
LogLike_vec[i]<-LogLikelihood.i+1000
}

Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
Weight
plot(Rsample,Weight, xlim=c(0.005,0.025), pch=20, xlab="Value of the parameter R")
abline(v=ParamVal[SigmaRvec==0.01])
title("B. Q=100                            ")

Npara<-1000

Rsample<-rnorm(Npara,0.02,0.01)
LogLike_vec<-rep(NA,Npara)


for (i in 1:Npara) {

Simul.i<-CarbonMod(Data$Time,Rsample[i],Data$Sc0)
LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,SigmaYm)))
LogLike_vec[i]<-LogLikelihood.i+1000
}

Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
Weight
plot(Rsample,Weight, xlim=c(0.005,0.025), pch=20, xlab="Value of the parameter R")
abline(v=ParamVal[SigmaRvec==0.01])
title("C. Q=1,000                           ")

Npara<-10000

Rsample<-rnorm(Npara,0.02,0.01)
LogLike_vec<-rep(NA,Npara)


for (i in 1:Npara) {

Simul.i<-CarbonMod(Data$Time,Rsample[i],Data$Sc0)
LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,SigmaYm)))
LogLike_vec[i]<-LogLikelihood.i+1000
}

Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
Weight
plot(Rsample,Weight, xlim=c(0.005,0.025), pch=20, xlab="Value of the parameter R")
abline(v=ParamVal[SigmaRvec==0.01])
title("D. Q=10,000                           ")


dev.new()
par(mfrow=c(1,1), oma=c(5,1,5,1))
Rsample2<-Rsample[Weight>0]
Weight2<-Weight[Weight>0]
hist(sample(Rsample2,length(Rsample2), prob=Weight2, replace=T), xlab="Value of the parameter R", main=" ", xlim=c(0.01,0.015))

summary(sample(Rsample2,length(Rsample2), prob=Weight2, replace=T))


############################################
####Importance sampling. Triangle prior#####
############################################


Npara<-10000

Rsample<-rtriangle(Npara,0,  0.04, 0.02)
LogLike_vec<-rep(NA,Npara)


for (i in 1:Npara) {

Simul.i<-CarbonMod(Data$Time,Rsample[i],Data$Sc0)
LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,SigmaYm)))
LogLike_vec[i]<-LogLikelihood.i+1000
}

Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
Weight
dev.new()
par(mfrow=c(1,2))
plot(Rsample,Weight, xlim=c(0.005,0.025), pch=20, xlab="Value of the parameter R")
abline(v=ParamVal[SigmaRvec==0.01])
title("Triangle, N=10,000                           ")

Rsample2<-Rsample[Weight>0]
Weight2<-Weight[Weight>0]
hist(sample(Rsample2,length(Rsample2), prob=Weight2, replace=T), xlab="Value of the parameter R", main=" ", xlim=c(0.01,0.015))

summary(sample(Rsample2,length(Rsample2), prob=Weight2, replace=T))


#############################################################
####Importance sampling. Gaussian prior. Unknown variance####
#############################################################

set.seed(1)

#Number of parameter values
Npara<-10000

#Sampling from prior distributions
Rsample<-rnorm(Npara,0.02,0.01)
Ssample<-runif(Npara,0,2000)

#Initialisation of the vector including log likelihood
LogLike_vec<-rep(NA,Npara)

#Loop running the model and calculating log likelihood values

for (i in 1:Npara) {

Simul.i<-CarbonMod(Data$Time,Rsample[i],Data$Sc0)
LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,Ssample[i])))
LogLike_vec[i]<-LogLikelihood.i+1000

}

#Weight calculation
Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
Weight

dev.new()
par(mfrow=c(2,2))

plot(Rsample,Weight, xlim=c(0.005,0.025), pch=20, xlab="Value of the parameter R")
abline(v=ParamVal[SigmaRvec==0.01])
title("A.                                   ")

plot(Ssample,Weight, pch=20, xlab="Value of the residual standard deviation")
title("B.                                   ")

Rsample2<-Rsample[Weight>0]
Ssample2<-Ssample[Weight>0]
Weight2<-Weight[Weight>0]

Rsel<-sample(Rsample2,length(Rsample2), prob=Weight2, replace=T)
Ssel<-sample(Ssample2,length(Ssample2), prob=Weight2, replace=T)

hist(Rsel, xlab="Value of the parameter R", main=" ", xlim=c(0.01,0.015))
title("C.                                   ")
hist(Ssel, xlab="Value of the residual standard deviation", main=" ")
title("D.                                   ")

summary(Rsel)
summary(Ssel)

dev.new()
par(mfrow=c(1,2), oma=c(5,1,5,1))

Npara<-length(Rsel)
MeanVecR<-rep(NA,Npara)
Q10VecR<-rep(NA,Npara)
Q90VecR<-rep(NA,Npara)
MeanVecS<-rep(NA,Npara)
Q10VecS<-rep(NA,Npara)
Q90VecS<-rep(NA,Npara)

for (i in 10:Npara) {
MeanVecR[i]<-mean(Rsel[1:i])	
Q10VecR[i]<-quantile(Rsel[1:i],0.1, na.rm=T)	
Q90VecR[i]<-quantile(Rsel[1:i],0.9, na.rm=T)	
MeanVecS[i]<-mean(Ssel[1:i])	
Q10VecS[i]<-quantile(Ssel[1:i],0.1, na.rm=T)	
Q90VecS[i]<-quantile(Ssel[1:i],0.9, na.rm=T)	
}

plot(1:Npara,MeanVecR, xlab="Number of values", ylab="Parameter R", type="l", lwd=2)
title("A.                                    ")
#lines(1:Npara,Q10VecR, lty=1)
#lines(1:Npara,Q90VecR, lty=1)
plot(1:Npara,MeanVecS, xlab="Number of values", ylab="Standard deviation", type="l", lwd=2)
title("B.                                    ")
#lines(1:Npara,Q10VecS, lty=1)
#lines(1:Npara,Q90VecS, lty=1)

######################################################
##Metropolis-Hasting. Known variance. Gaussian prior##
######################################################

set.seed(1)

Npara<-20000
SigMH<-0.0015

R.1<-0.02

Rsample<-rep(NA, Npara)
Rsample[1]<-R.1

LogLike_vec<-rep(NA,Npara)
Simul.1<-CarbonMod(Data$Time,Rsample[1],Data$Sc0)
LogLikelihood.1<-sum(log(dnorm(Ym,Simul.1,SigmaYm)))
LogLike_vec[1]<-LogLikelihood.1

Count<-1

for (i in 2:Npara) {
	
	Rsample[i]<-Rsample[i-1]
	LogLike_vec[i]<-LogLike_vec[i-1]
	R.i<-rnorm(1,Rsample[i-1],SigMH)
	Simul.i<-CarbonMod(Data$Time,R.i,Data$Sc0)
	LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,SigmaYm)))
	Test.i<-LogLikelihood.i+log(dnorm(R.i,0.02,0.01))-LogLike_vec[i-1]-log(dnorm(Rsample[i-1],0.02,0.01))

	if (Test.i>0 | Test.i>log(runif(1,0,1))) {
		
	LogLike_vec[i]<-LogLikelihood.i
	Rsample[i]<-R.i
	Count<-Count+1	
										}
	}

print(100*Count/Npara)

par(mfrow=c(1,2), oma=c(5,1,5,1))
plot(1:100,Rsample[1:100], xlab="Iteration", ylab="Value of the parameter R", type="l", ylim=c(0.008,0.02))
title("A.                           ")
plot(1:Npara,Rsample, xlab="Iteration", ylab="Value of the parameter R", type="l", ylim=c(0.008,0.02))
title("B.                           ")

dev.new()
par(mfrow=c(1,1), oma=c(5,1,5,1))
acf(Rsample, main=" ")
summary(Rsample[seq(Npara/2,Npara, by=20)])



########################################################
##Metropolis-Hasting. Unknown variance. Gaussian prior##
########################################################

set.seed(3)

Npara<-20000
SigMH.R<-0.0015
SigMH.S<-0.2
R.1<-0.02
lS.1<-log(1500)

Rsample<-rep(NA, Npara)
Rsample[1]<-R.1
lSsample<-rep(NA, Npara)
lSsample[1]<-lS.1

LogLike_vec<-rep(NA,Npara)
Simul.1<-CarbonMod(Data$Time,Rsample[1],Data$Sc0)
LogLikelihood.1<-sum(log(dnorm(Ym,Simul.1,exp(lS.1))))
LogLike_vec[1]<-LogLikelihood.1

Count<-1

for (i in 2:Npara) {
	
	Rsample[i]<-Rsample[i-1]
	lSsample[i]<-lSsample[i-1]
	LogLike_vec[i]<-LogLike_vec[i-1]
	
	R.i<-rnorm(1,Rsample[i-1],SigMH.R)
	lS.i<-rnorm(1,lSsample[i-1],SigMH.S)
	Simul.i<-CarbonMod(Data$Time,R.i,Data$Sc0)
	LogLikelihood.i<-sum(log(dnorm(Ym,Simul.i,exp(lS.i))))
	
	Test.i<-LogLikelihood.i+log(dnorm(R.i,0.02,0.01))+log(dunif(exp(lS.i), 0, 2000))-LogLike_vec[i-1]-log(dnorm(Rsample[i-1],0.02,0.01))-log(dunif(exp(lSsample[i-1]), 0, 2000))

	if (Test.i>0 | Test.i>log(runif(1,0,1))) {
		
	LogLike_vec[i]<-LogLikelihood.i
	Rsample[i]<-R.i
	lSsample[i]<-lS.i
	Count<-Count+1	
										}
	}

print(100*Count/Npara)

par(mfrow=c(2,2), oma=c(1,1,1,1))
plot(1:100,Rsample[1:100], xlab="Iteration", ylab="Value of the parameter R", type="l", ylim=c(0.008,0.02))
title("A.                           ")
plot(1:Npara,Rsample, xlab="Iteration", ylab="Value of the parameter R", type="l", ylim=c(0.008,0.02))
title("B.                           ")
plot(1:100,exp(lSsample[1:100]), xlab="Iteration", ylab=expression(paste("Value of ",sigma)), type="l", ylim=c(500,2000))
title("C.                           ")
plot(1:Npara,exp(lSsample), xlab="Iteration", ylab=expression(paste("Value of ",sigma)), type="l", ylim=c(500,2000))
title("D.                           ")

dev.new()
par(mfrow=c(1,2), oma=c(5,1,5,1))

acf(Rsample, main="A. Parameter R    ")
acf(exp(lSsample), main=expression(paste("B. ", sigma,"          ")))

summary(Rsample[seq((Npara/4),Npara, by=20)])
sd(Rsample[seq((Npara/4),Npara, by=20)])
summary(exp(lSsample[seq((Npara/4),Npara, by=20)]))
sd(exp(lSsample[seq((Npara/4),Npara, by=20)]))

dev.new()
par(mfrow=c(1,2), oma=c(5,1,5,1))

hist(Rsample[seq((Npara/4),Npara, by=20)], xlab="Parameter R", main="A.                      ", xlim=c(0.01, 0.015))
hist(exp(lSsample[seq((Npara/4),Npara, by=20)]), xlab=expression(sigma), main="B.                        ", xlim=c(500,1700))

