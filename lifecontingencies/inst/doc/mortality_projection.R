### R code from vignette source 'mortality_projection.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: mortality_projection.Rnw:66-67
###################################################
options(width=80, prompt='R> ')


###################################################
### code chunk number 2: load
###################################################
library(demography)
library(forecast)
library(lifecontingencies)


###################################################
### code chunk number 3: createDemogData
###################################################
#italyDemo<-hmd.mx(country="ITA", username="username@email.domain", 
#password="password", label="Italy")
load(file="mortalityDatasets.RData")


###################################################
### code chunk number 4: italyDemoFig
###################################################
par(mfrow=c(1,3))
plot(italyDemo,series="male",datatype="rate", main="Male rates")
plot(italyDemo,series="female",datatype="rate", main="Female rates")
plot(italyDemo,"total",datatype="rate", main="Total rates")


###################################################
### code chunk number 5: italyDemoFigTime
###################################################
par(mfrow=c(1,3))
plot(italyDemo,series="male",datatype="rate",
     plot.type="time", main="Male rates",xlab="Years")
plot(italyDemo,series="female",datatype="rate",
     plot.type="time", main="Female rates",xlab="Years")
plot(italyDemo,series="total",datatype="rate",
     plot.type="time", main="Total rates",xlab="Years")


###################################################
### code chunk number 6: fitLeeCarter
###################################################
italyLcaM<-lca(italyDemo,series="male",max.age=100)
italyLcaF<-lca(italyDemo,series="female",max.age=100)
italyLcaT<-lca(italyDemo,series="total",max.age=100)


###################################################
### code chunk number 7: leeCarterResultsFig
###################################################
  par(mfrow=c(1,3))
  plot(italyLcaT$ax, main="ax", xlab="Age",ylab="ax",type="l")
  lines(x=italyLcaF$age, y=italyLcaF$ax, main="ax", col="red")
  lines(x=italyLcaM$age, y=italyLcaM$ax, main="ax", col="blue")
  legend("topleft" , c("Male","Female","Total"),
  cex=0.8,col=c("blue","red","black"),lty=1);
  plot(italyLcaT$bx, main="bx", xlab="Age",ylab="bx",type="l")
  lines(x=italyLcaF$age, y=italyLcaF$bx, main="bx", col="red")
  lines(x=italyLcaM$age, y=italyLcaM$bx, main="bx", col="blue")
  legend("topright" , c("Male","Female","Total"),
  cex=0.8,col=c("blue","red","black"),lty=1);
  plot(italyLcaT$kt, main="kt", xlab="Year",ylab="kt",type="l")
  lines(x=italyLcaF$year, y=italyLcaF$kt, main="kt", col="red")
  lines(x=italyLcaM$year, y=italyLcaM$kt, main="kt", col="blue")
  legend("topright" , c("Male","Female","Total"),
  cex=0.8,col=c("blue","red","black"),lty=1);


###################################################
### code chunk number 8: ktProjections
###################################################
fM<-forecast(italyLcaM,h=110)
fF<-forecast(italyLcaF,h=110)
fT<-forecast(italyLcaT,h=110)


###################################################
### code chunk number 9: ktProjectionFig
###################################################
par(mfrow=c(1,3))
plot(fM$kt.f,main="Male")
plot(fF$kt.f,main="Female",)
plot(fT$kt.f,main="Total")


###################################################
### code chunk number 10: ktrates
###################################################
ratesM<-cbind(italyDemo$rate$male[1:100,],fM$rate$male[1:100,])
ratesF<-cbind(italyDemo$rate$female[1:100,],fF$rate$female[1:100,])
ratesT<-cbind(italyDemo$rate$total[1:100,],fT$rate$total[1:100,])


###################################################
### code chunk number 11: ktratesFig
###################################################
par(mfrow=c(1,1))
plot(seq(min(italyDemo$year),max(italyDemo$year)+110),ratesF[65,],
     col="red",xlab="Years",ylab="Death Rates",type="l")
lines(seq(min(italyDemo$year),max(italyDemo$year)+110),ratesM[65,],
      col="blue",xlab="Years",ylab="Death Rates")
lines(seq(min(italyDemo$year),max(italyDemo$year)+110),ratesT[65,],
      col="black",xlab="Years",ylab="Death Rates")
legend("topright" , c("Male","Female","Total"),
       cex=0.8,col=c("blue","red","black"),lty=1);


###################################################
### code chunk number 12: lifeTableProject
###################################################

createActuarialTable<-function(yearOfBirth,rate){

  mxcoh <- rate[1:nrow(rate),(yearOfBirth-min(italyDemo$year)+1):ncol(rate)]
  cohort.mx <- diag(mxcoh)
  cohort.px=exp(-cohort.mx)
  #get projected Px
  fittedPx=cohort.px #add px to table
	px4Completion=seq(from=cohort.px[length(fittedPx)], to=0, length=20)
	totalPx=c(fittedPx,px4Completion[2:length(px4Completion)])
	#create life table
	irate=1.04/1.02-1

	cohortLt=probs2lifetable(probs=totalPx, radix=100000,type="px", 
  name=paste("Cohort",yearOfBirth))
	cohortAct=new("actuarialtable",x=cohortLt@x, lx=cohortLt@lx, 
	interest=irate, name=cohortLt@name)
	return(cohortAct)
	}




###################################################
### code chunk number 13: annuityAPV
###################################################
	getAnnuityAPV<-function(yearOfBirth,rate) {
		actuarialTable<-createActuarialTable(yearOfBirth,rate)
		out=axn(actuarialTable,x=65,m=12)
		return(out)
	}
rate<-ratesM
for(i in seq(1920,2000,by=10)) {
		cat("For cohort ",i, "of males the e0 is",
		round(exn(createActuarialTable(i,rate)),2),
		" and the APV is :",round(getAnnuityAPV(i,rate),2),"\n")
		
	}
rate<-ratesF
for(i in seq(1920,2000,by=10)) {
  	cat("For cohort ",i, "of females the e0 at birth is",
	round(exn(createActuarialTable(i,rate)),2),
	" and the APV is :",round(getAnnuityAPV(i,rate),2),"\n")
		
	}
rate<-ratesT
for(i in seq(1920,2000,by=10)) {
    cat("For cohort ",i, "of total population the e0 is",
		round(exn(createActuarialTable(i,rate)),2),
		" and the APV is :",round(getAnnuityAPV(i,rate),2),"\n")
		
	}


