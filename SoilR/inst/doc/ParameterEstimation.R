### R code from vignette source 'ParameterEstimation.Rnw'

###################################################
### code chunk number 1: ParameterEstimation.Rnw:58-64
###################################################
library(SoilR)
library(FME)
library(MASS)
library(lattice)
BorealCO2=subset(eCO2, subset=Sample=="AK_T25", select=-Sample)
names(BorealCO2)<-c("time","eCO2","eCO2sd")


###################################################
### code chunk number 2: ParameterEstimation.Rnw:71-74
###################################################
plot(BorealCO2[,1:2], xlab="Days", ylab="Evolved CO2 (mgC g-1 soil)")
arrows(BorealCO2[,1],BorealCO2[,2]-BorealCO2[,3],BorealCO2[,1],
       BorealCO2[,2]+BorealCO2[,3],code=3,angle=90,length=0.1)


###################################################
### code chunk number 3: ParameterEstimation.Rnw:89-105
###################################################
days=seq(0,42)
Ctotal=7.7

eCO2func=function(pars){
  mod=TwopFeedbackModel(
  t=days,
  ks=pars[1:2],
  a21=pars[3]*pars[1],
  a12=pars[4]*pars[2], 
  C0=Ctotal*c(pars[5],1-pars[5]), 
  In=0,
  pass=TRUE
  )
  AccR=getAccumulatedRelease(mod)
  return(data.frame(time=days,eCO2=rowSums(AccR)))
}


###################################################
### code chunk number 4: ParameterEstimation.Rnw:112-116
###################################################
eCO2cost=function(pars){
  modelOutput=eCO2func(pars)
  return(modCost(model=modelOutput, obs=BorealCO2, err="eCO2sd"))
}


###################################################
### code chunk number 5: ParameterEstimation.Rnw:122-126
###################################################
inipars=c(k1=0.5,k2=0.05,alpha21=0.5,alpha12=0.1,gamma=0.5)

eCO2fit=modFit(f=eCO2cost,p=inipars,method="Marq",
               upper=c(Inf,Inf,1,1,1),lower=c(0,0,0,0,0))


###################################################
### code chunk number 6: ParameterEstimation.Rnw:130-131
###################################################
eCO2fit$par


###################################################
### code chunk number 7: ParameterEstimation.Rnw:136-137
###################################################
fitmod=eCO2func(eCO2fit$par)


###################################################
### code chunk number 8: ParameterEstimation.Rnw:142-146
###################################################
plot(BorealCO2[,1:2], xlab="Days", ylab="Evolved CO2 (mgC g-1 soil)")
arrows(BorealCO2[,1],BorealCO2[,2]-BorealCO2[,3],BorealCO2[,1],
       BorealCO2[,2]+BorealCO2[,3],code=3,angle=90,length=0.1)
lines(fitmod)


###################################################
### code chunk number 9: ParameterEstimation.Rnw:154-159
###################################################
var0=eCO2fit$var_ms_unweighted

eCO2mcmc=modMCMC(f=eCO2cost, p=eCO2fit$par, niter=1000, jump=var0,  
                 var0=var0, wvar0=0.1, updatecov=50, lower=c(0,0,0,0,0),
                 upper=c(1,1,1,1,1))


###################################################
### code chunk number 10: ParameterEstimation.Rnw:165-166
###################################################
summary(eCO2mcmc)


###################################################
### code chunk number 11: ParameterEstimation.Rnw:173-174
###################################################
pairs(eCO2mcmc)


###################################################
### code chunk number 12: ParameterEstimation.Rnw:184-190
###################################################
predRange=sensRange(func=eCO2func, parInput=eCO2mcmc$par)
plot(summary(predRange),ylim=c(0,9),xlab="Days",
     ylab="Evolved CO2 (mg C g-1 C)",main="")
points(BorealCO2)
arrows(BorealCO2[,1],BorealCO2[,2]-BorealCO2[,3],BorealCO2[,1],
       BorealCO2[,2]+BorealCO2[,3],code=3,angle=90,length=0.1)


###################################################
### code chunk number 13: ParameterEstimation.Rnw:204-206
###################################################
plot(D14C~Year,data=HarvardForest14CO2,
     ylab=expression(paste(Delta^14,"C ","(\u2030)")))


###################################################
### code chunk number 14: ParameterEstimation.Rnw:228-231
###################################################
time=C14Atm_NH$YEAR
t_start=min(time)
t_end=max(time)


###################################################
### code chunk number 15: ParameterEstimation.Rnw:236-241
###################################################
inputFluxes=BoundInFlux(
	function(t0){matrix(nrow=3,ncol=1,c(270,150,0))},
	t_start,
	t_end
)


###################################################
### code chunk number 16: ParameterEstimation.Rnw:246-247
###################################################
C0=c(390,220+390+1376,90+1800+560) 


###################################################
### code chunk number 17: ParameterEstimation.Rnw:253-272
###################################################
Fc=BoundFc(C14Atm_NH,lag=0,format="Delta14C")
Mod1<-function(ks,pass=TRUE){
  At=ConstLinDecompOp(
           matrix(nrow=3,ncol=3,byrow=TRUE,c(ks[1],0,0,
                                             ks[4],ks[2],0,
                                             ks[5],0,ks[3]))
  ) 
  mod=GeneralModel_14(
	t=time,
	A=At,
	ivList=C0,
	initialValF=ConstFc(rep(0,3),"Delta14C"),
	inputFluxes=inputFluxes,
	inputFc=Fc,
	pass=TRUE
  ) 
  R14t=getF14R(mod)
  return(data.frame(time=time,R14t=R14t))
}


###################################################
### code chunk number 18: ParameterEstimation.Rnw:276-279
###################################################
DataR14t=cbind(time=HarvardForest14CO2[,1],
               R14t=HarvardForest14CO2[,2],
               sd=sd(HarvardForest14CO2[,2]))


###################################################
### code chunk number 19: ParameterEstimation.Rnw:285-303
###################################################
#Create the cost function
R14tCost <- function(pars){
  R14t <- Mod1(pars)
  return(modCost(model=R14t,obs=DataR14t,err="sd"))
}

#Fit the model to the observed data given some initial value for the parameters
Fit <- modFit(f=R14tCost,p=c(-0.5,-0.9,-0.1,0.1,0.1))
#
# Run an MCMC using the variance and covariance results from the previous optimization
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled 
MCMC <- modMCMC(f=R14tCost, p = Fit$par, niter = 1000, jump = NULL, var0 = var0, wvar0 = 0, 
                lower=c(-3,-3,-1,0,0),upper=c(0,0,0,1,1))

#The sensitivity range is calculated from the output of the MCMC
sR=sensRange(func=Mod1, parInput=MCMC$par)



###################################################
### code chunk number 20: ParameterEstimation.Rnw:309-310
###################################################
pairs(MCMC,nsample=500)


###################################################
### code chunk number 21: ParameterEstimation.Rnw:319-324
###################################################
par(mar=c(5,5,4,1))
plot(summary(sR),xlim=c(1950,2010),ylim=c(0,1000),xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),main="")
points(DataR14t,pch=20)
lines(C14Atm_NH,col=4)


