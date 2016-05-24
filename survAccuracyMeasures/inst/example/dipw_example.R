

#example running the dipw code

## example  
##

## 

source("R/dipw.R"); source("R/subroutines.R")
bw.power = .3; subcohort = TRUE
N = 500
c0 = qnorm(c(0.25,0.5,0.75)); n.cuts = length(c0)
SE0=0.9;  n.u0=length(SE0)
t0=1; beta=log(3)
ncch0.z0=ncch0.z1=ncch1.z0=ncch1.z1 = 50
data0=SIM.CCH.Z.FUN(N, DepCen=F,ncch0.z0,ncch0.z1, ncch1.z0,ncch1.z1)
phati=P0HAT.cch.z.FUN.type2(data0) 
#names(data0) = c("xi","di","yi","vi","zi","si","psi")
#data0: cohort data, cbind(xi,di,yi,vi,zi,si,psi)
data0$wi= 1/phati

EstROC.DIPW.NP.FUN(data0,u0=SE0,type="TPR",c0=c0,rtn="ALL")

## that works but now I want to calculate using package cohort data SimData

data(SimData)
head(SimData)

SimData$vi = 1
SimData$zi = 0

N = nrow(SimData)
SimData$si = SimData$status + 1

psi = rep(0, N)
psi[SimData$status==1] <- sum(SimData$status==1)/N
psi[SimData$status==0] <- sum(SimData$status==0)/N
SimData$psi = psi

names(SimData) = c("xi","di","yi","vi","zi","si","psi")

SimData$wi = 1
subcohort = FALSE
EstROC.DIPW.NP.FUN(SimData,u0=SE0,type="TPR",c0=0,rtn="ALL")


##

data(SimData)


survAM.estimate(time  = SimData$survTime, 
                   event = SimData$status, 
                   marker = SimData$Y, 
                   ESTmethod = "NP",
                   SEmethod = "normal",
                   predict.time = 1, cutpoint=0 )







