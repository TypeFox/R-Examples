################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
################################################################################
# log of concentrate likelihood for LAI and Biomass
#maize.logConclikelihood.optim(c(1.80),c("RUE"),maize.define.param()["nominal",],maize.data_EuropeEU,list_n_sy, transf=log)
maize.logConclikelihood.optim<-function(param,param_opti,param_default,data_sy,list_sy, transf=function(x){x}){
    param_all = param_default
    param_all[param_opti] <- param
    sim=maize.multisy(param_all,list_sy, 100, 250)
    simobs=merge(sim,data_sy,by=c("sy","day"))
    n_Bobs=length(na.omit(simobs$Bobs))
    n_LAIobs=length(na.omit(simobs$LAIobs))
    MSE_B=sum((transf(simobs$Bobs)-transf(simobs$B))^2,na.rm=TRUE)/n_Bobs
    MSE_LAI=sum((transf(simobs$LAIobs)-transf(simobs$LAI))^2,na.rm=TRUE)/n_LAIobs
    Conclikelihood = log(MSE_B^(n_Bobs/2)) + log( MSE_LAI^(n_LAIobs/2))
    return(Conclikelihood)
    }
################################################################################
list_n_sy=unique(maize.data_EuropeEU$sy)
# 2) Use supplementary data with concentrate likelihood without any transformation.
# Why ? 2 variables, LAI and Bobs with quite different magnitude
# test the use of a log tranformation too
# up to 100 minutes !
param_default=maize.define.param()["nominal",]
# 1 parameter
param_opti=c("RUE")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(optim_ConcL1p<-optim(param_init, maize.logConclikelihood.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(optim_ConcL1p_a<-optim(param_init, maize.logConclikelihood.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))

if(optim_ConcL1p_a$value<optim_ConcL1p$value){optim_ConcL1p<-optim_ConcL1p_a}
optim_ConcL1p
save(optim_ConcL1p,file="maize.case_03_optim_ConcL1p.rda")

param_all = param_default
param_all[param_opti] <- optim_ConcL1p$par
simobs=merge(maize.multisy(param_all,list_n_sy, 100, 250),maize.data_EuropeEU,by=c("sy","day"))
AIC_B_ConcL1p=AICf(simobs$Bobs,simobs$B,length(param_opti)+1)
goodness.of.fit(simobs$Bobs,simobs$B)

AIC_LAI_ConcL1p=AICf(simobs$LAIobs,simobs$LAI,length(param_opti)+1)
goodness.of.fit(simobs$LAIobs,simobs$LAI)

simobs240=subset(simobs,day==240)
AIC_B240_ConcL1p=AICf((simobs240$Bobs),(simobs240$B),length(param_opti)+1)
goodness.of.fit(simobs240$Bobs,simobs240$B)

# analysis of residuals
resB=simobs$Bobs-simobs$B
par(mfrow=c(1,2))
plot(simobs$B,resB, xlab="predicted values of B",ylab="residuals of B")
abline(a=0,b=0)
resLAI=simobs$LAIobs-simobs$LAI
plot(simobs$LAI,resLAI, xlab="predicted values of LAI",ylab="residuals of LAI")
abline(a=0,b=0)

# Kolmogorov-Smirnov Tests of normality
ks.test(resB,"pnorm",mean(resB),sd(resB))
ks.test(na.omit(resLAI),"pnorm",mean(na.omit(resLAI)),sd(na.omit(resLAI)))
# Bartlett Test of Homogeneity of Variances
bartlett.test(resB[order(simobs$B)],cut(sort(simobs$B), breaks=seq(from=min(simobs$B),to=max(simobs$B),length.out=5),include.lowest =TRUE))
bartlett.test(resLAI[order(simobs$LAI)],cut(sort(simobs$LAI), breaks=seq(from=min(simobs$LAI),to=max(simobs$LAI),length.out=3),include.lowest =TRUE),na.action=na.omit())
################################################################################
# 2 parameters
param_opti=c("RUE","TTM")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(optim_ConcL2p<-optim(param_init,  maize.logConclikelihood.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(optim_ConcL2p_a<-optim(param_init,  maize.logConclikelihood.optim , method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))

if(optim_ConcL2p_a$value<optim_ConcL2p$value){optim_ConcL2p<-optim_ConcL2p_a}
optim_ConcL2p
save(optim_ConcL2p,file="maize.case_03_optim_ConcL2p.rda")

param_all = param_default
param_all[param_opti] <- optim_ConcL2p$par
simobs=merge(maize.multisy(param_all,list_n_sy, 100, 250),maize.data_EuropeEU,by=c("sy","day"))
AIC_B_ConcL2p=AICf(simobs$Bobs,simobs$B,length(param_opti)+1)
goodness.of.fit(simobs$Bobs,simobs$B)

AIC_LAI_ConcL2p=AICf(simobs$LAIobs,simobs$LAI,length(param_opti)+1)
goodness.of.fit(simobs$LAIobs,simobs$LAI)

simobs240=subset(simobs,day==240)
AIC_B240_ConcL2p=AICf((simobs240$Bobs),(simobs240$B),length(param_opti)+1)
goodness.of.fit(simobs240$Bobs,simobs240$B)
################################################################################
# 3 parameters
param_opti=c("RUE","TTM","alpha")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(optim_ConcL3p<-optim(param_init,  maize.logConclikelihood.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))
optim_ConcL3p

# second starting point
param_init <- param_default[param_opti]*1.25
system.time(optim_ConcL3p_a<-optim(param_init,  maize.logConclikelihood.optim , method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))
optim_ConcL3p_a

if(optim_ConcL3p_a$value<optim_ConcL3p$value){optim_ConcL3p<-optim_ConcL3p_a}
optim_ConcL3p
# criterion > criterion obtain for 2 param => pbl of local minimum

# third starting point : taking into account values obtained for 2 parameters
param_init <- param_default[param_opti]
param_init[c("RUE","TTM")] <- optim_ConcL2p$par

system.time(optim_ConcL3p_b<-optim(param_init,  maize.logConclikelihood.optim , method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy, transf=function(x){x} ))
optim_ConcL3p_b

if(optim_ConcL3p_b$value<optim_ConcL3p$value){optim_ConcL3p<-optim_ConcL3p_b}
optim_ConcL3p
save(optim_ConcL3p,file="maize.case_03_optim_ConcL3p.rda")

param_all = param_default
param_all[param_opti] <- optim_ConcL3p$par
simobs=merge(maize.multisy(param_all,list_n_sy, 100, 250),maize.data_EuropeEU,by=c("sy","day"))
AIC_B_ConcL3p=AICf(simobs$Bobs,simobs$B,length(param_opti)+1)
goodness.of.fit(simobs$Bobs,simobs$B)

AIC_LAI_ConcL3p=AICf(simobs$LAIobs,simobs$LAI,length(param_opti)+1)
goodness.of.fit(simobs$LAIobs,simobs$LAI)

simobs240=subset(simobs,day==240)
AIC_B240_ConcL3p=AICf((simobs240$Bobs),(simobs240$B),length(param_opti)+1)
goodness.of.fit(simobs240$Bobs,simobs240$B)

data.frame(param_opti=c("RUE", "RUE&TTM", "RUE&TTM&alpha"),
AIC_LAI=c(AIC_LAI_ConcL1p["AICcomplete"],AIC_LAI_ConcL2p["AICcomplete"],AIC_LAI_ConcL3p["AICcomplete"]),
AIC_B=c(AIC_B_ConcL1p["AICcomplete"],AIC_B_ConcL2p["AICcomplete"],AIC_B_ConcL3p["AICcomplete"]),
AIC_B240=c(AIC_B240_ConcL1p["AICcomplete"],AIC_B240_ConcL2p["AICcomplete"],AIC_B240_ConcL3p["AICcomplete"]))
# end of script