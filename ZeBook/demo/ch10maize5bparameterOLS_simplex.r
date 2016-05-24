################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
################################################################################
#maize.LS.B240.optim(c(1.80),c("RUE"),maize.define.param()["nominal",],maize.data_EuropeEU,list_n_sy)
maize.LS.B240.optim<-function(param,param_opti,param_default,data_sy,list_sy){
    param_all = param_default
    param_all[param_opti] <- param
    sim=maize.multisy(param_all,list_sy, 100, 250)
    simobs=merge(sim,data_sy,by=c("sy","day"))
    return(sum((simobs$Bobs-simobs$B)^2,na.rm=TRUE))
    }
################################################################################
list_n_sy=unique(maize.data_EuropeEU$sy)
# 1)  3.5.1.	Using OLS (ordinary least squares), based on the final biomass (day=240) data.
data_n_B240=subset(maize.data_EuropeEU,day==240)
data_n_B240$LAIobs=NA
################################################################################
# 1b) Nelder Nead (simplex) method (optim)
# gradient methods (as nls) are often not possible for this type of model with treshold (TTM, TTL).
param_default=maize.define.param()["nominal",]
# 1 parameter
param_opti=c("RUE")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(optim_OLS1p<-optim(param_init, maize.LS.B240.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(optim_OLS1p_b<-optim(param_init, maize.LS.B240.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))

if(optim_OLS1p_b$value<optim_OLS1p$value){optim_OLS1p<-optim_OLS1p_b}
optim_OLS1p
save(optim_OLS1p,file="maize.case_03_optim_OLS1p.rda")

# compute AIC
param_all = param_default
param_all[param_opti] <- optim_OLS1p$par
simobs240=merge(maize.multisy(param_all,list_n_sy, 100, 250),data_n_B240,by=c("sy","day"))
AIC_B_1p=AICf(simobs240$Bobs,simobs240$B,length(param_opti)+1)
goodness.of.fit(simobs240$Bobs,simobs240$B)

# 2 parameters
param_opti=c("RUE","TTM")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(optim_OLS2p<-optim(param_init,  maize.LS.B240.optim, method = "Nelder-Mead",control=list(trace=1, maxit=15000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(optim_OLS2p_a<-optim(param_init,  maize.LS.B240.optim , method = "Nelder-Mead",control=list(trace=1, maxit=15000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))

if(optim_OLS2p_a$value<optim_OLS2p$value){optim_OLS2p<-optim_OLS2p_a}
optim_OLS2p
save(optim_OLS2p,file="maize.case_03_optim_OLS2p.rda")

# compute AIC
param_all = param_default
param_all[param_opti] <- optim_OLS2p$par
simobs240=merge(maize.multisy(param_all,list_n_sy, 100, 250),data_n_B240,by=c("sy","day"))
AIC_B_2p=AICf(simobs240$Bobs,simobs240$B,length(param_opti)+1)
goodness.of.fit(simobs240$Bobs,simobs240$B)

# 3 parameters
param_opti=c("RUE","TTM", "alpha")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(optim_OLS3p<-optim(param_init,  maize.LS.B240.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))
optim_OLS3p
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(optim_OLS3p_a<-optim(param_init,  maize.LS.B240.optim , method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))
optim_OLS3p_a
if(optim_OLS3p_a$value<optim_OLS3p$value){optim_OLS3p<-optim_OLS3p_a}

# third starting point : taking into account values obtained for 2 parameters
param_init <- param_default[param_opti]
param_init[c("RUE","TTM")] <- optim_OLS2p$par
system.time(optim_OLS3p_b<-optim(param_init,  maize.LS.B240.optim , method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=data_n_B240,list_sy=list_n_sy ))
if(optim_OLS3p_b$value<optim_OLS3p$value){optim_OLS3p<-optim_OLS3p_b}
optim_OLS3p
save(optim_OLS3p,file="maize.case_03_optim_OLS3p.rda")

# compute AIC
param_all = param_default
param_all[param_opti] <- optim_OLS3p$par
simobs240=merge(maize.multisy(param_all,list_n_sy, 100, 250),data_n_B240,by=c("sy","day"))
AIC_B_3p=AICf(simobs240$Bobs,simobs240$B,length(param_opti)+1)
goodness.of.fit(simobs240$Bobs,simobs240$B)

data.frame(param_opti=c("RUE", "RUE&TTM", "RUE&TTM&alpha"),AIC=c(AIC_B_1p["AICcomplete"],AIC_B_2p["AICcomplete"],AIC_B_3p["AICcomplete"]))

# end of script