################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
################################################################################
maize.ValCalNLS<-function(param,param_opti,param_default,data_sy, list_sy,model=TRUE)
{
  if (length(param)==length(param_opti)){
    param_all = param_default
    param_all[param_opti] <- param
    simobs<-merge(data_sy,maize.multisy(param_all, list_sy,100,250),by=c("sy","day"),sort=FALSE)
    simobs=simobs[order(simobs$sy,simobs$day),]
    # for B
    if (model) return(B=simobs$B)
    else return(Bobs=simobs$Bobs)
   } else { print("Be carreful incoherence in the number of parameters to calibrate")}
}
# maize.ValCalNLS(2.176,c("RUE"), maize.define.param()["nominal",],data_n_B240,list_n_sy)
################################################################################
# Obs function for OLS with utilization of NLS function
################################################################################
list_n_sy=unique(maize.data_EuropeEU$sy)
# 1)  3.5.1.	Using OLS (ordinary least squares), based on the final biomass (day=240) data.
data_n_B240=subset(maize.data_EuropeEU,day==240)
data_n_B240$LAIobs=NA
################################################################################
# 1a) Gradient method (NLS)
nls_cont<-nls.control() #nls.control(maxiter = 150, tol = 1e-05, minFactor = 1/10000,printEval = FALSE, warnOnly = TRUE)
param_default=maize.define.param()["nominal",]
Obs=maize.ValCalNLS(2.176,c("RUE"), maize.define.param()["nominal",],data_n_B240,list_n_sy,model=FALSE)#maize.ValObsNLS(data_n_B240,list_sy=list_n_sy)
# 1 parameter "RUE"
param_opti=c("RUE")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(OLS1p <- try(nls((Obs~maize.ValCalNLS(param,param_opti, param_default,data_n_B240,list_sy=list_n_sy)), start = list(param = param_init), trace = T,control=nls_cont)))
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(OLS1p_b <- try(nls((Obs~maize.ValCalNLS(param,param_opti, param_default,data_n_B240,list_sy=list_n_sy)), start = list(param = param_init), trace = T,control=nls_cont)))
if(summary(OLS1p_b)$sigma<summary(OLS1p)$sigma){OLS1p<-OLS1p_b}

summary(OLS1p)
sum(summary(OLS1p)$residuals^2)
save(OLS1p,file="ch10maize5_OLSgrad_1p.rda")

#param_init = 2.176

# compute AIC and RMSE
param_all = param_default
param_all[param_opti] <- coefficients(OLS1p)
simobs240=merge(maize.multisy(param_all,list_n_sy, 100, 250),data_n_B240,by=c("sy","day"))
AICf(simobs240$Bobs,simobs240$B,length(param_opti)+1)
AIC(OLS1p)
goodness.of.fit(simobs240$Bobs,simobs240$B)

# 2 parameters "RUE","TTM" - will lead to an error  because TTM is a threshold, so the derivatives of the objective functions are discontinuous)
param_opti=c("RUE","TTM")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(OLS2p <- try(nls((Obs~maize.ValCalNLS(param,param_opti, param_default,data_n_B240,list_sy=list_n_sy)), start = list(param = param_init), trace = T,control=nls_cont)))

# 2 parameters "RUE","alpha"
param_opti=c("RUE","alpha")
# first starting point
param_init <- param_default[param_opti]*0.85
system.time(OLS2p <- try(nls((Obs~maize.ValCalNLS(param,param_opti, param_default,data_n_B240,list_sy=list_n_sy)), start = list(param = param_init), trace = T,control=nls_cont)))
# second starting point
param_init <- param_default[param_opti]*1.15
system.time(OLS2p_b <- try(nls((Obs~maize.ValCalNLS(param,param_opti, param_default,data_n_B240,list_sy=list_n_sy)), start = list(param = param_init), trace = T,control=nls_cont)))

if(summary(OLS2p_b)$sigma<summary(OLS2p)$sigma){OLS2p<-OLS2p_b}
summary(OLS2p)
sum(summary(OLS2p)$residuals^2)
save(OLS2p,file="ch10maize5_OLSgrad_2p.rda")

# comparison 1 to 2 param
anova(OLS1p,OLS2p)

# end of script