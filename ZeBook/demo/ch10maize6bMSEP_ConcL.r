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
################################################################################
# 2) MSEP by cross validation Using concentrated likelihood and a log transformation using all data.
list_n_sy=unique(maize.data_EuropeEU$sy)

param_default=maize.define.param()["nominal",]
param_opti=c("RUE")
best_param=param_default
best_param[param_opti]=2.183994
param_init <- best_param[param_opti]

#for each site-year
MSEPcvB240<-data.frame()
MSEPcvB<-data.frame()
MSEPcvLAI<-data.frame()
for (sy in  list_n_sy)
{
    # show site-year left out for computation of componenent of MSEP
    print(paste("site year left out ",paste(sy, collapse = " - ")) )
    list_n_sy_w_sy<- list_n_sy[list_n_sy%in%sy]
    # exclude this site-year for calibration
    list_n_sy_wo_sy<- list_n_sy[!list_n_sy%in%sy]
    #estimate parameter on sites-years without sy
    print(paste("estimate parameter on sites-years :", paste(list_n_sy_wo_sy, collapse = " ; ")))
    system.time(result_opti<-optim(param_init, maize.logConclikelihood.optim, method = "Nelder-Mead",control=list(trace=1, maxit=1000),param_opti=param_opti,param_default=param_default,data=maize.data_EuropeEU,list_sy=list_n_sy_wo_sy,transf=log ))
    # show results of estimation, to see if there are problems with convergence
    print(result_opti)
    best_param_sy<-param_default
    best_param_sy[param_opti]<-result_opti$par

    #simulation for year y and with param optim on other years
    simobs=merge(maize.multisy(best_param_sy,list_n_sy,100,250),maize.data_EuropeEU,by=c("sy","day"))
    simobs240=subset(simobs,day==240)

    # compute componnent of MSEPcv    
    MSEPcvB240_i=goodness.of.fit(simobs240$Bobs,simobs240$B)["MSE"]
    MSEPcvB_i=goodness.of.fit(simobs$Bobs,simobs$B)["MSE"]
    MSEPcvLAI_i=goodness.of.fit(simobs$LAIobs,simobs$LAI)["MSE"]

    MSEPcvB240<-rbind(MSEPcvB240, MSEPcvB240_i)
    MSEPcvB<-rbind(MSEPcvB, MSEPcvB_i)
    MSEPcvLAI<-rbind(MSEPcvLAI, MSEPcvLAI_i)
}

#show mean of MSE values
MSEPcvB240
mean(MSEPcvB240)

MSEPcvB
mean(MSEPcvB)

MSEPcvLAI
mean(MSEPcvLAI)

# end of file