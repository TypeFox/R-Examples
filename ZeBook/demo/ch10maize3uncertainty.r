################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-14
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
list_n_sy=unique(maize.data_EuropeEU$sy)

################################################################################
# a function to run the model on a matrix of parameter (virtual design) and give Biomass at day240 derived from maize.simule (ZeBook package)
maize.simule240<-function(X, weather, sdate, ldate, all=FALSE){
# output : biomass at day=240
Y <- apply(X,1,function(v) maize.model2(v[1:7],weather, sdate, ldate)[240-sdate+1,"B"])
if(all) Y = cbind(X,B = Y)
return(as.matrix(Y))
}

# 1/ definition of the distribution of parameter values
#read parameter value : nominal, minimum and maximum
maize.define.param()

# define distribution as independent uniform distribution and generate a random sample
set.seed(123)
N = 10000
param.mat=param.runif(maize.define.param(), N)

# 2/ Uncertainty on one single site-year as an illustration
sy="18-2006"
weather = maize.weather(working.year=strsplit(sy,"-")[[1]][2], working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_EuropeEU)
    
# run the model for all the values of parameters
system.time(simUnc240<-maize.simule240(param.mat,  weather=weather, sdate=100, ldate=250))

# analyse the result
mean(simUnc240)

q90=quantile(simUnc240, prob=c(0.05,0.95))
q90

hist(simUnc240, main=paste("site-year:",sy))
abline(v=q90, col="red", lwd=3,lty=2)
abline(v=mean(simUnc240),  lwd=4)

# 3/ Uncertainty all site-years
system.time(sumary.uncertainty.all.sy<-sapply(list_n_sy, function(sy) {
simUnc240=maize.simule240(param.mat,weather=maize.weather(working.year=strsplit(sy,"-")[[1]][2],working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_EuropeEU), sdate=100, ldate=250)
return(data.frame(sy=sy, mean=mean(simUnc240),q05=quantile(simUnc240, prob=0.05), q95=quantile(simUnc240, prob=0.95),  stringsAsFactors =FALSE))
}))
colnames(sumary.uncertainty.all.sy)<-NULL
t(sumary.uncertainty.all.sy)

# end of file