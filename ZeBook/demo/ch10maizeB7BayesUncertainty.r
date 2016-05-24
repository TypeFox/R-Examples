################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
############################### MAIN PROGRAM ###################################
# Chapter 10. Putting it all together in a case study
library(ZeBook)
library(coda)
library(mnormt)

list_n_sy=unique(maize.data_EuropeEU$sy)

################################################################################
# 1/ definition of the prior distribution of parameter values
# read parameter value : nominal, minimum and maximum
# Uniform
maize.define.param()
sdate=100
ldate=241 # 240 instead of 250, to make the model run faster (the data are at maximum at day 240)

################################################################################
# read distribution of parameters et residual variances
MCMC.R1.R2.garde_ssautocor=read.table("MCMC.R1.R2.garde_ssautocor.dat")

Xparam = MCMC.R1.R2.garde_ssautocor

param.mat=Xparam[,param.opti]
id_param=as.numeric(rownames(param.mat))

################################################################################
# 3/a Uncertainty on one single site-year "18-2006"
idsite = 18
year = 2006
weather <- maize.weather(working.year=year, working.site=idsite,weather_all=weather_EuropeEU)

# only parameters
simUnc240=maize.simule240(param.mat,  weather=weather, sdate=100, ldate=250)
mean(simUnc240)
q90=quantile(simUnc240, prob=c(0.05,0.95))
q90

# with residuals variance
B.epsilon.pred<-matrix(rnorm(length(id_param)*100,mean=0,sd=Xparam[,"B.sigma.new"]),length(id_param),100)
#plot(Xparam[1:100,"B.sigma.new"], apply(B.epsilon.pred[1:100,],1,sd))

B240.pred= B.epsilon.pred + as.vector(simUnc240)
B240.pred[B240.pred<0]=0

mean(B240.pred)
q90e=quantile(B240.pred, prob=c(0.05,0.95))
q90e

par(mfrow=c(1,1), mar=c(4,4,1,1))
hist(B240.pred,main="")
abline(v=q90e, col="red", lwd=3,lty=2)
abline(v=mean(B240.pred),  lwd=5)

################################################################################
# 3/b verification with data of credible intervals

data240 = maize.data_EuropeEU[maize.data_EuropeEU$day==240,]

set.seed(1234)
uncertainty.sy=function(sy) {
simUnc240=maize.simule240(param.mat,weather=maize.weather(working.year=strsplit(sy,"-")[[1]][2],working.site=strsplit(sy,"-")[[1]][1],weather_all=weather_EuropeEU), sdate=100, ldate=250);
B.epsilon.pred<-matrix(rnorm(dim(param.mat)[1]*100,mean=0,sd=Xparam[,"B.sigma.new"]),dim(param.mat)[1],100);
B240.pred= B.epsilon.pred + as.vector(simUnc240);
B240.pred[B240.pred<0]=0;
return(data.frame(sy=sy, mean=mean(B240.pred),q05=quantile(B240.pred, prob=0.05),q10=quantile(B240.pred, prob=0.10),q25=quantile(B240.pred, prob=0.25),q75=quantile(B240.pred, prob=0.75),q90=quantile(B240.pred, prob=0.90), q95=quantile(B240.pred, prob=0.95),Bobs240=data240$Bobs[data240$sy==sy],  stringsAsFactors =FALSE)) ;
}

system.time(sumary.uncertainty.all.sy<-lapply(list_n_sy,uncertainty.sy ))
sumary.uncertainty.all.sy=do.call(rbind,sumary.uncertainty.all.sy)
row.names(sumary.uncertainty.all.sy)=NULL

mean((sumary.uncertainty.all.sy$Bobs240>=sumary.uncertainty.all.sy$q05)&(sumary.uncertainty.all.sy$Bobs240<=sumary.uncertainty.all.sy$q95))
mean((sumary.uncertainty.all.sy$Bobs240>=sumary.uncertainty.all.sy$q10)&(sumary.uncertainty.all.sy$Bobs240<=sumary.uncertainty.all.sy$q90))
mean((sumary.uncertainty.all.sy$Bobs240>=sumary.uncertainty.all.sy$q25)&(sumary.uncertainty.all.sy$Bobs240<=sumary.uncertainty.all.sy$q75))


# end of file