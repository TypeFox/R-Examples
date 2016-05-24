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

for (p in .libPaths()){try(source(paste(p,"/ZeBook/demo/","ch10maizeB7Bayes_functionsMH.r",sep="")),silent=TRUE)}

list_n_sy=unique(maize.data_EuropeEU$sy)

################################################################################
# 1/ definition of the prior distribution of parameter values
# read parameter value : nominal, minimum and maximum
# Uniform
maize.define.param()
sdate=100
ldate=241 # 240 instead of 250, to make the model run faster (the data are at maximum at day 240)

################################################################################
# 2/ Bayesian estimation of parameter  : MCMC
# A/ RUE

param.opti=c("RUE")
# variance
LAI.sigma<-sqrt(10)
B.sigma<-sqrt(225000)
n_iter=500
# chaine 1.
set.seed(123)
coeff.teta0<-0.25
MetropolisHastings_Gibbs(param.apriori=maize.define.param(),param.opti=param.opti,Nb.iterations=n_iter,LAI.sigma=1,B.sigma=10,coeff.teta0=coeff.teta0,list_sy=list_n_sy,data=maize.data_EuropeEU,sdate,ldate, NomFichierSortie="mcmc_1")
# chaine 2.
set.seed(321)
coeff.teta0<-0.75
MetropolisHastings_Gibbs(param.apriori=maize.define.param(),param.opti=param.opti,Nb.iterations=n_iter,LAI.sigma=1,B.sigma=10,coeff.teta0=coeff.teta0,list_sy=list_n_sy,data=maize.data_EuropeEU,sdate,ldate, NomFichierSortie="mcmc_2")
############## analyse CHAINES

NomFichierSortie<-c("mcmc_2.csv","mcmc_1.csv")
l<-Lecture.chaines(NomChaines=NomFichierSortie)
MCMC.result1<-l$MCMC.result1
MCMC.result2<-l$MCMC.result2
N = dim(MCMC.result1)[1]

### taux de rejet
  Sequence<-data.frame(chaine1=MCMC.result1$rmq,chaine2=MCMC.result2$rmq)
  compt<-apply(Sequence,2,table)
  nb_rejet1<-compt[2,]  #cas 1 pas ds l'apriori
  nb_rejet2<-compt[4,]  # cas 2.21
  nb_rejet<-nb_rejet1+nb_rejet2
  cas1_temp<- nb_rejet1/nb_rejet *100
  taux_rejet_temp<- nb_rejet/length(MCMC.result1$rmq) *100


### Resultats de convergence, correlations
MCMC.R1<-mcmc(MCMC.result1[,param.opti])
MCMC.R2<-mcmc(MCMC.result2[,param.opti])
MCMC.L<-mcmc.list(MCMC.R1,MCMC.R2)
gelman.diag(MCMC.L, confidence = 0.95, transform=FALSE, autoburnin=TRUE)
plot(MCMC.L)
MCMC.R1.garde<-mcmc(MCMC.result1[floor(N/2):N,c(param.opti,"LAI.MSE","B.MSE","LAI.sigma.new", "B.sigma.new")])
MCMC.R2.garde<-mcmc(MCMC.result2[floor(N/2):N,c(param.opti,"LAI.MSE","B.MSE","LAI.sigma.new", "B.sigma.new")])
MCMC.L.garde<-mcmc.list(MCMC.R1.garde,MCMC.R2.garde)
summary(MCMC.L.garde)
MCMC.R1.R2.garde=rbind((MCMC.result1[floor(N/2):N,c(param.opti,"LAI.MSE","B.MSE","LAI.sigma.new", "B.sigma.new")]),(MCMC.result2[floor(N/2):N,c(param.opti,"LAI.MSE","B.MSE","LAI.sigma.new", "B.sigma.new")]))
autocorr.plot(MCMC.R1.R2.garde,lag.max=1000,ask=TRUE)
#autocorr(mcmc(MCMC.R1.R2.garde), lags = c(1, 10, 50, 100, 200,250,300,350,400,500,1000,5000,10000), relative=TRUE)

effectiveSize(MCMC.L)
#geweke.diag(MCMC.R1,MCMC.R2)

# end of file