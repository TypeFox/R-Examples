## analyse data from JAGS
library(coda);
sim.dat1<-read.coda("outchain1.txt","outindex.txt");
sim.dat2<-read.coda("outchain2.txt","outindex.txt");
Loc.x<-mcmc.list(sim.dat1[,"Loc.x"],sim.dat2[,"Loc.x"]);
plot(Loc.x);
Loc.y<-mcmc.list(sim.dat1[,"Loc.y"],sim.dat2[,"Loc.y"]);
plot(Loc.y);
Year<-mcmc.list(sim.dat1[,"Year"],sim.dat2[,"Year"]);
plot(Year);
# convergence seems fine.

#check lag for residual correlation - want independent realisations
par(mfrow=c(2,2));
acf(sim.dat1[,"Loc.x"]);
acf(sim.dat1[,"Loc.y"]);
acf(sim.dat1[,"Year"]);#indep samples
#seems fine

#####################################################################################
############ Generated data seems fine so now create a new data set
#####################################################################################
#load original pigs data
library(abn);#provides pigs.1par
#we want only the same number of samples as in pigs dataset
get.these<-sample(x=1:dim(sim.dat1)[1],size=dim(pigs.1par)[1],replace=FALSE);

pigs.boot<-sim.dat1[get.these,];
pigs.boot<-as.data.frame(pigs.boot);
#now coerce to factors if need be and set levels - NOTES setting levels works as
# "0" "1" is in the same order as "absent" "present" from original data
for(i in 1:dim(pigs.1par)[2]){if(is.factor(pigs.1par[,i])){pigs.boot[,i]<-as.factor(pigs.boot[,i]);
                                                 levels(pigs.boot[,i])<-levels(pigs.1par[,i]);}}
# pigs.boot is a SINGLE bootstrap data set
save(pigs.boot,file="pigs_boot1.RData",compress=TRUE);


