bal.ks.psa<-function(continuous,treatment=NULL,strata=NULL){

#Computes two sample Kolmogorov-Smirnov statistics and p-values for
#the comparison of the two distrubition of a continuous covariate
#within each defined stratum in a PSA.

#If "continuous" has three columns, treat as m, t, s.
if(dim(as.data.frame(continuous))[2]==3){ treatment   <- continuous[,2]
                                           strata      <- continuous[,3]
                                           continuous <- continuous[,1]}
nstrat<-dim(table(strata))
strat.f<-as.factor(strata)
levels(strat.f)<-1:nstrat

ks<-NULL
for(j in 1:nstrat){
kol.sm<-ks.test(continuous[treatment==unique(treatment)[1]&strat.f==j],continuous[treatment==unique(treatment)[2]&strat.f==j])
ks[j]<-kol.sm$p
}
return(ks)
}
