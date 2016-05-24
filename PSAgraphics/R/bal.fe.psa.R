bal.fe.psa<-function(categorical, treatment=NULL, strata=NULL, FB=2000){

#This function makes repeated calls to fisher.test, Fisher's Exact test,
#a test of whether the distribution of categorical is 
#independent of treatment within each stratum.
#p-values for the test for each stratum are returned.

#If "categorical" has three columns, treat as c, t, s.
if(dim(as.data.frame(categorical))[2]==3){ treatment   <- categorical[,2]
                                           strata      <- categorical[,3]
                                           categorical <- categorical[,1]}

if(is.factor(treatment)){treatment<-factor(treatment,levels=levels(treatment),labels=c(0,1))}
categorical<-as.factor(categorical)

n<-length(categorical)
nstrat<-dim(table(strata))



contingency<-table(categorical,treatment,strata)
fe<-NULL
for(i in 1:nstrat){ fe<-c(fe,fisher.test(contingency[,,i], B = FB)$p) }

return(fe)              }

