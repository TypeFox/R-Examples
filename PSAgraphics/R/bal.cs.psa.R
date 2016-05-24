bal.cs.psa<-function(categorical, treatment=NULL, strata=NULL, B=1000, 
    eps = .02, main=NULL,...){

#This function provides an ad-hoc measure of the balance achieved between control and treatment 
#groups for a categorical variable from user defined strata. This is compared to the same measure 
#for randomly generated strata. The proportion of each category for a given strata and treatment 
#is calculated. The proportions are transformed slightly so as to give more weight to differences in 
#small categories.  Then the absolute difference between (transformed) treatment proportions calculated for each 
#category and stratum is found and summed, both within and across strata.  This constitutes 
#the basic statistic.  For comparison purposes, strata are randomly generated 'B' times and the same statistic
#recorded for the random strata.  The percentile rank of the true measure is found in comparison with the randomly 
#generated distribution.  A histogram of the distribution and with the true value highlighted is generated.
#Output is a list with the true measure, the percentile of the true measure, and the 'B' random replicates. 

#If "categorical" has three columns, treat as c, t, s.
if(dim(as.data.frame(categorical))[2]==3){ treatment   <- categorical[,2]
                                           strata      <- categorical[,3]
                                           categorical <- categorical[,1]}

n<-length(categorical)
nstrat<-dim(table(strata))

#Calculation of the statistic of the original strata
sum.abs.diff.original<-NULL
table.cts.original<-table(categorical,treatment,strata)
sum.strata.01<-NULL
for(j in 1:nstrat){sum.strata.01<-rbind(sum.strata.01,apply(table.cts.original[,,j],2,sum))}

prop.bt<-NULL
  for(j in 1:nstrat){
    # prop.i is a vector of the proportion of 0/1 by categorical level
    # level for stratum j
        prop.0 <- table.cts.original[,1,j]/sum.strata.01[j,1]
        prop.1 <- table.cts.original[,2,j]/sum.strata.01[j,2]
        prop.bt <- sum(prop.bt,sum(abs((prop.0 + eps)^.5 - (prop.1 + eps)^.5)))  
                    }
sum.abs.diff.original<-prop.bt              
                   
sum.abs.diff<-NULL
for(i in 1:B){

rand.strata<-sample(strata)
table.cts<-table(categorical,treatment,rand.strata)
sum.strata.01<-NULL
for(j in 1:nstrat){sum.strata.01<-rbind(sum.strata.01,apply(table.cts[,,j],2,sum))}

prop.bt<-NULL
  for(j in 1:nstrat){
    # prop.i is a vector of the proportion of 0/1 by categorical level
    # level for stratum j
        prop.0 <- table.cts[,1,j]/sum.strata.01[j,1]
        prop.1 <- table.cts[,2,j]/sum.strata.01[j,2]
        prop.bt <- sum(prop.bt,sum(abs((prop.0 + .02)^.5 - (prop.1 + .02)^.5)))
        }
sum.abs.diff<-c(sum.abs.diff,prop.bt)              
                  }
rnk<-rank(c(sum.abs.diff.original,sum.abs.diff))[1]
              
hist(c(sum.abs.diff.original,sum.abs.diff),xlab=paste("Balance of",B,"Randomly Generated Strata"),main=main,...)
points(sum.abs.diff.original,0,col="red",cex=4,pch=20)
legend(x="topleft", legend=list(paste("Original Balance:",round(sum.abs.diff.original,2)),
           paste("Rank of Original:",round(rnk,2),"of",B+1)),pch=c(19,19),col=c(2,0),bty="n")
out<-list(balance.orig=sum.abs.diff.original,rank.orig=rnk)
return(out)              }

