bal.ms.psa<-function(continuous, treatment = NULL, strata = NULL, trim = 0, B = 1000, 
main = NULL){

#Compares means within randomly generated strata for a continuous covariate.
#The analogue of bal.cs.psa


n<-length(continuous)
nstrat<-dim(table(strata))

#If "continuous" has three columns, treat as c, t, s.
if(dim(as.data.frame(continuous))[2]==3){ treatment   <- continuous[,2]
                                           strata      <- continuous[,3]
                                           continuous <- continuous[,1]}
                                          

meas.mns<-tapply(continuous,list(treatment,strata),mean,trim=trim)
sum.abs.diff.original<-sum((abs(meas.mns[1,]-meas.mns[2,]))*table(strata)/n)

sum.abs.diff<-NULL
for(i in 1:B){
rstrat<-sample(strata,n)
meas.mns<-tapply(continuous,list(treatment,rstrat),mean,trim=trim)
sum.abs.diff<-c(sum.abs.diff,sum((abs(meas.mns[1,]-meas.mns[2,]))*table(strata)/n))
             }

res<-c(sum.abs.diff.original,sum.abs.diff)
rnk<-NULL
rnk<-rank(c(sum.abs.diff.original[1],sum.abs.diff))[1]
hist(c(sum.abs.diff.original,sum.abs.diff),xlab=paste("Balance Statistics for", B, 
"Randomly Permuted Stratifications"), main=main)
points(sum.abs.diff.original,0,col="red",cex=4,pch=20)
legend(x="topright", legend=list(paste("Original Balance:", 
           round(sum.abs.diff.original,2)),
           paste("Rank of Original:",round(rnk,2),"of",B+1)),pch=c(19,19),col=c(2,0))

out<-list(balance.orig=sum.abs.diff.original,rank.orig=rnk)
return(out)         
     }



