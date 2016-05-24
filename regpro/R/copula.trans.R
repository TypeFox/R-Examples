copula.trans<-function(dendat,marginal=rep("gauss",dim(dendat)[2]),remna=TRUE)
{

n<-dim(dendat)[1]
d<-dim(dendat)[2]
copdat<-dendat

for (ii in 1:d){
   if ((marginal[ii]=="gauss")||(marginal[ii]=="uniform")){
      or<-order(dendat[,ii])
      mones<-matrix(0,n,1)
      for (i in 1:n) mones[or[i]]<-i  
      copdat[,ii]<-mones/(n+1)  # copdat[or,ii]<-seq(1:n)/(n+1)
   }
   if (marginal[ii]=="gauss"){ 
      copdat[,ii]<-qnorm(copdat[,ii])
      for (jj in 1:d){
           copdat[(copdat[,jj]==Inf),jj]<-NA
           copdat[(copdat[,jj]==-Inf),jj]<-NA
      }
   }
}

# we remove those rows where there is at least one NA
if (remna){
   for (ii in 1:d){
      indexes<-!is.na(copdat[,ii])
      copdat<-matrix(copdat[indexes,],sum(indexes),d)
   }
}
return(copdat)
}


