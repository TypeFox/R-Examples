bootworpl<-function(dendat,seed,scatter=0)
{
set.seed(seed)
n<-dim(dendat)[1]
d<-dim(dendat)[2]
nout<-ceiling(n/2)
dendatout<-matrix(0,nout,d)

blacklist<-matrix(0,n,1)
found<-0

while (found<nout){
   res<-ceiling(runif(1)*n)
   if (blacklist[res]==0){
       found<-found+1
       addi<-scatter*2*(runif(d)-0.5)
       dendatout[found,]<-dendat[res,]+addi
       blacklist[res]<-1
   }

}

return(dendatout)
}
