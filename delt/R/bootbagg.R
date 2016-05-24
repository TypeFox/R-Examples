bootbagg<-function(dendat,seed,scatter=0)
{
set.seed(seed)
n<-dim(dendat)[1]
d<-dim(dendat)[2]
dendatout<-matrix(0,n,d)

for (i in 1:n){
   res<-ceiling(runif(1)*n)
   addi<-scatter*2*(runif(d)-0.5)
   dendatout[i,]<-dendat[res,]+addi
}

return(dendatout)
}
