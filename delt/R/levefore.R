levefore<-function(dendat,B,leaf,minlkm=5,seed=1,lambda=0.01,thres=0.5,
sample="bagg",prune="off",
splitscan=0,seedf=1,
scatter=0,
src="c",method="loglik")
{
#B number of bootstrap samples
#leaf number of leaves in the trees to be grown

n<-dim(dendat)[1]
d<-dim(dendat)[2]
suppo<-supp(dendat)+scatter*rep(c(-1,1),d)
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/(n+1)
}

if (sample=="bagg"){
  dendatB<-bootbagg(dendat,seed,scatter) 
}
else{  # scheme=="baggworpl"
  dendatB<-bootworpl(dendat,seed,scatter)
}

tr<-densplit(dendatB,minlkm,suppo,leaf=0,
             method=method,splitscan=splitscan,seedf=seedf)
treeseq<-prunelev(tr,lambda=lambda,n=n)
approleaf<-roundlnum(treeseq$leafs,leaf)
tr<-picktreelev(treeseq,approleaf)

bi<-2
while (bi<=B){
   
   if (sample=="bagg"){
      dendatB<-bootbagg(dendat,seed+bi-1)  
   }
   else{
      dendatB<-bootworpl(dendat,seed+bi-1)
   }

   trnew<-densplit(dendatB,minlkm,suppo,leaf=0,
                   method=method,splitscan=splitscan,seedf=seedf)
   treeseq<-prunelev(trnew,lambda=lambda,n=n)
   approleaf<-roundlnum(treeseq$leafs,leaf)
   trnew<-picktreelev(treeseq,approleaf)

   tr<-treeadd(tr,trnew,bi-1)
   bi<-bi+1
   
}

tr<-c(tr,list(suppo=suppo,step=step))

for (i in 1:length(tr$mean)){
    if (tr$mean[i]>=thres) tr$mean[i]<-1
    else tr$mean[i]<-0
}
#tr$mean<-round(tr$mean)

return(tr)
}







