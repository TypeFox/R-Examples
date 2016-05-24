eval.bagg<-function(dendat,B,leaf,minobs=NULL,seed=1,
sample="bagg",prune="off",
splitscan=0,seedf=1,
scatter=0,
src="c",method="loglik")
{
#B number of bootstrap samples
#leaf number of leaves in the trees to be grown

n<-dim(dendat)[1]
d<-dim(dendat)[2]
if (is.null(minobs)) minobs<-ceiling(sqrt(n)/2)

suppo<-supp(dendat)+scatter*rep(c(-1,1),d)
step<-matrix(0,d,1)
for (i in 1:d){
   step[i]<-(suppo[2*i]-suppo[2*i-1])/(n+1)
}
{
if (sample=="bagg"){
  dendatB<-bootbagg(dendat,seed,scatter) 
}
else{  # scheme=="baggworpl"
  dendatB<-bootworpl(dendat,seed,scatter)
}
}
{
if (prune=="off"){
   if (leaf==0){
        tr<-densplit(dendatB,minobs,leaf=0,
        method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
   }
   else{
        tr<-eval.greedy(dendatB,leaf,method=method,suppo=suppo)
   }
}
else{  #prune == on
  if (src=="c"){
    tr<-densplit(dendatB,minobs,leaf=0,
        method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
    treeseq<-prune(tr)
    approleaf<-roundlnum(treeseq$leafs,leaf)
    tr<-eval.pick(treeseq,approleaf)
  }
  else{
    tr<-densplitF(dendatB,0,minobs,suppo)
    treeseq<-prune(tr)
    approleaf<-roundlnum(treeseq$leafs,leaf)
    tr<-eval.pick(treeseq,approleaf)
  }
}
}

bi<-2
while (bi<=B){
   {
   if (sample=="bagg"){
      dendatB<-bootbagg(dendat,seed+bi-1)  
   }
   else{
      dendatB<-bootworpl(dendat,seed+bi-1)
   }
   }
   if (prune=="off"){
       if (leaf==0){
          trnew<-densplit(dendatB,minobs,leaf=0,
                 method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
       }
       else{
          trnew<-eval.greedy(dendatB,leaf,method=method,suppo=suppo)
       }
   }
   else{  #prune == on
      if (src=="c"){
        trnew<- densplit(dendatB,minobs,leaf=0,
                method=method,splitscan=splitscan,seedf=seedf,suppo=suppo)
      treeseq<-prune(trnew)
      approleaf<-roundlnum(treeseq$leafs,leaf)
      trnew<-eval.pick(treeseq,approleaf)
      }
      else{
      trnew<-densplitF(dendatB,0,minobs,suppo,splitscan=splitscan,seedf=seedf)
      treeseq<-prune(trnew)
      approleaf<-roundlnum(treeseq$leafs,leaf)
      trnew<-eval.pick(treeseq,approleaf)
      }
   }

   tr<-treeadd(tr,trnew,bi-1)
   bi<-bi+1
   
}

tr<-c(tr,list(support=suppo,step=step))
tr$N<-rep(dim(dendatB)[1],d)

tayd<-partigen.disc(tr)
tr$value<-tayd$value
tr$down<-tayd$down
tr$high<-tayd$high

return(tr)
}







