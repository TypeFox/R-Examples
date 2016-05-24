
GSA.xl.plot=function(GSA.obj,  fac=.10, FDRcut=1){

eps=.001
a=GSA.listsets(GSA.obj, geneset.names=NULL,  FDRcut=FDRcut)
neg=matrix(as.numeric(as.character(a$negative[,-2])),ncol=4)+eps
pos=matrix(as.numeric(as.character(a$positive[,-2])),ncol=4)+eps

ymax=max(c(neg[,4],pos[,4]))
xmax=max(c(neg[,3],pos[,3]))
ymin=min(c(neg[,4],pos[,4]))
xmin=min(c(neg[,3],pos[,3]))

o1=order(neg[,3])
o2=order(pos[,3])
amount1=runif(nrow(neg),min=-fac*neg[o1,3],max=fac*neg[o1,3])
amount2=runif(nrow(pos),min=-fac*pos[o2,3],max=fac*pos[o2,3])

return(list(xneg=neg[o1,3]+amount1,yneg=neg[o1,4],xpos=pos[o2,3]+amount2,ypos=pos[o2,4]))
}


GSA.xl.summary.genesets=function(GSA.genesets.obj){
  nsets=length(GSA.genesets.obj$genesets)
  ngenes=unlist(lapply(GSA.genesets.obj$genesets,length))
  allgenes=unlist(GSA.genesets.obj$genesets)
  tt=table(ngenes,dnn=NULL)
if(length(tt)>20){
qq=round(quantile(ngenes,(0:20)/20),0)
  tt=table(cut(ngenes,breaks=unique(qq)),dnn=NULL)
}

tab.ngenes=matrix(tt, nrow=1)
dimnames(tab.ngenes)[[2]]=names(tt)
  return(list(nsets=nsets, totgenes=sum(ngenes), totunique.genes=length(unique(allgenes)), tab.ngenes=tab.ngenes))
}


GSA.xl.correlate=function(GSA.genesets.obj, genenames){
  nsets=length(GSA.genesets.obj$genesets)
  ngenes=unlist(lapply(GSA.genesets.obj$genesets,length))
  allgenes=unlist(GSA.genesets.obj$genesets)
  sets.in.exp=match(unique(allgenes),genenames)
  exp.in.sets=match(genenames,allgenes)
  

  nn=rep(NA,nsets)
  for(i in 1:nsets){
    nn[i]=sum(!is.na(match(GSA.genesets.obj$genesets[[i]],genenames)))
  }
qq=quantile(nn/ngenes, seq(0,1,by=.1))
quant=matrix(qq,nrow=1)
dimnames(quant)[[2]]=names(qq)

  return(list(totgenes=length(genenames),tot.unique.genes=length(unique(genenames)),
 num.both=sum(!is.na(sets.in.exp)), quant=quant))
}

