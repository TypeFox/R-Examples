
#################
RNASeq.Data=function(Count,Normalizer=NULL,Treatment,GeneID=NULL){
  Count=data.matrix(Count)
  Count=matrix(Count,ncol=length(Treatment))
  if(is.null(Normalizer)){
     n=Count
     n[n<=0]=NA
     q3=apply(n,2,quantile,.75,na.rm=TRUE)
     Normalizer=log(q3/mean(q3))
    }
  if(is.vector(Normalizer)) Normalizer=matrix(rep(Normalizer,each=nrow(Count)),ncol=length(Treatment))
  Count=Count[,order(Treatment)]
  Normalizer=Normalizer[,order(Treatment)]
  Treatment=Treatment[order(Treatment)]
  Treatment=cumsum(!duplicated(Treatment))
  if(is.null(GeneID)) GeneID=paste(1:nrow(Count))
  colnames(Count)=colnames(Normalizer)=Treatment
  rownames(Count)=GeneID
  if(nrow(Normalizer)==1) Normalizer=c(Normalizer)

  fc=sumRow(Count,by=Treatment)/sumRow(exp(Normalizer),by=Treatment)
  fc[fc<=0]=min(fc[fc>0])/10
  logFC=log(fc)-rowMeans(log(fc))
  colnames(logFC)=unique(Treatment)
  expr=rowMeans(fc)
  disp=est.nb.v.QL(Count,Normalizer,Treatment)

  data=list(GeneID=GeneID,Count=Count,Normalizer=Normalizer,Treatment=Treatment,Aver.Expr=expr,logFC=logFC,NB.Dispersion=disp)
  return(data)
}
