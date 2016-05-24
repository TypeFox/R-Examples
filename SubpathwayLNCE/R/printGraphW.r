printGraphW<-function(ann,detail=FALSE){
	  if(detail==FALSE){
	  pathwayId<-sapply(ann,function(x) x$pathwayId)
      pathwayName<-sapply(ann,function(x) x$pathwayName)
      annMoleculeRatio<-sapply(ann,function(x) paste(x$annMoleculeNumber,x$moleculeNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
	  weight<-sapply(ann,function(x) x$weight)
      #qvalue<-sapply(ann,function(x) x$qvalue)
	  #lfdr<-sapply(ann,function(x) x$lfdr)
	  fdr<-sapply(ann,function(x) x$fdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annMoleculeRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annMoleculeRatio=annMoleculeRatio,
	  annBgRatio=annBgRatio,annWeight=weight,pvalue=pvalue,fdr=fdr,stringsAsFactors=FALSE)							 
	  }
	  else{	 
      pathwayId<-sapply(ann,function(x) x$pathwayId)	  
	  pathwayName<-sapply(ann,function(x) x$pathwayName)
	  annMoleculeList<-sapply(ann, function(x){ paste(x$annMoleculeList,collapse=";") })
      annBgMoleculeList<-sapply(ann, function(x){ paste(x$annBgMoleculeList,collapse=";")})
	  annMoleculeRatio<-sapply(ann,function(x) paste(x$annMoleculeNumber,x$moleculeNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
	  weight<-sapply(ann,function(x) x$weight)
      #qvalue<-sapply(ann,function(x) x$qvalue)
	  #lfdr<-sapply(ann,function(x) x$lfdr)
	  fdr<-sapply(ann,function(x) x$fdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annMoleculeRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr,annMoleculeList,annBgMoleculeList))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annMoleculeRatio=annMoleculeRatio,
	  annBgRatio=annBgRatio,annWeight=weight,pvalue=pvalue,fdr=fdr,annMoleculeList=annMoleculeList,
	  annBgMoleculeList=annBgMoleculeList,stringsAsFactors=FALSE)								 
	  }
      return(ann.data.frame)
}