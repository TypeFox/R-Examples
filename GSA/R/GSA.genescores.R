GSA.genescores=function(geneset.number, genesets,  GSA.obj,  genenames, negfirst=FALSE){
 
  tt=GSA.obj$gene.scores

  oo=match(genesets[[geneset.number]],genenames)
oo=oo[!is.na(oo)]
ooo=order(tt[oo],decreasing=!negfirst)
  res=cbind(genenames[oo], round(tt[oo],3))[ooo,,drop=F]

      dimnames(res)=list(NULL,c("Gene","Score"))
return(res)
}

