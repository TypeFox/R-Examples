GSA.make.features=function(GSA.func.obj, x, genesets, genenames){
np=length(GSA.func.obj$gene.ind)
xs=t(scale(t(x),center=GSA.func.obj$mean,scale=GSA.func.obj$sd))

val=matrix(NA,nrow=np,ncol=ncol(x))
for(i in 1:np){
if(!is.null(GSA.func.obj$gene.ind[[i]])){
 gene.set=match(genesets[[i]],genenames)
 gene.set=gene.set[!is.na(gene.set)]
 geneind=gene.set[GSA.func.obj$gene.ind[[i]]]
 val[i,]=(colSums(xs[geneind,,drop=F])/length(gene.set))
}
}
return(val)
}

