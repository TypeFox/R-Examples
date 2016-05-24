MAPE_I_KS <-
function(MAP_GENE.obj, MAP_SET.obj, study){


set.common=intersect(rownames(MAP_GENE.obj$qvalue.meta),rownames(MAP_SET.obj$qvalue.meta) )

nperm=  ncol(MAP_SET.obj$pvalue.meta.B)

pvalue.set.B.array=array(data=NA,dim=c(length(set.common),nperm,2),dimnames=names(study))
pvalue.set.0.mtx=matrix(NA,length(set.common),2)

pvalue.set.B.array[,,1]=MAP_SET.obj$pvalue.meta.B[set.common,]
pvalue.set.B.array[,,2]=MAP_GENE.obj$pvalue.meta.B[set.common,]

pvalue.set.0.mtx[,1]=MAP_SET.obj$pvalue.meta[set.common,]
pvalue.set.0.mtx[,2]=MAP_GENE.obj$pvalue.meta[set.common,]

rownames(pvalue.set.0.mtx)=set.common

## minP statistics
minP.0=as.matrix(apply(pvalue.set.0.mtx,1,min))
rownames(minP.0)=set.common
minP.B=apply(pvalue.set.B.array,c(1,2),min)
rownames(minP.B)=set.common

## pqvalue calculation
meta.out=pqvalues.compute(minP.0,minP.B,Stat.type='Pvalue')


return(list(pvalue.meta=meta.out$pvalue.0, qvalue.meta=meta.out$qvalue.0))


}
