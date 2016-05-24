MAPE_P_gene_KS <-
function(study,label,censoring.status,DB.matrix,size.min=15,size.max=500,nperm=500,stat=NULL,rth.value=NULL,resp.type){

if (is.null(names(study))) names(study)=paste('study.',1:length(study),sep="")

out=list()
for(t1 in 1:length(study)){
	madata=study[[t1]]
	testlabel=madata[[label]]
	out[[t1]]=list()
	
	if (resp.type=="survival") {
		censoring=madata[[censoring.status]]
	}
	out[[t1]]=Enrichment_KS_gene(madata=madata,label=testlabel,censoring=censoring,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,resp.type=resp.type)

}



set.common=rownames(out[[1]]$pvalue.set.0)
for(t1 in 2:length(study)){
	set.common=intersect(set.common,rownames(out[[t1]]$pvalue.set.0))
}
	
if (is.null(names(study))) names(study)=paste('study.',1:length(study),sep="")
pvalue.B.array=array(data=NA,dim=c(length(set.common),nperm,length(study)),dimnames=names(study))
pvalue.0.mtx=matrix(NA,length(set.common),length(study))
qvalue.0.mtx=matrix(NA,length(set.common),length(study))
	
for(t1 in 1:length(study)){
pvalue.B.array[,,t1]=out[[t1]]$pvalue.set.B[set.common,]
pvalue.0.mtx[,t1]=out[[t1]]$pvalue.set.0[set.common,]
qvalue.0.mtx[,t1]=out[[t1]]$qvalue.set.0[set.common,]
}
rownames(qvalue.0.mtx)=set.common
rownames(pvalue.0.mtx)=set.common

rm(out)

## statistics for meta-analysis
if(stat=='maxP'){
	## maxP statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,max))
	rownames(P.0)=rownames(qvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),max)
	rownames(P.B)=rownames(qvalue.0.mtx)
} else if (stat=='minP'){
	## minP statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,min))
	rownames(P.0)=rownames(qvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),min)
	rownames(P.B)=rownames(qvalue.0.mtx)
} else if (stat=='rth'){
	## rth statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,function(x) sort(x)[ceiling(rth.value*ncol(pvalue.0.mtx))]))
	rownames(P.0)=rownames(qvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),function(x) sort(x)[ceiling(rth.value*dim(pvalue.B.array)[3])])
	rownames(P.B)=rownames(qvalue.0.mtx)
} else if (stat=='Fisher'){
	DF=2*length(study)
	## rth statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,function(x)  pchisq(-2*sum(log(x)),DF,lower.tail=T) ))
	rownames(P.0)=rownames(qvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),function(x) pchisq(-2*sum(log(x)),DF,lower.tail=T) )
	rownames(P.B)=rownames(qvalue.0.mtx)
}   else { stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher") }


## pvalues and qvalues calculation
meta.out=pqvalues.compute(P.0,P.B,Stat.type='Pvalue')

return(list(pvalue.meta=meta.out$pvalue.0, qvalue.meta=meta.out$qvalue.0,  pvalue.meta.B=meta.out$pvalue.B,qvalue.set.study=qvalue.0.mtx,pvalue.set.study=pvalue.0.mtx))

}
