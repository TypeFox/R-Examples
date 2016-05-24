MAPE_P_sample_KS <-
function(study,label,censoring.status=NULL,DB.matrix,size.min=15,size.max=500,nperm=500,stat,rth.value=NULL,resp.type){
	
if (is.null(names(study))) names(study)=paste('study.',1:length(study),sep="")

out=list()
for(t1 in 1:length(study)){
	madata=study[[t1]]
	testlabel=madata[[label]]
	out[[t1]]=list()
	
	if (resp.type=="survival") {
		censoring=madata[[censoring.status]]
	}
	out[[t1]]=Enrichment_KS_sample(madata=madata,label=testlabel,censoring=censoring,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,resp.type=resp.type)

}

common.pathway=rownames(out[[1]]$pvalue.set.0)
for(t1 in 1:length(study)) {
	common.pathway=intersect(common.pathway,rownames(out[[t1]]$pvalue.set.0))
}

pvalue.B.array=array(data=NA,dim=c(length(common.pathway),nperm,length(study)))
pvalue.0.mtx=matrix(NA,length(common.pathway),length(study))
qvalue.0.mtx=matrix(NA,length(common.pathway),length(study))
rownames(pvalue.0.mtx)=common.pathway
colnames(pvalue.0.mtx)=names(study)
rownames(qvalue.0.mtx)=common.pathway
colnames(qvalue.0.mtx)=names(study)
dimnames(pvalue.B.array)=list(common.pathway,paste('perm',1:nperm,sep=''),names(study))

for(t1 in 1:length(study)){
	pvalue.B.array[,,t1]=out[[t1]]$pvalue.set.B[common.pathway,]
	pvalue.0.mtx[,t1]=out[[t1]]$pvalue.set.0[common.pathway,]
	qvalue.0.mtx[,t1]=out[[t1]]$qvalue.set.0[common.pathway,]
}


## statistics for meta-analysis
if(stat=='maxP'){
	## maxP statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,max))
	rownames(P.0)=rownames(pvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),max)
	rownames(P.B)=rownames(pvalue.0.mtx)
} else if (stat=='minP'){
	## minP statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,min))
	rownames(P.0)=rownames(pvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),min)
	rownames(P.B)=rownames(pvalue.0.mtx)
} else if (stat=='rth'){
	## rth statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,function(x) sort(x)[ceiling(rth.value*ncol(pvalue.0.mtx))]))
	rownames(P.0)=rownames(pvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),function(x) sort(x)[ceiling(rth.value*dim(pvalue.B.array)[3])])
	rownames(P.B)=rownames(pvalue.0.mtx)
} else if (stat=='Fisher'){
	DF=2*length(study)
	## rth statistics
	P.0=as.matrix(apply(pvalue.0.mtx,1,function(x)  pchisq(-2*sum(log(x)),DF,lower.tail=T) ))
	rownames(P.0)=rownames(pvalue.0.mtx)
	P.B=apply(pvalue.B.array,c(1,2),function(x) pchisq(-2*sum(log(x)),DF,lower.tail=T) )
	rownames(P.B)=rownames(pvalue.0.mtx)
}   else { stop("Please check: the selection of stat should be one of the following options: maxP,minP,rth and Fisher") }

colnames(P.0)='perm0'
colnames(P.B)=paste('perm',1:ncol(P.B),sep='')


## pqvalues calculation
meta.out=pqvalues.compute(P.0,P.B,Stat.type='Pvalue')

colnames(meta.out$pvalue.0)='MAPE_P_sample'
colnames(meta.out$qvalue.0)='MAPE_P_sample'

return(list(pvalue.meta=meta.out$pvalue.0,qvalue.meta=meta.out$qvalue.0,pvalue.meta.B=meta.out$pvalue.B,
		pvalue.set.study=pvalue.0.mtx,qvalue.set.study=qvalue.0.mtx))

}
