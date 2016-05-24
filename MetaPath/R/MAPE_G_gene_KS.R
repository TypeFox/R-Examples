MAPE_G_gene_KS <-
function(study,label,censoring.status,DB.matrix,size.min=15,size.max=500,nperm=500,stat=NULL,rth.value=NULL,resp.type){

gene.common=featureNames(study[[1]])
for(t1 in 2:length(study)){
	gene.common=intersect(gene.common,featureNames(study[[t1]]))
}
	
if (is.null(names(study))) names(study)=paste('study.',1:length(study),sep="")

madata=list()
for(t1 in 1:length(study)){
	madata[[t1]]=study[[t1]][gene.common,]
}



Tstat.perm=list()
out=list()
#################################
### meta analysis on genes

for(t1 in 1:length(study)){
	Tstat.perm[[t1]]=list()
	x=exprs(madata[[t1]])
	testlabel=madata[[t1]][[label]]
	
	if(resp.type=="twoclass") {	
		Tstat.perm[[t1]]=Tperm.sample(x=x, fac=as.factor(testlabel), nperm=nperm)
		
   	} else if (resp.type=="multiclass") {	
		Tstat.perm[[t1]]=F.perm.sample(x=x, fac=as.factor(testlabel), nperm=nperm)
		
	  } else  if (resp.type=="survival") {
	     if (is.null(censoring.status)) {
			stop("Error: censoring.status is null")
		  }
		  
		censoring=madata[[t1]][[censoring.status]] 
		Tstat.perm[[t1]]=cox.perm.sample(expr=x, testgroup=testlabel,censoring=censoring,nperm=nperm)
	
   }   else  if (resp.type=="continuous") {
	Tstat.perm[[t1]]=reg.perm.sample(expr=x, testgroup=testlabel,nperm=nperm)
		
  } else {
  	stop("Error: Wrong input augument for resp.type")
  }
 
	out[[t1]]=list()
	out[[t1]]=pqvalues.compute(Tstat.perm[[t1]]$obs,Tstat.perm[[t1]]$perms,Stat.type='Tstat')
}



pvalue.B.array=array(data=NA,dim=c(length(gene.common),nperm,length(study)),dimnames=names(study))
pvalue.0.mtx=matrix(NA,length(gene.common),length(study))
qvalue.0.mtx=matrix(NA,length(gene.common),length(study))
	
for(t1 in 1:length(study)){
pvalue.B.array[,,t1]=out[[t1]]$pvalue.B
pvalue.0.mtx[,t1]=out[[t1]]$pvalue.0
qvalue.0.mtx[,t1]=out[[t1]]$qvalue.0
}
rownames(qvalue.0.mtx)=gene.common
rownames(pvalue.0.mtx)=gene.common

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



## pqvalues calculation
meta.out=pqvalues.compute(P.0,P.B,Stat.type='Pvalue')

gene.qvalues=as.vector(meta.out$qvalue.0)
names(gene.qvalues)=rownames(meta.out$qvalue.0)

gene.pvalues=as.vector(meta.out$pvalue.0)
names(gene.pvalues)=rownames(meta.out$pvalue.0)


## enrichment analysis
meta.set=Enrichment_KS_gene(madata=NULL,label=NULL,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,gene.pvalues=gene.qvalues)


return(list(pvalue.meta=meta.set$pvalue.set.0, qvalue.meta=meta.set$qvalue.set.0,  
			pvalue.meta.B=meta.set$pvalue.set.B ,gene.meta.qvalues=gene.qvalues,gene.meta.pvalues=gene.pvalues,
			gene.study.qvalues=qvalue.0.mtx,gene.study.pvalues=pvalue.0.mtx ))

}
