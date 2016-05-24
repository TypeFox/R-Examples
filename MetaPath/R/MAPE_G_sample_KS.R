MAPE_G_sample_KS <-
function(study,label,censoring.status=NULL,DB.matrix,size.min=15,size.max=500,nperm=500,stat,rth.value=NULL,resp.type){

gene.common=featureNames(study[[1]])
for(t1 in 2:length(study)){
	gene.common=intersect(gene.common,featureNames(study[[t1]]))
}
	
if (is.null(names(study))) names(study)=paste('study.',1:length(study),sep="")

madata=list()
for(t1 in 1:length(study)){
	madata[[t1]]=study[[t1]][gene.common,]
}


######################

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
		
	 }	else  if (resp.type=="survival") {
	     if (is.null(censoring.status)) {
			stop("Error: censoring.aus is null")
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

colnames(P.0)='perm0'
colnames(P.B)=paste('perm',1:ncol(P.B),sep='')


## pqvalues calculation
meta.out=pqvalues.compute(P.0,P.B,Stat.type='Pvalue')

gene.qvalues=as.vector(meta.out$qvalue.0)
names(gene.qvalues)=rownames(meta.out$qvalue.0)

gene.pvalues=as.vector(meta.out$pvalue.0)
names(gene.pvalues)=rownames(meta.out$pvalue.0)


##################  Enrichment

score.0=P.0
score.B=P.B

genes.in.study=rownames(score.0)
set2allgenes.mtx=DB.matrix
gene.common=intersect(genes.in.study,colnames(set2allgenes.mtx))
set2allgenes.mtx=set2allgenes.mtx[,gene.common,drop=F]
score.0=score.0[gene.common,,drop=F]
score.B=score.B[gene.common,,drop=F]


## select the gene sets whose size >=size.min and <=size.max
set.length=apply(set2allgenes.mtx,1,sum)
idx.1=which(set.length>=size.min)
idx.2=which(set.length<=size.max)
set.idx=intersect(idx.1,idx.2)
if (length(set.idx)<1) {
	stop('no gene sets satisfying size.min<=size<=size.max')
	} else {
	set2allgenes.mtx=set2allgenes.mtx[set.idx,,drop=F]
} 


gene.name.sort=names(sort(abs(score.0[,1]), decreasing = F))  ###???
	

## KS test
## 3 order matrix

order.mtx.1=set2allgenes.mtx[,gene.name.sort,drop=F]
order.mtx.0=( 1-order.mtx.1 )
order.mtx.1=t(apply(order.mtx.1,1,function(x) x/sum(x)))
order.mtx.0=-t(apply(order.mtx.0,1,function(x) x/sum(x)))
order.mtx=order.mtx.0+order.mtx.1


## 4 ES score 
ES.0=as.matrix(apply(t(apply(order.mtx,1,cumsum)),1,max))

## 5 permute gene labels  
ES.B=matrix(NA,nrow(ES.0),ncol(score.B))
for(t1 in 1:ncol(score.B)){
	gene.name.sort=names(sort(abs(score.B[,t1]), decreasing = F))  ###???
	order.mtx.1=set2allgenes.mtx[,gene.name.sort,drop=F]
	order.mtx.0=( 1-order.mtx.1 )
	order.mtx.1=t(apply(order.mtx.1,1,function(x) x/sum(x)))
	order.mtx.0=-t(apply(order.mtx.0,1,function(x) x/sum(x)))
	order.mtx=order.mtx.0+order.mtx.1
	order.cumsum=t(apply(order.mtx,1,cumsum))
	ES.B[,t1]=apply(order.cumsum,1,max)
}

rownames(ES.B)=rownames(order.mtx)


### 6 pvalues calculation
N.X=apply(set2allgenes.mtx,1,sum)
N.Y=ncol(set2allgenes.mtx)-N.X
N=N.X* N.Y/(N.X + N.Y)

ES.pval.0=exp(-2 * N * ES.0^2)
ES.pval.B=exp(-2 * N * ES.B^2)


## 7 qvalues and qvalues calculation
enrich.out=pqvalues.compute(ES.pval.0,ES.pval.B,Stat.type='Pvalue')

colnames(enrich.out$pvalue.0)='MAPE_G_sample'
colnames(enrich.out$qvalue.0)='MAPE_G_sample'
colnames(enrich.out$pvalue.B)=paste('perm',1:ncol(enrich.out$pvalue.B),sep='')


return(list(pvalue.meta=enrich.out$pvalue.0,qvalue.meta=enrich.out$qvalue.0,pvalue.meta.B=enrich.out$pvalue.B))

}
