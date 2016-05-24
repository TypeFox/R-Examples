Enrichment_KS_sample <-
function(madata,label,censoring=NULL,DB.matrix,size.min=15,size.max=500,nperm=500,resp.type=NULL){

genes.in.study=featureNames(madata)
set2allgenes.mtx=DB.matrix
gene.common=intersect(featureNames(madata),colnames(set2allgenes.mtx))
set2allgenes.mtx=as.matrix(set2allgenes.mtx[,gene.common],drop=F)
madata=madata[gene.common,]
	
## select the gene sets whose size >=15 and <=500
set.length=apply(set2allgenes.mtx,1,sum)
idx.1=which(set.length>=size.min)
idx.2=which(set.length<=size.max)
set.idx=intersect(idx.1,idx.2)
if (length(set.idx)<1) {
	stop('no gene sets satisfying size.min<=size<=size.max')
	} else {
	set2allgenes.mtx=set2allgenes.mtx[set.idx,,drop=F]
} 
	
## permutation	
	x=exprs(madata)
	testlabel=label
	
	if(resp.type=="twoclass") {	
		Tstat.perm=Tperm.sample(x=x, fac=as.factor(testlabel), nperm=nperm)

   	} else if(resp.type=="multiclass") {	
		Tstat.perm=F.perm.sample(x=x, fac=as.factor(testlabel), nperm=nperm)

   	} else  if (resp.type=="survival") {
	     if (is.null(censoring)) {
			stop("Error: censoring.status is null")
		  }
		  
	Tstat.perm=cox.perm.sample(expr=x, testgroup=testlabel,censoring=censoring,nperm=nperm)
	
   }   else  if (resp.type=="continuous") {
	Tstat.perm=reg.perm.sample(expr=x, testgroup=testlabel,nperm=nperm)
		
  } else {
  	stop("Error: Wrong input augument for resp.type")
  }
 

score.0=Tstat.perm$obs
score.B=Tstat.perm$perms
gene.name.sort=names(sort(abs(score.0), decreasing = T))  

## KS test
## 3 order matrix

order.mtx.1=set2allgenes.mtx[,gene.name.sort,drop=F]
order.mtx.0=( 1-order.mtx.1 )
order.mtx.1=t(apply(order.mtx.1,1,function(x) x/sum(x)))
order.mtx.0=-t(apply(order.mtx.0,1,function(x) x/sum(x)))
order.mtx=order.mtx.0+order.mtx.1

## 4 ES score 
ES.0=as.matrix(apply(t(apply(order.mtx,1,cumsum)),1,max))

## 5 ES score by permutation 
ES.B=matrix(NA,nrow(ES.0),ncol(score.B))
for(t1 in 1:ncol(score.B)){
	gene.name.sort=names(sort(abs(score.B[,t1]), decreasing = T))  
	order.mtx.1=set2allgenes.mtx[,gene.name.sort,drop=F]
	order.mtx.0=( 1-order.mtx.1 )
	order.mtx.1=t(apply(order.mtx.1,1,function(x) x/sum(x)))
	order.mtx.0=-t(apply(order.mtx.0,1,function(x) x/sum(x)))
	order.mtx=order.mtx.0+order.mtx.1
	order.cumsum=t(apply(order.mtx,1,cumsum))
	ES.B[,t1]=apply(order.cumsum,1,max)
}

rownames(ES.B)=rownames(order.mtx)

# 6 pvalues calculation
N.X=apply(set2allgenes.mtx,1,sum)
N.Y=ncol(set2allgenes.mtx)-N.X
N=N.X* N.Y/(N.X + N.Y)
 

## 7 qvalues and qvalues calculation

enrich.out=pqvalues.compute(ES.0,ES.B,Stat.type='Tstat')  


return(list(pvalue.set.0=enrich.out$pvalue.0,pvalue.set.B=enrich.out$pvalue.B,qvalue.set.0=enrich.out$qvalue.0))


}
