Enrichment_KS_gene <-
function(madata,label,censoring=NULL,DB.matrix,size.min=15,size.max=500,nperm=500,gene.pvalues=NULL,resp.type=NULL){

##########  madata=NULL;label=NULL;DB.matrix=DB.matrix;size.min=size.min;size.max=size.max;nperm=nperm;gene.pvalues=gene.pvalues
##########  madata=madata;label=testlabel;censoring=NULL;DB.matrix=DB.matrix;size.min=size.min;size.max=size.max;nperm=5;gene.pvalues=NULL;resp.type='continuous'
##########  censoring=sample(c(0,1),length(testlabel),replace=T); resp.type='survival'
############ madata=madata;label=group;censoring=NULL;DB.matrix;size.min=15;size.max=500;nperm=50;gene.pvalues=NULL;resp.type='discrete'

### 1 calculate the pvalues of each genes and ## 2 sort genes according their pvalues#### cor func



if(!is.null(madata)){
	genes.in.study=featureNames(madata)
	set2allgenes.mtx=DB.matrix
	gene.common=intersect(featureNames(madata),colnames(set2allgenes.mtx))
	
	if (nrow(set2allgenes.mtx)>1) {
	set2allgenes.mtx=as.matrix(set2allgenes.mtx[,gene.common])
	} else {
	set2allgenes.mtx=t(as.matrix(set2allgenes.mtx[,gene.common]))
	}
	
	madata=madata[gene.common,]

	
	if(resp.type=="twoclass") {	
		tstat=genefilter::rowttests(exprs(madata), as.factor(label), tstatOnly= F)
		p.values=(tstat$p.value)
		names(p.values)=rownames(tstat)
		gene.name.sort=names(sort(p.values, decreasing = F))
	
	} else if(resp.type=="multiclass") {	
		tstat=genefilter::rowFtests(exprs(madata), as.factor(label), var.equal = TRUE)
		p.values=(tstat$p.value)
		names(p.values)=rownames(tstat)
		gene.name.sort=names(sort(p.values, decreasing = F))
	
	} else  if (resp.type=="survival") {
		 
		 if (is.null(censoring)) {
			stop("Error: censoring status is null")
		  }
	
		cox.out<- coxfunc(exprs(madata), label, censoring)
		scores<-abs(cox.out$tt)
		gene.name.sort=names(sort(scores, decreasing = T))
	
   }   else  if (resp.type=="continuous") {
		cor.out<- cor.func(exprs(madata), label)
		scores<-abs(cor.out$tt[,1])
		gene.name.sort=names(sort(scores, decreasing = T))
	
  } else {
  	stop("Error: Wrong input augument for resp.type")
  }
	
	
} 

 if (!is.null(gene.pvalues)) {

	if( (!is.vector(gene.pvalues)) | (is.null(names(gene.pvalues)))) stop('gene.pvalues should be a vector with gene names')
	genes.in.study=names(gene.pvalues)
	set2allgenes.mtx=DB.matrix
	gene.common=intersect(genes.in.study,colnames(set2allgenes.mtx))
	
	if (nrow(set2allgenes.mtx)>1) {
	set2allgenes.mtx=as.matrix(set2allgenes.mtx[,gene.common])
	} else {
	set2allgenes.mtx=t(as.matrix(set2allgenes.mtx[,gene.common]))
	}
	
	gene.pvalues=gene.pvalues[gene.common]
	gene.name.sort=names(sort(gene.pvalues, decreasing = F))
	
}  

## select the gene sets whose size >=15 and <=500
set.length=apply(set2allgenes.mtx,1,sum)
idx.1=which(set.length>=size.min)
idx.2=which(set.length<=size.max)
set.idx=intersect(idx.1,idx.2)
if (length(set.idx)<1) stop('no gene sets satisfying size.min<=size<=size.max')

if(length(set.idx)>1) {
set2allgenes.mtx=set2allgenes.mtx[set.idx,]
} else {
set2allgenes.mtx=t(as.matrix(set2allgenes.mtx[set.idx,]))
}

## KS test
## 3 order matrix
if (nrow(set2allgenes.mtx)>1) {
order.mtx.1=( set2allgenes.mtx[,gene.name.sort] )
} else {
order.mtx.1=t(as.matrix( set2allgenes.mtx[,gene.name.sort] ))
}
order.mtx.0=( 1-order.mtx.1 )
order.mtx.1=t(apply(order.mtx.1,1,function(x) x/sum(x)))
order.mtx.0=-t(apply(order.mtx.0,1,function(x) x/sum(x)))

order.mtx=order.mtx.0+order.mtx.1


## 4 ES score 
ES.0=as.matrix(apply(t(apply(order.mtx,1,cumsum)),1,max))

## 5 permute gene labels  
ES.B=matrix(NA,nrow(ES.0),nperm)
for(t1 in 1:nperm){
	if (nrow(order.mtx)>1){
	order.mtx.perm=order.mtx[,sample(ncol(order.mtx))]
	} else {
	order.mtx.perm=t(as.matrix(order.mtx[,sample(ncol(order.mtx))]))
	}
order.cumsum=t(apply(order.mtx.perm,1,cumsum))
ES.B[,t1]=apply(order.cumsum,1,max)
}
rownames(ES.B)=rownames(order.mtx)

## 6 pvalues calculation
N.X=apply(set2allgenes.mtx,1,sum)
N.Y=ncol(set2allgenes.mtx)-N.X
N=N.X* N.Y/(N.X + N.Y)


## 7 qvalues and qvalues calculation
enrich.out=pqvalues.compute(ES.0,ES.B,Stat.type='Tstat')


#### return(list(Tstat.set.0=ES.pval.0,Tstat.set.B=ES.pval.B,pvalue.set.0=enrich.out$pvalue.0,pvalue.set.B=enrich.out$pvalue.B,qvalue.set.0=enrich.out$qvalue.0))
return(list(pvalue.set.0=enrich.out$pvalue.0,pvalue.set.B=enrich.out$pvalue.B,qvalue.set.0=enrich.out$qvalue.0))

}
