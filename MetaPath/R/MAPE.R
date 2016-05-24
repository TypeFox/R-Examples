MAPE <-
function(arraydata,pathway.DB,resp.type=c('twoclass','multiclass','continuous', 'survival'),stat=c('maxP','minP','rth','Fisher'),rth.value=NULL,
	permutation=c('sample','gene'),nperm=500,size.min=15,size.max=500,knn.neighbors=10,qvalue.cal=c('permute','estimate')){


## load required packages
##packages.meta=c('Biobase','genefilter','GSEABase','limma')
##for(package in packages.meta) do.call(require, list(package))


## 1. check the input auguments (eset or list)
###  list
# x--exprs data
# y-- label
# z-- censoring.status
# geneid
# samplename

### eset
### eset$label
### eset$censoring.status


## check the input arguments

stat=match.arg(stat)
resp.type=match.arg(resp.type)
permutation=match.arg(permutation)
qvalue.cal=match.arg(qvalue.cal)

data.type=class(arraydata[[1]])
if( length(arraydata)>1) {
	for(t1 in 2:length(arraydata)){
		if (data.type!=class(arraydata[[t1]])) stop('Please check: the data stucture of each study should be same')
	}
}

## eset
if (is(arraydata[[1]],'ExpressionSet')) {
	for(t1 in 1:length(arraydata)){
		if(is.null(arraydata[[t1]]$label)) stop('Please check: the group label can not be NULL')
		if(resp.type=='survival' & is.null(arraydata[[t1]]$censoring.status)) stop('Please check: the censoring.status can not be NULL if the resp.type is survival')
	}

	study=arraydata
		
	if(is.null(names(study))){
	names(study)=paste('study',1:length(study),sep='')
	} else {
	names(study)=names(arraydata)
	}		
}

## data list -> eset
if (is(arraydata[[1]],'list')) {
	study=list()
	for(t1 in 1:length(arraydata)){
	
	    if(is.null(arraydata[[t1]]$x)) stop('Please check: the expression values can not be NULL')
		if(is.null(arraydata[[t1]]$y)) stop('Please check: the group lable can not be NULL')
		if(resp.type=='survival' & is.null(arraydata[[t1]]$z)) stop('the censoring.status can not be NULL if the resp.type is survival')
		if(is.null(arraydata[[t1]]$geneid)) stop('Please check: the gene IDs can not be NULL')
		if(is.null(arraydata[[t1]]$samplename)) stop('Please check: the sample names can not be NULL')
	    
	    
		exprs=as.matrix(arraydata[[t1]]$x)
		rownames(exprs)=arraydata[[t1]]$geneid
		colnames(exprs)=arraydata[[t1]]$samplename

			if(resp.type!='survival'){
				pheno=as.data.frame(matrix(arraydata[[t1]]$y,length(arraydata[[t1]]$y),1))
				rownames(pheno)=arraydata[[t1]]$samplename
				colnames(pheno)='label'
			} else {
			    yy=as.data.frame(matrix(arraydata[[t1]]$y,length(arraydata[[t1]]$y),1))
			    zz=as.data.frame(matrix(arraydata[[t1]]$z,length(arraydata[[t1]]$z),1))
				pheno=as.data.frame(cbind(yy,zz))
				rownames(pheno)=arraydata[[t1]]$samplename
				colnames(pheno)=c('label','censoring.status')
			}

		phenoData=new("AnnotatedDataFrame", data =pheno)
		study[[t1]]<- new("ExpressionSet", exprs = exprs ,phenoData=phenoData )
	}

	if(is.null(names(arraydata))){
	names(study)=paste('study',1:length(arraydata),sep='')
	} else {
	names(study)=names(arraydata)
	}

}

	label='label';censoring.status='censoring.status'
	rm(arraydata)

# impute missing data
for(study.no in 1:length(study)){
	if(sum(is.na(exprs(study[[study.no]])))>0){
##	require(impute)
	exprs(study[[study.no]])=impute::impute.knn(exprs(study[[study.no]]),k=knn.neighbors)$data
    }
}	

gene.in.array=featureNames(study[[1]])
if(length(study)>1) {
	for(t1 in 2:length(study)){
		gene.in.array=unique(gene.in.array,featureNames(study[[t1]]))
	}
}

gene.in.DB=unique(unlist(geneIds(pathway.DB)))
set.name=names(pathway.DB)

gene.common=intersect(gene.in.array,gene.in.DB)

DB.matrix=matrix(0,length(set.name),length(gene.common))
rownames(DB.matrix)=set.name
colnames(DB.matrix)=gene.common

for(t1 in 1:length(set.name)){
gene=intersect(geneIds(pathway.DB[[t1]]),gene.common)
DB.matrix[set.name[t1],gene]=1
}
colnames(DB.matrix)=toupper(colnames(DB.matrix))

keep.idx=(apply(DB.matrix,1,sum)>=size.min & apply(DB.matrix,1,sum)<=size.max)
DB.matrix=DB.matrix[keep.idx,]


	
if(length(study)==1) {
	######### enrichment analysis for single study	
	madata=study[[1]]
	testlabel=madata[[label]]
	if (resp.type=="survival") 	censoring=madata[[censoring.status]]

	if (permutation=='gene'){
		enrich=Enrichment_KS_gene(madata=madata,label=testlabel,censoring=censoring,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,resp.type=resp.type)
	} else if (permutation=='sample'){
		enrich=Enrichment_KS_sample(madata=madata,label=testlabel,censoring=censoring,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,resp.type=resp.type)
	} else {
	stop('Please check: Wrong permutation methods.')
	}
	
	pvalue=as.data.frame(enrich$pvalue.set.0)
	colnames(pvalue)='pvalue'
	
	qvalue=as.data.frame(enrich$qvalue.set.0)
	colnames(qvalue)='qvalue'

	return(list(qvalue=qvalue,pvalue=pvalue))
	
	
  } else {
	########## enrichment analysis for multipe studies
	if (permutation=='gene'){

		cat("Performing MAPE_P analysis...\n")
		MAP_SET.obj=MAPE_P_gene_KS(study=study,label=label,censoring.status=censoring.status,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,stat=stat,rth.value=rth.value,resp.type)
		cat("Performing MAPE_G analysis...\n")
		MAP_GENE.obj=MAPE_G_gene_KS(study=study,label=label,censoring.status=censoring.status,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,stat=stat,rth.value=rth.value,resp.type)
		cat("Performing MAPE_I analysis...\n")
		MAP_I.obj=MAPE_I_KS(MAP_GENE.obj=MAP_GENE.obj, MAP_SET.obj=MAP_SET.obj, study=study) 

	} else if (permutation=='sample'){

		cat("Performing MAPE_P analysis...\n")
		MAP_SET.obj=MAPE_P_sample_KS(study=study,label=label,censoring.status=censoring.status,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,stat=stat,rth.value=rth.value,resp.type)
		cat("Performing MAPE_G analysis...\n")
		MAP_GENE.obj=MAPE_G_sample_KS(study=study,label=label,censoring.status=censoring.status,DB.matrix=DB.matrix,size.min=size.min,size.max=size.max,nperm=nperm,stat=stat,rth.value=rth.value,resp.type)
		cat("Performing MAPE_I analysis...\n")
		MAP_I.obj=MAPE_I_KS(MAP_GENE.obj=MAP_GENE.obj, MAP_SET.obj=MAP_SET.obj, study=study) 	

	} else {
	stop('Please check: Wrong permutation methods.')
	}

	## store qvalues of all gene sets in a matrix
	study.no=length(study)
	study.name=names(study)
	common.set.name=rownames(MAP_I.obj$qvalue.meta)
	qvalue.all=matrix(NA,length(common.set.name),(study.no+3))
	rownames(qvalue.all)=common.set.name
	colnames(qvalue.all)=c(names(study),'MAPE_P','MAPE_G','MAPE_I')
	qvalue.all[,1:study.no]=MAP_SET.obj$qvalue.set.study[common.set.name,]
	qvalue.all[,(study.no+1)]=MAP_SET.obj$qvalue.meta[common.set.name,]
	qvalue.all[,(study.no+2)]=MAP_GENE.obj$qvalue.meta[common.set.name,]
	qvalue.all[,(study.no+3)]=MAP_I.obj$qvalue.meta[common.set.name,]
	qvalue.all=as.data.frame(qvalue.all)

	############
	pvalue.all=matrix(NA,length(common.set.name),(study.no+3))
	rownames(pvalue.all)=common.set.name
	colnames(pvalue.all)=c(names(study),'MAPE_P','MAPE_G','MAPE_I')
	pvalue.all[,1:study.no]=MAP_SET.obj$pvalue.set.study[common.set.name,]
	pvalue.all[,(study.no+1)]=MAP_SET.obj$pvalue.meta[common.set.name,]
	pvalue.all[,(study.no+2)]=MAP_GENE.obj$pvalue.meta[common.set.name,]
	pvalue.all[,(study.no+3)]=MAP_I.obj$pvalue.meta[common.set.name,]
	pvalue.all=as.data.frame(pvalue.all)

if(qvalue.cal=='estimate'){
	for(t1 in 1:ncol(pvalue.all)){
		qvalue.all[,t1]=p.adjust(pvalue.all[,t1],"BH") 
	}
	
}

	return(list(qvalue=qvalue.all,pvalue=pvalue.all))
  }
 
 }
