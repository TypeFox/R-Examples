RNAgamma <-
function(formula, data, rnaprefix="ENSG", gammaShape=1,STT=NULL,
         pheno.type=c("case.control"),tagwise=F,perm=T,n.perm=1000,seed=12212012) 
{  
	set.seed(seed)
        rnaseqcounts <- data[,names(data)[substring(names(data),1,nchar(rnaprefix))==rnaprefix]]
	
	if(is.null(STT) & is.null(gammaShape)) stop("need either STT or gammaShape specified")
        if(!is.null(STT) & !is.null(gammaShape)) stop("specify only one of STT or gammaShape")
        if(!is.null(STT))  gammaShape<-STTtoShapeParameter(STT)
        phenotype <- as.character(formula)[2]
        pheno <- data[,phenotype]
        if(class(pheno) != "list") 
		pheno<-list(pheno)
        if(!perm) 
		n.perm<-0
		i<-2
        while(i<=n.perm+1) 
	{
		pheno[[i]]<-sample(pheno[[1]],replace=F)
		i<-i+1
        }
        rm(i)

  	dgeList <- DGEList(counts=t(rnaseqcounts),group=pheno[[1]])
        keep <- rowSums(cpm(dgeList)>100) >=2
        dgeList <- dgeList[keep,]
        dgeList$samples$lib.size <- colSums(dgeList$counts)
        dgeList <- calcNormFactors(dgeList)
        genenames <- row.names(dgeList$counts)
        
        rnafunc <- function(pheno.tmp) {
         dgeList$samples$group <- pheno.tmp
         if(tagwise==F) {
          if(as.character(formula)[3]!=".") {
             design <- model.matrix(as.formula(paste("~",as.character(formula)[3])),data=data)
             dgeList <- estimateGLMCommonDisp(dgeList, design)
          } else {
             dgeList <- estimateCommonDisp(dgeList)
          }
         } else {
          if(as.character(formula)[3]!=".") {
             design <- model.matrix(as.formula(paste("~",as.character(formula)[3])),data=data)
             dgeList <- estimateGLMTagwiseDisp(dgeList, design)
          } else {
             dgeList <- estimateTagwiseDisp(dgeList, trend="none")
          }
         }
	 et <- exactTest(dgeList)
	 pvalues <- et$table$PValue
       }
       pvalues <- sapply(pheno,rnafunc)
       row.names(pvalues) <- genenames
       
    gamma.stat<-rowSums(qgamma(1-signif(pvalues),shape=gammaShape),na.rm=T)
    if(perm) {
	gamma.perm.pvalue<-mean(gamma.stat[-1]>=gamma.stat[1])     
    } else {
	gamma.perm.pvalue<-NULL
    }
	gamma.obs.pvalue<-pgamma(gamma.stat[1],shape=gammaShape*dim(pvalues)[2],lower.tail=F)
    return(list(gamma.pvalue=gamma.obs.pvalue,perm.pvalue=gamma.perm.pvalue,gene.pvalues=pvalues))
}
