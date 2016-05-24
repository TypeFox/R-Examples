#' @title Run SKAT on data from a single cohort, conditional on specified SNP 
#'   effects
#'   
#' @description This function works exactly as \code{\link{prepScores}}, but 
#'   with the additional argument `adjustments' specifying genes for which 
#'   conditional analyses are desired, and which SNPs to condition on.
#'   
#' @param adjustments A data frame of the same format at SNPInfo, pairing genes 
#'   to analyze with snp
#' @inheritParams prepScores
#'   
#' @details This function has the same syntax as \code{\link{prepCondScores}},
#'   but requires an extra argument `adjustments`. This is a data frame of the
#'   same format as the SNPInfo, i.e. with a `snpNames` and `aggregateBy`
#'   columns. The function works by looping through the genes in the adjustment
#'   file, adding the corresponding SNPs to the null model.  For instance, if
#'   one wants to adjuste `gene1` for SNPs a and b (which need not be in gene
#'   1), and `gene2' for SNPs c, the adjustments would be something like
#'   \code{adjustments = data.frame(Name = c("a","b","c"), gene =
#'   c("gene1","gene1","gene2"))} 
#'   
#'   See the examples for an illustration.
#' 
#' @return 	an object of class 'seqMeta'. Note that unlike output from the
#'   function \code{\link{prepScores}}, the null models in each element of the
#'   list may be different. When meta analyzing these, it may be good to subset
#'   the SNPInfo file to the genes of interest.
#'   
#' @author Arie Voorman, Jennifer Brody
#' 
#' @seealso \code{\link{prepScores}}
#' \code{\link{skatMeta}}
#' \code{\link{burdenMeta}}
#' \code{\link{singlesnpMeta}}
#' 
#' @examples
#' ###load example data for two studies:
#' ### see ?seqMetaExample	
#' data(seqMetaExample)
#' 
#' #specify adjustment variables
#' adjustments <- SNPInfo[c(1:3, 20,100), ]
#' adjustments
#' 
#' ####run on each study:
#' cohort1.adj <- prepCondScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, 
#'                 adjustments=adjustments, data =pheno1)
#' cohort2.adj <- prepCondScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
#'                 adjustments=adjustments, kins=kins, data=pheno2)
#' 
#' SNPInfo.sub <- subset(SNPInfo, (SNPInfo$gene \%in\% adjustments$gene) & 
#'                                 !(SNPInfo$Name \%in\% adjustments$Name) )
#' 
#' #skat
#' out.skat <- skatMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub)
#' head(out.skat)
#' 
#' ##T1 test
#' out.t1 <- burdenMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub, mafRange = c(0,0.01))
#' head(out.t1)
#' 
#' ##single snp tests:
#' out.ss <- singlesnpMeta(cohort1.adj,cohort2.adj, SNPInfo = SNPInfo.sub)
#' head(out.ss)
#' 
#' @export
prepCondScores <- function(Z, formula, family = stats::gaussian(), SNPInfo=NULL, adjustments= NULL, snpNames = "Name", aggregateBy = "gene",kins=NULL, sparse = TRUE, data=parent.frame()){
	if(is.null(SNPInfo)) stop("SNPInfo file must be provided!")
	if(is.null(adjustments)) stop("adjustments must be provided!")
	
	re <- NULL
	#update formula
	null.formula.adj <- update(formula, ".~.+ geno.adj")
	
	genes <- unique(adjustments[,aggregateBy])
	for(gene in genes){
		#subset to the gene:
		minisnpinfo <- SNPInfo[SNPInfo[,aggregateBy] == gene,]
		
		#remove adjustment snp
		adjustmentsnps <- subset(adjustments[,snpNames],adjustments[,aggregateBy] ==gene)
		minisnpinfo <- subset(minisnpinfo, !(minisnpinfo[,snpNames] %in% adjustmentsnps) , drop=FALSE)
		
		#if there are more snps in this gene....
		if(nrow(minisnpinfo) > 0){
			#add genotype to adjustment variables
			zz <- subset(Z,select = colnames(Z) %in% adjustmentsnps,drop=FALSE)
			zz <- apply(zz,2,function(z){
				z[is.na(z)] <- mean(z,na.rm=T)
				z
			})
			data$geno.adj <- zz
			
			Z.sub <- subset(Z, select = colnames(Z) %in% minisnpinfo[,snpNames],drop=FALSE)
			
			#if our sample has more snps in this gene...
			if(ncol(Z.sub) > 0){	
				re.tmp <- prepScores(Z= Z.sub, formula = null.formula.adj, family = family, 
					SNPInfo = minisnpinfo, aggregateBy = aggregateBy, snpNames = snpNames, kins=kins, sparse=sparse, data=data)
				class(re.tmp) <- "list"
				if(is.null(re)){
					re <- re.tmp
				} else {
					re <- c(re,re.tmp)
				}		
			} 
		}	
	}
	class(re) <- "seqMeta"
	return(re)
}