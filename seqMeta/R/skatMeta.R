#' @title Combine SKAT analyses from one or more studies
#'   
#' @description Takes as input `seqMeta` objects (from the
#'   \code{\link{prepScores}} function), and meta-analyzes them.
#'
#' @inheritParams singlesnpMeta
#' @inheritParams burdenMeta
#' @param wts Either a function to calculate testing weights, or a character
#'   specifying a vector of weights in the SNPInfo file. For skatMeta the
#'   default are the `beta' weights.
#' @param method p-value calculation method. Default is 'saddlepoint',
#'   'integration' is the Davies method used in the SKAT package. See
#'   pchisqsum() for more details.
#' 
#' @details \code{skatMeta} implements an efficient SKAT meta analysis by
#'   meta-analyzing scores statistics and their variances. 
#'   
#'   Note: all studies must use coordinated SNP Info files - that is, the SNP
#'   names and gene definitions must be the same.
#'   
#'   Please see the package vignette for more details.
#'   
#' @return a data frame with the following columns:
#'   \item{gene}{the name of the gene or unit of aggregation being meta analyzed}
#'   \item{p}{p-value of the SKAT test.}
#'   \item{Q}{The SKAT Q-statistic, defined as sum_j w_jS_j, where S_j is the
#'   squared score for SNP j, and w_j is a weight.}
#'   \item{cmaf}{The cumulative minor allele frequency.}
#'   \item{nmiss}{The number of `missing` SNPs. For a gene with a single SNP
#'   this is the number of individuals which do not contribute to the analysis,
#'   due to studies that did not report results for that SNP. For a gene with
#'   multiple SNPs, is totalled over the gene.}
#'   \item{nsnps}{The number of SNPs in the gene.}
#'   
#' @references Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X.
#'   (2011) Rare Variant Association Testing for Sequencing Data Using the
#'   Sequence Kernel Association Test (SKAT). American Journal of Human
#'   Genetics.
#'   
#' @author Arie Voorman, Jennifer Brody
#' @seealso 
#' \code{\link{prepScores}}
#' \code{\link{burdenMeta}} 
#' \code{\link{singlesnpMeta}}
#' \code{\link{skatOMeta}}
#' 
#' @examples 
#' ###load example data for two studies:
#' ### see ?seqMetaExample	
#' data(seqMetaExample)
#' 
#' ####run on each study:
#' cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
#' cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo=SNPInfo, kins=kins, data=pheno2)
#' 
#' #### combine results:
#' ##skat
#' out <- skatMeta(cohort1, cohort2, SNPInfo = SNPInfo)
#' head(out)
#' 
#' \dontrun{
#' ##T1 test
#' out.t1 <- burdenMeta(cohort1,cohort2, SNPInfo = SNPInfo, mafRange = c(0,0.01))
#' head(out.t1)
#' 
#' ##single snp tests:
#' out.ss <- singlesnpMeta(cohort1,cohort2, SNPInfo = SNPInfo)
#' head(out.ss)
#' 
#' ########################
#' ####binary data
#' 
#' cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
#' out.bin <- skatMeta(cohort1, SNPInfo=SNPInfo)
#' head(out.bin)
#' 
#' ####################
#' ####survival data
#' cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo=SNPInfo, data=pheno1)
#' out.surv <- skatMeta(cohort1, SNPInfo=SNPInfo)
#' head(out.surv)
#' 
#' ##### Compare with SKAT on full data set
#' require(SKAT)
#' n <- nrow(pheno1)
#' bigZ <- matrix(NA,2*n,nrow(SNPInfo))
#' colnames(bigZ) <- SNPInfo$Name
#' 
#' for(gene in unique(SNPInfo$gene)) {
#'  snp.names <- SNPInfo$Name[SNPInfo$gene == gene]
#'    bigZ[1:n,SNPInfo$gene == gene][ , snp.names \%in\% colnames(Z1)] <- 
#'                    Z1[ , na.omit(match(snp.names,colnames(Z1)))]
#'    bigZ[(n+1):(2*n),SNPInfo$gene == gene][ , snp.names \%in\% colnames(Z2)] <- 
#'                    Z2[ , na.omit(match(snp.names,colnames(Z2)))]
#' }
#' 
#' pheno <- rbind(pheno1[,c("y","sex","bmi")], pheno2[,c("y","sex","bmi")])
#' 
#' obj <- SKAT_Null_Model(y~sex+bmi+gl(2,nrow(pheno1)), data=pheno)
#' skat.pkg.p <- c(by(SNPInfo$Name, SNPInfo$gene, function(snp.names) {
#'            inds <- match(snp.names,colnames(bigZ))
#'            if(sum(!is.na(inds)) ==0 ) return(1)
#'            SKAT(bigZ[,na.omit(inds)],obj, is_check=TRUE, missing=1)$p.value
#'            }))
#' 
#' head(cbind(out$p,skat.pkg.p))
#' 
#' #Note: SKAT ignores family strucutre, resulting in p-values that are systematically too small: 
#' plot(y=out$p,x=skat.pkg.p, ylab = "SKAT meta p-values", xlab = "SKAT p-values")
#' abline(0,1)
#' 
#' ignore family structure:
#' cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
#' cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo=SNPInfo, data=pheno2)
#'
#' out.nofam <- skatMeta(cohort1,cohort2,SNPInfo=SNPInfo)
#' plot(y=out.nofam$p,x=skat.pkg.p, ylab = "SKAT meta p-values", xlab = "SKAT p-values")
#' abline(0,1)
#' }
#' 
#' @export
skatMeta <- function(..., SNPInfo=NULL, wts=function(maf){ stats::dbeta(maf,1,25) }, method="saddlepoint", snpNames="Name", aggregateBy="gene", mafRange=c(0,0.5), verbose=FALSE) {
	cl <- match.call(expand.dots = FALSE)
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("seqMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	} else {
	  SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy, wt1=wts)
	}
	
	genelist <- stats::na.omit(unique(SNPInfo[,aggregateBy]))
	cohortNames <- lapply(cl[[2]],as.character)
	ncohort <- length(cohortNames)
	
	ev <- parent.frame()
	classes <- unlist(lapply(cohortNames,function(name){class(get(name,envir=ev))})) 
	if(!all(classes == "seqMeta" | classes == "skatCohort") ){
	 	stop("an argument to ... is not a seqMeta object!")
	}
	
	res.strings <- data.frame("gene"=genelist,stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","Qmeta","cmaf","nmiss", "nsnps")))
	colnames(res.numeric) <- c("p","Qmeta","cmaf","nmiss", "nsnps")		
	
    if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- utils::txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }
    ri <- 0
    snp.names.list <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	for(gene in genelist){		
		ri <- ri+1
		nsnps.sub <- length(snp.names.list[[gene]])
		
		n.miss <- n.total <- mscores <- maf <- numeric(nsnps.sub)
		big.cov <- matrix(0, nsnps.sub,nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- get(cohortNames[[cohort.k]],envir=ev)[[gene]]
			
			if(!is.null(cohort.gene)){
				#cohort.gene <- lapply(cohort.gene,function(x){replace(x,is.nan(x),0)})
				sub <- match(snp.names.list[[gene]],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
							#if(any(is.na(sub))) warning("Some SNPs were not in SNPInfo file for gene ", gene," and cohort ",names(cohorts)[cohort.k])
							cohort.gene$cov <- as.matrix(cohort.gene$cov)[sub,sub,drop=FALSE]
							cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0
							
							cohort.gene$maf <- cohort.gene$maf[sub]
							cohort.gene$maf[is.na(sub)] <- -1
							
							cohort.gene$scores <- cohort.gene$scores[sub]
							cohort.gene$scores[is.na(sub)] <- 0
					}				
					
					n.total <- n.total + (cohort.gene$maf >= 0)*cohort.gene$n
					n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
					cohort.gene$maf[cohort.gene$maf < 0] <- 0
					
					mscores <- mscores + cohort.gene$scores/cohort.gene$sey^2
					maf <- maf + 2*cohort.gene$maf*(cohort.gene$n)
					big.cov <- big.cov + cohort.gene$cov/cohort.gene$sey^2
					vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
			}else{
				n.miss <- n.miss + get(cohortNames[[cohort.k]],envir=parent.frame())[[1]]$n
			} 
		}
		if(any(maf >0)){ 
			maf <- maf/(2*n.total)
			maf[is.nan(maf)] <- 0

			flip <- maf > 0.5
			mscores[flip] <- -mscores[flip]
			big.cov[flip,!flip] <- -big.cov[flip,!flip]
			big.cov[!flip,flip] <- -big.cov[!flip,flip]
			maf <- pmin(maf,1-maf)
		}
		if(is.function(wts)){
			tmpwts <- ifelse(maf > 0, wts(maf),0)
		} else if(is.character(wts)){
			tmpwts <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,wts])
		} else {
			tmpwts <- rep(1,length(maf))
		}
		tmpwts <- as.numeric(tmpwts)
		
		if( !all(mafRange == c(0,0.5))){
		  keep <- (maf >= min(mafRange)) & (maf <= max(mafRange))
      	  tmpwts[!keep] <- 0
		}
    
		if(length(maf) > 0){
		  Qmeta <- sum((tmpwts*mscores)^2, na.rm=TRUE)
		  wcov <- tmpwts*t(t(big.cov)*tmpwts)
			lambda<-eigen(zapsmall(wcov),symmetric=TRUE)$values
    } else {
    	Qmeta <- 0
      lambda <- 0
    }
    	
    if(any(lambda > 0)){
    		p<-pchisqsum2(Qmeta,lambda,method=method)$p
		} else {
			p<-1
		}
		res.numeric[ri,"p"] = p
		res.numeric[ri,"Qmeta"] = Qmeta
		res.numeric[ri,"cmaf"] = sum(maf[tmpwts > 0],na.rm=TRUE)
		res.numeric[ri,"nsnps"] = sum(maf[tmpwts > 0] != 0, na.rm =T)
		res.numeric[ri,"nmiss"] = sum(n.miss, na.rm =T)
		if(verbose){
			pb.i <- pb.i+1
			utils::setTxtProgressBar(pb, pb.i)
		}
	}	
	if(verbose) close(pb)
	return(cbind(res.strings,res.numeric))
}
