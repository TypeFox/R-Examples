#' @title Meta analyze burden tests from multiple studies
#'   
#' @description Takes as input `seqMeta` objects (from the
#'   \code{\link{prepScores}} function), and meta-analyzes the corresponding
#'   burden test.
#'
#' @inheritParams singlesnpMeta   
#' @param wts weights for the burden test, as a function of maf, or a character string specifying weights in the SNP Info file.
#' @param mafRange Range of MAF's to include in the analysis (endpoints included).  Default is all SNPs (0 <= MAF <= 0.5).
#' 
#' @details This function uses the scores and their variances available in a
#'   seqMeta object to perform burden tests. Though coefficients are reported,
#'   the tests are formally score tests, and the coefficients can be thought of
#'   as one-step approximations to those reported in a Wald test.
#'   
#' @return a data frame with the following columns:
#'   \item{gene}{the name of the gene or unit of aggregation being meta analyzed}
#'   \item{p}{the p-value from the burden tests}
#'   \item{beta}{approximate coefficient for the effect of genotype}
#'   \item{se}{approximate standard error for the effect of genotype}
#'   \item{cmafTotal}{the cumulative minor allele frequency of the gene}
#'   \item{cmafUsed}{the cumulative minor allele frequency of snps used in the
#'   analysis}
#'   \item{nsnpsTotal}{the number of snps in the gene}
#'   \item{nsnpsUsed}{the number of snps used in the analysis}
#'   \item{nmiss}{The number of `missing` SNPs. For a gene with a single SNP
#'   this is the number of individuals which do not contribute to the analysis,
#'   due to studies that did not report results for that SNP. For a gene with
#'   multiple SNPs, is totalled over the gene.}
#'
#' @author Arie Voorman, Jennifer Brody
#' @seealso 
#' \code{\link{skatMeta}}
#' \code{\link{skatOMeta}}
#' \code{\link{prepScores}} 
#' 
#' @examples 
#' ###load example data for two studies:
#' ### see ?seqMetaExample	
#' data(seqMetaExample)
#' 
#' ####run on each cohort:
#' cohort1 <- prepScores(Z=Z1, y~1, SNPInfo=SNPInfo, data=pheno1)
#' cohort2 <- prepScores(Z=Z2, y~1, SNPInfo=SNPInfo, data=pheno2)
#' 
#' #### combine results:
#' out <- burdenMeta(cohort1, cohort2, SNPInfo = SNPInfo, mafRange=c(0,.01))
#' head(out)
#' 
#' \dontrun{
#' ##### Compare with analysis on full data set:
#' bigZ <- matrix(NA,2*n,nrow(SNPInfo))
#' colnames(bigZ) <- SNPInfo$Name
#' 
#' for(gene in unique(SNPInfo$gene)) {
#'    snp.names <- SNPInfo$Name[SNPInfo$gene == gene]
#'      bigZ[1:n,SNPInfo$gene == gene][ , snp.names \%in\% colnames(Z1)] <- 
#'                    Z1[ , na.omit(match(snp.names,colnames(Z1)))]
#'      bigZ[(n+1):(2*n),SNPInfo$gene == gene][ , snp.names \%in\% colnames(Z2)] <- 
#'                    Z2[ , na.omit(match(snp.names,colnames(Z2)))]
#' }
#' 
#' pheno <- rbind(pheno1[, c("y","sex","bmi")],pheno2[,c("y","sex","bmi")])
#' burden.p <- c(by(SNPInfo$Name, SNPInfo$gene, function(snp.names) {
#'    inds <- match(snp.names,colnames(bigZ)) burden <- rowSums(bigZ[,na.omit(inds)],na.rm=TRUE)
#'    mod <- lm(y~burden + gl(2,nrow(pheno1)),data=pheno)
#'    summary(mod)$coef[2,4]
#' }))
#' 
#' head(cbind(out$p,burden.p))
#' 
#' #will be slightly different:
#' plot(y=out$p,x=burden.p, ylab = "burden meta p-values", xlab = "complete data p-values")
#' }
#' 
#' @export
burdenMeta <- function(..., SNPInfo=NULL, wts=1, snpNames="Name", aggregateBy="gene", mafRange=c(0,0.5), verbose=FALSE) {
	cl <- match.call(expand.dots = FALSE)
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
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
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","beta","se","cmafTotal", "cmafUsed","nsnpsTotal","nsnpsUsed","nmiss")))
	colnames(res.numeric) <- c("p","beta","se","cmafTotal", "cmafUsed","nsnpsTotal","nsnpsUsed","nmiss")

	if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- utils::txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }

	snp.names.list <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	
	ri <- 0
	for(gene in genelist){
		
		ri <- ri+1
		nsnps.sub <- length(snp.names.list[[gene]])
		
		maf <- vcount <- numeric(nsnps.sub)
		n.total <- numeric(nsnps.sub)
		n.miss <- numeric(nsnps.sub)
		big.cov <- matrix(0, nsnps.sub,nsnps.sub)
		mscores <- numeric(nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- get(cohortNames[[cohort.k]],envir=ev)[[gene]]
			
			if(!is.null(cohort.gene)){
				#cohort.gene <- lapply(cohort.gene,function(x){replace(x,is.nan(x),0)})
				sub <- match(snp.names.list[[gene]],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
						cohort.gene$cov <- as.matrix(cohort.gene$cov)[sub,sub,drop=FALSE]
						cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0
						
						cohort.gene$maf <- cohort.gene$maf[sub]
						cohort.gene$maf[is.na(sub)] <- -1
						
						cohort.gene$scores <- cohort.gene$scores[sub]
						cohort.gene$scores[is.na(sub)] <- 0
				}
				#imputation
				n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
				n.total[cohort.gene$maf >= 0] <- n.total[cohort.gene$maf >= 0]+cohort.gene$n
				cohort.gene$maf[cohort.gene$maf < 0] <- 0
				
				maf <- maf + 2*cohort.gene$maf*cohort.gene$n
				mscores <- mscores + cohort.gene$scores/cohort.gene$sey^2
				big.cov <- big.cov + cohort.gene$cov/cohort.gene$sey^2
				vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
			} else {
				n.miss <- n.miss + get(cohortNames[[cohort.k]],envir=parent.frame())[[1]]$n
			}
		}

		maf <- maf/(2*n.total)
		maf[is.nan(maf)] <- 0
		
		flip <- maf > 0.5
		mscores[flip] <- -mscores[flip]
		big.cov[flip,!flip] <- -big.cov[flip,!flip]
		big.cov[!flip,flip] <- -big.cov[!flip,flip]
		maf <- pmin(maf,1-maf)
		
		if(is.function(wts)){
			tmpwts <- ifelse(maf > 0, wts(maf),0)
		} else if(is.character(wts)){
			tmpwts <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,wts])
		} else {
			tmpwts <- rep(1,length(maf))
		}
		tmpwts <- as.numeric(tmpwts)
		
		if( !all(mafRange == c(0,0.5))){
			keep <- (maf > min(mafRange)) & (maf <= max(mafRange))
			tmpwts[!keep] <- 0
		}
		
		vary.ave <- vary.ave/max(n.total,na.rm=TRUE)
		bscorevar <- as.numeric(tmpwts%*%big.cov%*%tmpwts)
		bscore <- as.numeric(tmpwts%*%mscores)
		
		res.numeric[ri,"beta"] = ifelse(bscorevar != 0, bscore/bscorevar, NA)
		res.numeric[ri,"se"] = ifelse(bscorevar != 0, sqrt(1/bscorevar), NA)
		res.numeric[ri,"cmafTotal"] = sum(maf,na.rm=TRUE)
		res.numeric[ri,"cmafUsed"] = sum(maf[tmpwts != 0],na.rm = TRUE)
		res.numeric[ri,"nsnpsTotal"] = length(maf)
		res.numeric[ri,"nmiss"] = sum(n.miss[tmpwts != 0], na.rm = T)
		res.numeric[ri,"nsnpsUsed"] = sum(tmpwts != 0)
		res.numeric[ri,"p"] = ifelse(bscorevar !=0,stats::pchisq(bscore^2/bscorevar,lower.tail=FALSE,df=1),NA)
		if(verbose){
			pb.i <- pb.i+1
			utils::setTxtProgressBar(pb, pb.i)
		}
	
	}
	if(verbose) close(pb)
	return(cbind(res.strings,res.numeric))
}
