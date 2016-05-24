#' @title Meta analyze single snp effects from multiple studies
#'   
#' @description Takes as input `seqMeta` objects (from the
#'   \code{\link{prepScores}} function), and meta analyzes them.
#'   
#' @param ... \code{seqMeta} objects
#' @param SNPInfo The SNP Info file.  This should contain the fields listed in
#'   snpNames and aggregateBy. Only SNPs in this table will be meta analyzed, so
#'   this may be used to restrict the analysis.
#' @param snpNames The field of SNPInfo where the SNP identifiers are found. 
#'   Default is 'Name'
#' @param aggregateBy     The field of SNPInfo on which the skat results were
#'   aggregated.  Default is 'gene'. Though gene groupings are not explicitely
#'   required for single snp analysis, it is required to find where single snp
#'   information is stored in the seqMeta objects.
#' @param studyBetas Whether or not to include study-level effects in the
#'   output.
#' @param verbose logical. Whether progress bars should be printed.
#'   
#' @details This function meta analyzes score tests for single variant effects.
#'   Though the test is formally a score test, coefficients and standard errors
#'   are also returned, which can be interpreted as a `one-step` approximation
#'   to the maximum likelihood estimates.
#'   
#' @return a data frame with the gene, snp name, meta analysis.
#' @references Lin, DY and Zeng, D. On the relative efficiency of using summary statistics versus individual-level data in meta-analysis. Biometrika. 2010.
#' 
#' @author Arie Voorman, Jennifer Brody
#' @seealso 
#' \code{\link{prepScores}} 
#' \code{\link{burdenMeta}}
#' \code{\link{skatMeta}}
#' \code{\link{skatOMeta}}
#' 
#' @examples 
#' ###load example data for two studies:
#' ### see ?seqMetaExample
#' data(seqMetaExample)
#' 
#' ####run on each study:
#' cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
#' cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, data =pheno2)
#' 
#' #### combine results:
#' out <- singlesnpMeta(cohort1, cohort2, SNPInfo = SNPInfo)
#' head(out)
#' 
#' \dontrun{
#' ##compare
#' bigZ <- matrix(NA,2*n,nrow(SNPInfo))
#' colnames(bigZ) <- SNPInfo$Name
#' for(gene in unique(SNPInfo$gene)) {
#'    snp.names <- SNPInfo$Name[SNPInfo$gene == gene]
#'      bigZ[1:n,SNPInfo$gene == gene][, snp.names \%in\% colnames(Z1)] <- 
#'              Z1[, na.omit(match(snp.names,colnames(Z1)))]
#'      bigZ[(n+1):(2*n),SNPInfo$gene == gene][, snp.names \%in\% colnames(Z2)] <- 
#'              Z2[, na.omit(match(snp.names,colnames(Z2)))]
#' }
#' 
#' pheno <- rbind(pheno1[ ,c("y","sex","bmi")], pheno2[ , c("y","sex","bmi")])
#' out3 <- apply(bigZ, 2, function(z) {
#'          if(any(!is.na(z))) {
#'            z[is.na(z)] <- mean(z,na.rm=TRUE)
#'            mod <- lm(y ~ sex+bmi+c(rep(1,n),rep(0,n))+z, data=pheno)
#'            beta <- mod$coef["z"]
#'            se <- sqrt(vcov(mod)["z", "z"])
#'            p <- pchisq( (beta/se)^2,df=1,lower=F)
#'            return(c(beta,se,p))
#'          } else {
#'            return(c(0,Inf,1))
#'          }
#'  }) 
#'  out3 <- t(out3[,match(out$Name,colnames(out3))])
#'  
#'  ##plot
#'  par(mfrow=c(2,2))
#'  plot(x=out3[,1],y=out$beta, xlab="complete data (Wald)", 
#'       ylab="meta-analysis (Score)", main="coefficients"); abline(0,1)
#'  plot(x=out3[,2],y=out$se, xlab="complete data (Wald)", 
#'       ylab="meta-analysis (Score)", main="standard errors"); abline(0,1)
#'  plot(x=out3[,3],y=out$p, xlab="complete data (Wald)", 
#'       ylab="meta-analysis (Score)", main="p-values"); abline(0,1)
#'  }
#' @export
singlesnpMeta <- function(..., SNPInfo=NULL, snpNames="Name", aggregateBy="gene", studyBetas=TRUE, verbose=FALSE) {
	cl <- match.call(expand.dots=FALSE)
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("seqMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	} else {
	  SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy)
	}
	
	genelist <- stats::na.omit(unique(SNPInfo[,aggregateBy]))
	cohortNames <- lapply(cl[[2]],as.character)
	ncohort <- length(cohortNames)
	
  ev <- parent.frame()
	classes <- unlist(lapply(cohortNames,function(name){class(get(name,envir=ev))})) 
	if(!all(classes == "seqMeta" | classes == "skatCohort") ){
	 	stop("an argument to ... is not a seqMeta object!")
	}

	cohort.names <- unlist(cohortNames)
	#for(i in 1:ncohort) cohort.names <- c(cohort.names,as.character(cl[[2]][[i]]))
	
	res.strings <- data.frame("gene"=as.character(SNPInfo[,aggregateBy]),
					"Name" = as.character(SNPInfo[,snpNames]),stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","maf","caf","nmiss","ntotal", "beta" ,"se" ))) # v 1.6 added caf
	colnames(res.numeric) <- c("p","maf","caf","nmiss","ntotal", "beta" ,"se" ) # v 1.6 added caf
					
	if(studyBetas){
		resdf.cohort = matrix(NA,nrow=nrow(SNPInfo),ncol=2*ncohort)
		oo <- as.numeric(gl(ncohort,2))
		oo[1:ncohort*2]<- oo[1:ncohort*2]+ncohort
		colnames(resdf.cohort) <- c(paste(c("beta"),cohort.names,sep="."),paste(c("se"),cohort.names ,sep="."))[oo]
	}
	
	if(verbose){
    	cat("\n Meta Analyzing... Progress:\n")
    	pb <- utils::txtProgressBar(min = 0, max = length(genelist), style = 3)
    	pb.i <- 0
    }

	#gene.namelist <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	ris <- split(1:nrow(SNPInfo),SNPInfo[,aggregateBy])
	snp.names.list <- split(SNPInfo[,snpNames],SNPInfo[,aggregateBy])
	for(gene in genelist){
		ri <- ris[[gene]]
		nsnps.sub <- length(snp.names.list[[gene]])
		
		n.total <- n.miss <- maf <- caf <- vcount <- scorevar <- mscore <- numeric(nsnps.sub)
		big.cov <- matrix(0, nsnps.sub,nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- get(cohortNames[[cohort.k]],envir=ev)[[gene]]
			
			if(!is.null(cohort.gene)){
				#cohort.gene <- lapply(cohort.gene,function(x){replace(x,is.nan(x),0)})
				sub <- match(snp.names.list[[gene]],colnames(cohort.gene$cov))
				if(any(is.na(sub)) | any(sub != 1:length(sub), na.rm=TRUE) | length(cohort.gene$maf) > nsnps.sub){
						#if(any(is.na(sub))) warning("Some SNPs were not in SNPInfo file for gene ", gene," and cohort ",cohort.names[[cohort.k]])		
						cohort.gene$cov <-  as.matrix(cohort.gene$cov)[sub,sub,drop=FALSE]
						cohort.gene$cov[is.na(sub),] <- cohort.gene$cov[,is.na(sub)] <- 0
						
						cohort.gene$maf <- cohort.gene$maf[sub]
						cohort.gene$maf[is.na(sub)] <- -1
						
						cohort.gene$scores <- cohort.gene$scores[sub]
						cohort.gene$scores[is.na(sub)] <- 0
					}				
					
					n.total <- n.total + (cohort.gene$maf >= 0)*cohort.gene$n				
					n.miss[cohort.gene$maf < 0] <- n.miss[cohort.gene$maf < 0] + cohort.gene$n
					cohort.gene$maf[cohort.gene$maf < 0] <- 0
	
					mscore <- mscore + cohort.gene$score/cohort.gene$sey^2	
					scorevar <- scorevar + diag(cohort.gene$cov)/cohort.gene$sey^2
					maf <- maf + 2*cohort.gene$maf*(cohort.gene$n)
					vary.ave <- vary.ave + max(cohort.gene$n,na.rm=T)*cohort.gene$sey^2
		
					if(studyBetas){
						resdf.cohort[ri,paste(c("beta"),cohort.names[cohort.k] ,sep=".")] <- ifelse(diag(cohort.gene$cov) >0,
							cohort.gene$score/diag(cohort.gene$cov),
							NA)
						resdf.cohort[ri,paste(c("se"), cohort.names[cohort.k] ,sep=".")] <- ifelse(diag(cohort.gene$cov) >0,
							cohort.gene$sey/sqrt(diag(cohort.gene$cov)),
							Inf)
					}
			} else {
				n.miss <- n.miss + get(cohortNames[[cohort.k]],envir=parent.frame())[[1]]$n
			} 
		}
		
		vary.ave <- vary.ave/n.total

		maf <- maf/(2*n.total)
		maf[is.nan(maf)] <- 0
		caf = maf ## JB v 1.5
		maf <- sapply(maf, function(x){min(x,1-x)})
				
		res.numeric[ri,"beta"] <- ifelse(scorevar != 0, mscore/scorevar, NA)
		res.numeric[ri,"se"] <- ifelse(scorevar != 0, sqrt(1/scorevar), NA)
		res.numeric[ri,"maf"] <- maf
		res.numeric[ri,"caf"] <- caf #JB v 1.5
		res.numeric[ri,"nmiss"] <- n.miss
		res.numeric[ri,"ntotal"] <- n.total
		res.numeric[ri,"p"] <- ifelse(scorevar != 0, stats::pchisq(mscore^2/scorevar, lower.tail = FALSE,df = 1), NA)

		if(verbose){
			pb.i <- pb.i+1
			utils::setTxtProgressBar(pb, pb.i)
		}
	}
	if(verbose) close(pb)
	
	if(studyBetas){
		return(cbind(res.strings,res.numeric,resdf.cohort))
	} else {
		return(cbind(res.strings,res.numeric))
	}
}
