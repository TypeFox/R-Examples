#' @title Prepare scores for region based (meta) analysis
#'   
#' @description This function computes and organizes the neccesary output to 
#'   efficiently meta-analyze SKAT and other tests. Note that the tests are 
#'   *not* computed by these functions. The output must be passed to one of 
#'   \code{\link[seqMeta]{skatMeta}}, \code{\link[seqMeta]{burdenMeta}}, or 
#'   \code{\link[seqMeta]{singlesnpMeta}}.
#'   
#'   Unlike the SKAT package which operates on one gene at a time, these 
#'   functions are intended to operate on many genes, e.g. a whole exome, to 
#'   facilitate meta analysis of whole genomes or exomes.
#'   
#' @param Z A genotype matrix (dosage matrix) - rows correspond to individuals 
#'   and columns correspond to SNPs. Use 'NA' for missing values. The column 
#'   names of this matrix should correspond to SNP names in the SNP information 
#'   file.
#' @param formula Base formula, of the kind used in glm() - typically of the 
#'   form y~covariate1 + covariate2. For Cox models, the formula follows that of
#'   the coxph() function.
#' @param family either gaussian(), for continuous data, or binomial() for 0/1
#'   outcomes. Binary outcomes are not currently supported for family data.
#' @param SNPInfo SNP Info file - must contain fields given in 'snpName' and 
#'   'aggregateBy'.
#' @param snpNames The field of SNPInfo where the SNP identifiers are found. 
#'   Default is 'Name'.  See Details.
#' @param aggregateBy The field of SNPInfo on which the skat results were 
#'   aggregated. Default is 'gene'. For single snps which are intended only for 
#'   single variant analyses, it is recomended that they have a unique 
#'   identifier in this field.
#' @param data  data frame in which to find variables in the formula
#' @param kins  the kinship matrix for related individuals. Only supported for 
#'   family=gaussian(). See lmekin in the kinship2 package for more details.
#' @param sparse  whether or not to use a sparse Matrix approximation for dense 
#'   kinship matrices (defaults to TRUE).
#' @param verbose  logical. whether or not to print the progress bar.
#'   
#' @details This function computes the neccesary information to meta analyze 
#'   SKAT analyses: the individual SNP scores, their MAF, and a covariance 
#'   matrix for each unit of aggregation. Note that the SKAT test is *not* 
#'   calculated by this function. The output must be passed to one of 
#'   \code{\link[seqMeta]{skatMeta}}, \code{\link[seqMeta]{burdenMeta}}, or 
#'   \code{\link[seqMeta]{singlesnpMeta}}.
#'   
#'   A crucial component of SKAT and other region-based tests is a common unit 
#'   of aggregation accross studies. This is given in the SNP information file 
#'   (argument \code{SNPInfo}), which pairs SNPs to a unit of aggregation 
#'   (typically a gene). The additional arguments \code{snpNames} and 
#'   \code{aggregateBy} specify the columns of the SNP information file which 
#'   contain these pairings. Note that the column names of the genotype matrix 
#'   \code{Z} must match the names given in the \code{snpNames} field.
#'   
#'   Using \code{prepScores}, users are strongly recommended to use all SNPs, 
#'   even if they are monomorphic in your study. This is for two reasons; 
#'   firstly, monomorphic SNPs provide information about MAF across all studies;
#'   without providing the information we are unable to tell if a missing SNP 
#'   data was monomorphic in a study, or simply failed to genotype adequately in
#'   that study. Second, even if some SNPs will be filtered out of a particular 
#'   meta-analysis (e.g., because they are intronic or common) constructing 
#'   seqMeta objects describing all SNPs will reduce the workload for subsequent
#'   follow-up analyses.
#'   
#'   Note: to view results for a single study, one can pass a single seqMeta 
#'   object to a function for meta-analysis.
#'   
#' @return an object of class 'seqMeta'. This is a list, not meant for human 
#'   consumption, but to be fed to \code{skatMeta()} or another function. The 
#'   names of the list correspond to gene names. Each element in the list 
#'   contains
#'   
#'   \item{scores}{The scores (y-yhat)^t g}
#'   \item{cov}{The variance of the scores. When no covariates are used, this is the LD matrix.}
#'   \item{n}{The number of subjects}
#'   \item{maf}{The alternate allele frequency}
#'   \item{sey}{The residual standard error.}
#'   
#' @note For \code{prepCox}, the signed likelihood ratio statistic is used 
#'   instead of the score, as the score test is anti-conservative for 
#'   proportional hazards regression. The code for this routine is based on the 
#'   \code{coxph.fit} function from the \code{survival} package.
#'   
#'   Please see the package vignette for more details.
#'   
#' @author Arie Voorman, Jennifer Brody
#' 
#' @references Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011) Rare Variant Association Testing for Sequencing Data Using the Sequence Kernel Association Test (SKAT). American Journal of Human Genetics.
#' 
#' Chen H, Meigs JB, Dupuis J. Sequence Kernel Association Test for Quantitative Traits in Family Samples. Genetic Epidemiology. (To appear)
#' 
#' Lin, DY and Zeng, D. On the relative efficiency of using summary statistics versus individual-level data in meta-analysis. Biometrika. 2010.
#' 
#' @seealso \code{\link[seqMeta]{skatMeta}} \code{\link[seqMeta]{burdenMeta}} \code{\link[seqMeta]{singlesnpMeta}} \code{\link[seqMeta]{skatOMeta}} \code{\link[survival]{coxph}}
#' @examples
#' ###load example data for two studies:
#' ### see ?seqMetaExample
#' data(seqMetaExample)
#' 
#' ####run on each cohort:
#' cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
#' cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, kins=kins, data=pheno2)
#' 
#' #### combine results:
#' ##skat
#' out <- skatMeta(cohort1, cohort2, SNPInfo = SNPInfo)
#' head(out)
#' 
#' ##T1 test
#' out.t1 <- burdenMeta(cohort1,cohort2, SNPInfo = SNPInfo, mafRange = c(0,0.01))
#' head(out.t1)
#' 
#' ##single snp tests:
#' out.ss <- singlesnpMeta(cohort1,cohort2, SNPInfo = SNPInfo)
#' head(out.ss)
#' \dontrun{
#' ########################
#' ####binary data
#' cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo = SNPInfo, data =pheno1)
#' out <- skatMeta(cohort1, SNPInfo = SNPInfo)
#' head(out)
#' 
#' ####################
#' ####survival data
#' cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo = SNPInfo, data =pheno1)
#' out <- skatMeta(cohort1, SNPInfo = SNPInfo)
#' head(out)
#' }
#' @name prepScores
#' @export
prepScores <- function(Z, formula, family=stats::gaussian(), SNPInfo=NULL, snpNames="Name", aggregateBy="gene", kins=NULL, sparse=TRUE, data=parent.frame(), verbose=FALSE) {
	#fit Null model
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	} else {
	  SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy)
	}
	

	if(!is.null(kins)){
		n <- dim(kins)[1]
		stopifnot(nrow(kins) == nrow(Z))
		stopifnot(family$family == "gaussian")
		if(sparse){
			kins[kins < 2 * 2^{-6}] <- 0
			kins <- Matrix::forceSymmetric(kins)
			#oo <- order(getcc(kins))
		} else {
			#oo <- 1:n
		}
		
		#ioo <- match(1:n,oo)
		
		if(is.null(colnames(kins))){
			data$id <- 1:ncol(kins)
		} else {
			data$id <- colnames(kins)
		}
		
		nullmodel <- coxme::lmekin(formula=update(formula, '~.+ (1|id)'), data=data, varlist = 2*kins,method="REML")	
		nullmodel$theta <- c(nullmodel$vcoef$id, nullmodel$sigma^2)
	
		SIGMA <- nullmodel$theta[1]*2*kins+nullmodel$theta[2]*Diagonal(n)
		X1 <- stats::model.matrix(stats::lm(formula, data=data))
	
		s2 <- sum(nullmodel$theta)
		Om_i <- solve(SIGMA/s2)
	
		#rotate data:
		res <- as.vector(nullmodel$res)* s2 / nullmodel$theta[2]	
		nullmodel$family$var <- function(x){1}
	} else {
		nullmodel <- stats::glm(formula=formula, family = family, data=data)
		res <- stats::residuals(nullmodel, type = "response")
		X1 <- stats::model.matrix(nullmodel)
		n <- nrow(X1)
	}
 
	env <- environment()
	##check format:
	invisible(check_format_skat(Z, SNPInfo, nullmodel,aggregateBy, snpNames, formula))
	
	##match snps in Z with master list in SNPInfo file 
	mysnps <- colnames(Z)
	
	SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
	
	SItoZ <- which(colnames(Z) %in% SNPInfo[,snpNames])

	which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
	ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])

	nsnps <- sum(!is.na(ZtoSI)) #sum(which.snps.Z)
	if(nsnps == 0){ 
		stop("no column names in Z match SNP names in the SNP Info file!")
	}
	
	if(verbose){
    	cat("\n Scoring... Progress:\n")
    	pb <- utils::txtProgressBar(min = 0, max = nsnps, style = 3)
    	pb.i <- 0
    }
	
	##fit individual betas/se's
	maf0 <- colMeans(Z,na.rm=TRUE)[which.snps.Z]/2
  maf0[is.nan(maf0)] <- -1
  
	maf <- maf0[ZtoSI]
	names(maf) <- SNPInfo[,snpNames]
	
	scores <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
		if(any(is.na(z))){
      if(all(is.na(z))) z <- rep(0,length(z))
      mz <- mean(z, na.rm=TRUE)
			z[is.na(z)] <- mz
		}
        if (verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(nsnps/100) == 0) utils::setTxtProgressBar(get("pb",env),get("pb.i",env))
			}
		sum(res*z)
		})[ZtoSI]
	scores[is.na(scores)] <- 0
	names(scores) <- SNPInfo[,snpNames]
	
	if(verbose) close(pb)
	
	#deal with monomorphic SNPs (Fix issue #2 bd)
	monos <- monomorphic_snps(Z)
	scores[names(scores) %in% monos] <- 0
# 	scores[maf == 0] <- 0
	
	#differentiate missing from monomorphic:
	maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1

	#split into genes
	scores 	<- 	split(scores, SNPInfo[,aggregateBy])
	maf 	<- 	split(maf, SNPInfo[,aggregateBy])
	
	##get matrices for projection
	X1 <- sqrt(nullmodel$family$var(nullmodel$fitted))*X1
	if(is.null(kins)){
    AX1 <- with(svd(X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(u[, d > 0,drop=FALSE]))) 
	} else {
	  AX1 <- with(svd(t(X1)%*%Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%t(X1)%*%Om_i
	}
	ngenes <- length(unique(SNPInfo[,aggregateBy]))
	if(verbose){
    	cat("\n Calculating covariance... Progress:\n")
    	pb <- utils::txtProgressBar(min = 0, max = ngenes, style = 3)
    	pb.i <- 0
    }
	##get covariance matrices:
	re <- tapply(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
		inds <- match(snp.names,colnames(Z))
		mcov <- matrix(0,length(snp.names),length(snp.names))
		if(length(stats::na.omit(inds)) > 0){
		  Z0 <- as.matrix(Z[,stats::na.omit(inds),drop=FALSE])
#			Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*as.matrix(Z[,na.omit(inds),drop=FALSE])
			if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
			  if(all(is.na(z))) z <- rep(0,length(z))
       			mz <- mean(z, na.rm=TRUE)
				z[is.na(z)] <- mz
				z
			})
			Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*Z0
			if(!is.null(kins)){
				mcov[!is.na(inds), !is.na(inds)] <- as.matrix(t(Z0)%*%Om_i%*%Z0 - (t(Z0)%*%Om_i%*%X1)%*%(AX1%*%Z0))
			} else {
				mcov[!is.na(inds), !is.na(inds)] <- crossprod(Z0) - (t(Z0)%*%X1)%*%(AX1%*%Z0)
			}
		}
		rownames(mcov) <- colnames(mcov) <- snp.names
		if(verbose){
				assign("pb.i", get("pb.i",env)+1,env)
				if(get("pb.i", env)%%ceiling(ngenes/100) == 0) utils::setTxtProgressBar(get("pb",env),get("pb.i",env))		  
		}
    # set monomorphic inds to 0
    mono_snps <- intersect(snp.names, monos)
		mcov[mono_snps , ] <- 0
		mcov[ , mono_snps] <- 0
		return(Matrix::forceSymmetric(Matrix(mcov,sparse=TRUE)))
	},simplify = FALSE)
	sey = sqrt(stats::var(res)*(nrow(X1)-1)/(nrow(X1)-ncol(X1)) )
	if(family$family == "binomial") sey = 1
	if(!is.null(kins)) 	sey = sqrt(s2)

	##aggregate
	for(k in 1:length(re)){
		re[[k]] <- list("scores" = scores[[k]], "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = sey ) 
	}
	if(verbose) close(pb)
	
	attr(re,"family") <-  family$family
	class(re) <- "seqMeta"
	return(re)
}


### internal function to get connected components (pedigrees) of a sparse kinship matrix - reordering to get a block diagonal matrix can improve speed dramatically

getcc <- function(M){
	membership <- rep(0,ncol(M))
	clust <- 1
	ii <- 1
	visited <- 1
	tovisit <- NULL
	membership[ii] <- clust
	
	while(any(membership == 0)){
		ne <- which(M[ii,] != 0 & membership != clust)
		if(length(ne) > 0){
			membership[ne] <- clust
			tovisit <- union(tovisit,ne)
		}
		if(length(tovisit) >= 1){
			ii <- tovisit[1]
			tovisit <- tovisit[-1]
		} else {
			clust <- clust+1
			ii <- which(membership == 0)[1]
		}
	}	
	return(membership)
}

#' @rdname prepScores
#' @param male For analyzing the X chromosome, with prepScoresX, `male' is the
#'   gender (0/1 or F/T) indicating female/male. See details.
#' @export
prepScoresX <- function(Z, formula, male, family = stats::gaussian(), SNPInfo=NULL, snpNames = "Name", aggregateBy = "gene", kins = NULL, sparse= TRUE, data=parent.frame(), verbose = FALSE){
  #fit Null model
  if(is.null(SNPInfo)){ 
    warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
    load(paste(find.package("skatMeta"), "data", "SNPInfo.rda",sep = "/"))
    aggregateBy = "SKATgene"
  } else {
    SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy)
  }
  
  
  cl <- match.call()
  
  
  if(!is.null(kins)){
    n <- dim(kins)[1]
    stopifnot(nrow(kins) == nrow(Z))
    stopifnot(family$family == "gaussian")
    if(sparse){
      kins[kins < 2 * 2^{-6}] <- 0
      kins <- Matrix::forceSymmetric(kins)
      #oo <- order(getcc(kins))
    } else {
      #oo <- 1:n
    }
    
    #ioo <- match(1:n,oo)
    
    if(is.null(colnames(kins))){
      data$id <- 1:ncol(kins)
    } else {
      data$id <- colnames(kins)
    }
    
    nullmodel <- coxme::lmekin(formula=update(formula, '~.+ (1|id)'), data=data, varlist = 2*kins,method="REML")	
    nullmodel$theta <- c(nullmodel$vcoef$id, nullmodel$sigma^2)
    
    SIGMA <- nullmodel$theta[1]*2*kins+nullmodel$theta[2]*Diagonal(n)
    X1 <- stats::model.matrix(stats::lm(formula, data=data))
    
    s2 <- sum(nullmodel$theta)
    Om_i <- solve(SIGMA/s2)
    
    #rotate data:
    res <- as.vector(nullmodel$res)* s2 / nullmodel$theta[2]	
    nullmodel$family$var <- function(x){1}
  } else {
    nullmodel <- stats::glm(formula=formula, family = family, data=data)
    res <- stats::residuals(nullmodel, type = "response")
    X1 <- stats::model.matrix(nullmodel)
    n <- nrow(X1)
  }

  env <- environment()
  ##check format:
  invisible(check_format_skat(Z, SNPInfo, nullmodel,aggregateBy, snpNames, formula))
  
  male <- eval(cl$male,data)
  if(length(male) != length(res)) stop("`male' not the same length as phenotype")
  if(!all(male %in% c(0,1))) stop("`male' must be coded as 0/1 or T/F")
  male <- as.logical(male)
  
  ##match snps in Z with master list in SNPInfo file 
  mysnps <- colnames(Z)
  
  SNPInfo[,aggregateBy] <- as.character(SNPInfo[,aggregateBy])
  
  SItoZ <- which(colnames(Z) %in% SNPInfo[,snpNames])
  
  which.snps.Z <- colnames(Z) %in% SNPInfo[,snpNames]
  ZtoSI <- match(SNPInfo[,snpNames], mysnps[which.snps.Z])
  
  nsnps <- sum(!is.na(ZtoSI)) #sum(which.snps.Z)
  if(nsnps == 0){ 
    stop("no column names in Z match SNP names in the SNP Info file!")
  }
  
  if(verbose){
    cat("\n Scoring... Progress:\n")
    pb <- utils::txtProgressBar(min = 0, max = nsnps, style = 3)
    pb.i <- 0
  }
  
  ##calculate af, separated by gender
  maf0 <- (colSums(Z[male,which.snps.Z],na.rm=TRUE)/2 + colSums(Z[!male,which.snps.Z],na.rm=TRUE))/
    (colSums(!is.na(Z[male,which.snps.Z])) + 2*colSums(!is.na(Z[!male,which.snps.Z])))

  maf0[is.nan(maf0)] <- -1
  
  maf <- maf0[ZtoSI]
  names(maf) <- SNPInfo[,snpNames]
  
  scores <- apply(Z[,which.snps.Z, drop = FALSE],2,function(z){
    naz <- is.na(z)
    if(any(naz)){
      if(all(naz)){ 
        z <- rep(0,length(z))
      } else {
        mz <- (sum(z[male], na.rm=TRUE)/2+ sum(z[!male],na.rm=TRUE))/
        sum( (2-as.numeric(male))[!naz])
        z[naz] <- 2*mz
      }
    }
    if (verbose){
      assign("pb.i", get("pb.i",env)+1,env)
      if(get("pb.i", env)%%ceiling(nsnps/100) == 0) utils::setTxtProgressBar(get("pb",env),get("pb.i",env))
    }
    sum(res*z)
  })[ZtoSI]
  scores[is.na(scores)] <- 0
  names(scores) <- SNPInfo[,snpNames]
  
  if(verbose) close(pb)
  
  #deal with monomorphic SNPs
  scores[maf == 0] <- 0
  
  #differentiate missing from monomorphic:
  maf[!(SNPInfo[,snpNames] %in% colnames(Z))] <- -1
  
  #split into genes
  scores 	<- 	split(scores, SNPInfo[,aggregateBy])
  maf 	<- 	split(maf, SNPInfo[,aggregateBy])
  
  ##get matrices for projection
  X1 <- sqrt(nullmodel$family$var(nullmodel$fitted))*X1
  if(is.null(kins)){
    AX1 <- with(svd(X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(u[, d > 0,drop=FALSE]))) 
  } else {
    AX1 <- with(svd(t(X1)%*%Om_i%*%X1),  v[,d > 0,drop=FALSE]%*%( (1/d[d>0])*t(v[, d > 0,drop=FALSE])))%*%t(X1)%*%Om_i
  }
  ngenes <- length(unique(SNPInfo[,aggregateBy]))
  if(verbose){
    cat("\n Calculating covariance... Progress:\n")
    pb <- utils::txtProgressBar(min = 0, max = ngenes, style = 3)
    pb.i <- 0
  }
  ##get covariance matrices:
  re <- tapply(SNPInfo[,snpNames], SNPInfo[,aggregateBy],function(snp.names){
    inds <- match(snp.names,colnames(Z))
    mcov <- matrix(0,length(snp.names),length(snp.names))
    if(length(stats::na.omit(inds)) > 0){
      Z0 <- as.matrix(Z[,stats::na.omit(inds),drop=FALSE])
#       Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*as.matrix(Z[,na.omit(inds),drop=FALSE])
      if(any(is.na(Z0))) Z0 <- apply(Z0,2,function(z){
        naz <- is.na(z)
        if(any(naz)){
          if(all(naz)){ 
            z <- rep(0,length(z))
          } else {
            mz <- (sum(z[male], na.rm=TRUE)/2+ sum(z[!male],na.rm=TRUE))/
              sum( (2-as.numeric(male))[!naz])
            z[naz] <- 2*mz
          }
        }
        z
      })
      Z0 <- sqrt(nullmodel$family$var(nullmodel$fitted))*Z0
      if(!is.null(kins)){
        mcov[!is.na(inds), !is.na(inds)] <- as.matrix(t(Z0)%*%Om_i%*%Z0 - (t(Z0)%*%Om_i%*%X1)%*%(AX1%*%Z0))
      } else {
        mcov[!is.na(inds), !is.na(inds)] <- crossprod(Z0) - (t(Z0)%*%X1)%*%(AX1%*%Z0)
      }
    }
    rownames(mcov) <- colnames(mcov) <- snp.names
    if(verbose){
      assign("pb.i", get("pb.i",env)+1,env)
      if(get("pb.i", env)%%ceiling(ngenes/100) == 0) utils::setTxtProgressBar(get("pb",env),get("pb.i",env))		  
    }
    return(Matrix::forceSymmetric(Matrix(mcov,sparse=TRUE)))
  },simplify = FALSE)
  sey = sqrt(stats::var(res)*(nrow(X1)-1)/(nrow(X1)-ncol(X1)) )
  if(family$family == "binomial") sey = 1
  if(!is.null(kins)) 	sey = sqrt(s2)
  
  ##aggregate
  for(k in 1:length(re)){
    re[[k]] <- list("scores" = scores[[k]], "cov" = re[[k]], "n" =n, "maf" = maf[[k]], "sey" = sey ) 
  }
  if(verbose) close(pb)
  
  attr(re,"family") <-  family$family
  class(re) <- "seqMeta"
  return(re)
}
