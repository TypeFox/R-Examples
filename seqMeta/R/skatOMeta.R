#' @title Combine SKAT-O analyses from one or more studies.
#'   
#' @description Takes as input `seqMeta` objects (from e.g.
#'   \code{\link{prepScores}}), and meta analyzes them, using SKAT-O. See the
#'   package vignette for more extensive documentation.
#'
#' @inheritParams singlesnpMeta
#' @inheritParams burdenMeta
#' @param skat.wts Either a function to calculate testing weights for SKAT, or a
#'   character specifying a vector of weights in the SNPInfo file. For skatOMeta
#'   the default are the `beta' weights.
#' @param burden.wts Either a function to calculate weights for the burden test,
#'   or a character specifying a vector of weights in the SNPInfo file. For
#'   skatOMeta the default are the T1 weights.
#' @param rho A sequence of values that specify combinations of SKAT and a burden test to be considered. Default is c(0,1), which considers SKAT and a burden test.
#' @param method p-value calculation method. Should be one of 'saddlepoint', 'integration', or 'liu'.
#' 
#' @details \code{skatOMeta()} implements the SKAT-Optimal test, which picks the
#'   `best' combination of SKAT and a burden test, and then corrects for the
#'   flexibility afforded by this choice. Specifically, if the SKAT statistic is
#'   Q1, and the squared score for a burden test is Q2, SKAT-O considers tests
#'   of the form (1-rho)*Q1 + rho*Q2, where rho between 0 and 1. The values of
#'   rho are specified by the user using the argument \code{rho}. In the
#'   simplest form, which is the default, SKAT-O computes a SKAT test and a T1
#'   test, and reports the minimum p-value, corrected for multiple testing. See
#'   the vignette or the accompanying references for more details.
#'   
#'   If there is a single variant in the gene, or the burden test is undefined
#'   (e.g. there are no rare alleles for the T1 test), SKAT is reported (i.e.
#'   rho=0).
#'   
#'   Note 1: the SKAT package uses the same weights for both SKAT and the burden
#'   test, which this function does not.
#'   
#'   Note 2: all studies must use coordinated SNP Info files - that is, the SNP
#'   names and gene definitions must be the same.
#'   
#'   Note 3: The method of p-value calculation is much more important here than
#'   in SKAT.  The `integration' method is fast and typically accurate for
#'   p-values larger than 1e-9. The saddlepoint method is slower, but has higher
#'   relative accuracy.
#'   
#'   Note 4: Since p-value calculation can be slow for SKAT-O, and less accurate
#'   for small p-values, a reasonable alternative would be to first calculate
#'   SKAT and a burden test, and record the minimum p-value, which is a lower
#'   bound for the SKAT-O p-value. This can be done quickly and accurately.
#'   Then, one would only need to perform SKAT-O on the small subset of genes
#'   that are potentially interesting.
#'   
#'   Please see the package vignette for more details.
#'   
#' @return a data frame with the following columns:
#'   \item{gene}{Name of the gene or unit of aggregation being meta analyzed}
#'   \item{p}{p-value of the SKAT-O test.}
#'   \item{pmin}{The minimum of the p-values considered by SKAT-O (not corrected for multiple testing!).}
#'   \item{rho}{The value of rho which gave the smallest p-value.}
#'   \item{cmaf}{The cumulative minor allele frequency.}
#'   \item{nmiss}{The number of `missing` SNPs. For a gene with a single SNP
#'   this is the number of individuals which do not contribute to the analysis,
#'   due to studies that did not report results for that SNP. For a gene with
#'   multiple SNPs, is totalled over the gene. }
#'   \item{nsnps}{The number of SNPs in the gene.}
#'   \item{errflag}{An indicator of possible error: 0 suggests no error, > 0
#'   indicates probable loss of accuracy.}
#'   
#' @references Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X.
#'   (2011) Rare Variant Association Testing for Sequencing Data Using the
#'   Sequence Kernel Association Test (SKAT). American Journal of Human
#'   Genetics.
#'   
#'   Lee, S. and Wu, M.C. and Lin, X. (2012) Optimal tests for rare variant
#'   effects in sequencing association studies. Biostatistics.
#'   
#' @author Arie Voorman, Jennifer Brody
#' @seealso 
#' \code{\link{skatOMeta}}
#' \code{\link{prepScores}}
#' \code{\link{burdenMeta}} 
#' \code{\link{singlesnpMeta}}
#' 
#' @examples 
#' \dontrun{
#' ### load example data for 2 studies
#' data(seqMetaExample)
#' 
#' ####run on each study:
#' cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, data =pheno1)
#' cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, kins=kins, data=pheno2)
#' 
#' #### combine results:
#' ##skat-O with default settings:
#' out1 <- skatOMeta(cohort1, cohort2, SNPInfo = SNPInfo, method = "int")
#' head(out1)
#' 
#' ##skat-O, using a large number of combinations between SKAT and T1 tests:
#' out2 <- skatOMeta(cohort1, cohort2, rho=seq(0,1,length=11), SNPInfo=SNPInfo, method="int")
#' head(out2)
#' 
#' #rho = 0 indicates SKAT gave the smaller p-value (or the T1 is undefined) 
#' #rho=1 indicates the burden test was chosen
#' # 0 < rho < 1 indicates some other value was chosen
#' #notice that most of the time either the SKAT or T1 is chosen
#' table(out2$rho)
#' 
#' ##skat-O with beta-weights used in the burden test:
#' out3 <- skatOMeta(cohort1,cohort2, burden.wts = function(maf){dbeta(maf,1,25) }, 
#'                   rho=seq(0,1,length=11),SNPInfo = SNPInfo, method="int")
#' head(out3)
#' table(out3$rho)
#' 
#' ########################
#' ####binary data
#' cohort1 <- prepScores(Z=Z1, ybin~1, family=binomial(), SNPInfo=SNPInfo, data=pheno1)
#' out.bin <- skatOMeta(cohort1, SNPInfo = SNPInfo, method="int")
#' head(out.bin)
#' 
#' ####################
#' ####survival data
#' cohort1 <- prepCox(Z=Z1, Surv(time,status)~strata(sex)+bmi, SNPInfo=SNPInfo, 
#'                    data=pheno1)
#' out.surv <- skatOMeta(cohort1, SNPInfo = SNPInfo, method="int")
#' head(out.surv)
#' 
#' ##########################################
#' ###Compare with SKAT and T1 tests on their own:
#' cohort1 <- prepScores(Z=Z1, y~sex+bmi, SNPInfo=SNPInfo, data=pheno1)
#' cohort2 <- prepScores(Z=Z2, y~sex+bmi, SNPInfo=SNPInfo, kins=kins, data=pheno2)
#' 
#' out.skat <- skatMeta(cohort1,cohort2,SNPInfo=SNPInfo)
#' out.t1 <- burdenMeta(cohort1,cohort2, wts= function(maf){as.numeric(maf <= 0.01)}, 
#'                      SNPInfo=SNPInfo)
#'            
#' #plot results 
#' #We compare the minimum p-value of SKAT and T1, adjusting for multiple tests 
#' #using the Sidak correction, to that of SKAT-O.
#' 
#' par(mfrow=c(1,3))
#' pseq <- seq(0,1,length=100)
#' plot(y=out.skat$p, x=out1$p,xlab="SKAT-O p-value", ylab="SKAT p-value", main ="SKAT-O vs SKAT")
#' lines(y=pseq,x=1-(1-pseq)^2,col=2,lty=2, lwd=2)
#' abline(0,1)
#' 
#' plot(y=out.t1$p, x=out1$p,xlab="SKAT-O p-value", ylab="T1 p-value", main ="SKAT-O vs T1")
#' lines(y=pseq,x=1-(1-pseq)^2,col=2,lty=2, lwd=2)
#' abline(0,1)
#' 
#' plot(y=pmin(out.t1$p, out.skat$p,na.rm=T), x=out1$p,xlab="SKAT-O p-value", 
#'      ylab="min(T1,SKAT) p-value", main ="min(T1,SKAT) vs SKAT-O")	
#' lines(y=pseq,x=1-(1-pseq)^2,col=2,lty=2, lwd=2)
#' abline(0,1)
#' legend("bottomright", lwd=2,lty=2,col=2,legend="Bonferroni correction")	
#' }
#' 
#' @export
skatOMeta <- function(..., SNPInfo=NULL, skat.wts=function(maf){stats::dbeta(maf,1,25)}, burden.wts=function(maf){as.numeric(maf <= 0.01) }, rho=c(0,1), method=c("integration", "saddlepoint", "liu"), snpNames="Name", aggregateBy="gene", mafRange=c(0,0.5), verbose=FALSE) {
	cl <- match.call(expand.dots = FALSE)
	if(is.null(SNPInfo)){ 
		warning("No SNP Info file provided: loading the Illumina HumanExome BeadChip. See ?SNPInfo for more details")
		load(paste(find.package("seqMeta"), "data", "SNPInfo.rda",sep = "/"))
		aggregateBy = "SKATgene"
	} else {
	  SNPInfo <- prepSNPInfo(SNPInfo, snpNames, aggregateBy, wt1=skat.wts, wt2=burden.wts)
	}
	
	
	if(any(rho >1 | rho < 0 ) ) stop("rho must be between 0 and 1")
	method <- match.arg(method)
   #if( !(method %in% c("davies","farebrother","imhof","liu")) ) stop("Method specified is not valid! See documentation")
	
	genelist <- stats::na.omit(unique(SNPInfo[,aggregateBy]))
	cohortNames <- lapply(cl[[2]],as.character)
	ncohort <- length(cohortNames)
	
	ev <- parent.frame()
	classes <- unlist(lapply(cohortNames,function(name){class(get(name,envir=ev))})) 	
  if(!all(classes == "seqMeta" | classes == "skatCohort") ){
	 	stop("an argument to ... is not a seqMeta object!")
	}
		
	res.strings <- data.frame("gene"=genelist,stringsAsFactors=F)
	res.numeric <- matrix(NA, nrow= nrow(res.strings),ncol =  length(c("p","pmin","rho","cmaf","nmiss", "nsnps", "errflag")))
	colnames(res.numeric) <- c("p","pmin","rho","cmaf","nmiss", "nsnps","errflag")		
	
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
		
		mscores <- maf <- numeric(nsnps.sub)
		big.cov <- Matrix(0, nsnps.sub,nsnps.sub)
		n.total <- numeric(nsnps.sub)
		n.miss <- numeric(nsnps.sub)
		
		vary.ave <- 0
		for(cohort.k in 1:ncohort){
			cohort.gene <- get(cohortNames[[cohort.k]],envir=ev)[[gene]]
			
			if(!is.null(cohort.gene)){
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
					
					n.total[cohort.gene$maf >= 0] <- n.total[cohort.gene$maf >= 0]+cohort.gene$n
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
			maf <- sapply(maf, function(x){min(x,1-x)})
		
			if( !all(mafRange == c(0,0.5))){
				keep <- (maf >= min(mafRange)) & (maf <= max(mafRange))
				
				big.cov <- big.cov[keep,keep]
				mscores <- mscores[keep]
				maf <- maf[keep]
			}
		}
		if(length(maf)> 0){
		  if(is.function(skat.wts)){
		    w1 <- skat.wts(maf)
		  } else if(is.character(skat.wts)){
		    w1 <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,skat.wts])
		  } else {
		    w1 <- rep(1,length(maf))
		  }
		  
		  if(is.function(burden.wts)){
		    w2 <- burden.wts(maf)
		  } else if(is.character(burden.wts)){
		    w2 <- as.numeric(SNPInfo[SNPInfo[,aggregateBy]==gene,burden.wts])
		  } else {
		    w2 <- rep(1,length(maf))
		  }
		  
		  w1 <- ifelse(maf >0, w1,0)
		  w2 <- ifelse(maf >0, w2,0)
		  
		  ##
		  Q.skat <- sum((w1*mscores)^2, na.rm=TRUE)
		  V.skat <- (w1)*t(t(big.cov)*as.vector(w1))
		  
		  Q.burden <- sum(w2*mscores, na.rm=TRUE)^2
		  V.burden <- as.numeric(t(w2)%*%big.cov%*%w2)
		  
		  #If burden test is 0, or only 1 SNP in the gene, do SKAT:
		  if(sum(maf > 0) ==1 | V.burden ==0){
		    lambda <- eigen(zapsmall(V.skat), symmetric = TRUE)$values
		    if(any(lambda > 0) & length(lambda) >1) {
		      tmpP <- pchisqsum2(Q.skat,lambda=lambda,method=method, acc=1e-7)
		      if(tmpP$errflag !=0 ){ 
		        res.numeric[ri,"errflag"] = 1
		      } else {
		        res.numeric[ri,"errflag"] = 0
		      }
		      p <- tmpP$p
		    } else {
		      p <- ifelse(length(lambda) == 1 & all(lambda > 0), stats::pchisq(Q.skat/lambda,df=1,lower.tail=FALSE),1)
		      res.numeric[ri,"errflag"] = 0
		    }
		    res.numeric[ri,"pmin"] = res.numeric[ri,"p"] = p
		    res.numeric[ri,"rho"] = 0
		    
		    
		    #Else do SKAT-O
		  } else {
		    skato.res <- skatO_getp(mscores, big.cov, diag(w1), w2, rho, method= method, gene=gene)
		    
		    res.numeric[ri,"p"] <- skato.res$actualp
		    res.numeric[ri,"pmin"] = skato.res$minp
		    res.numeric[ri,"rho"] = skato.res$rho    	
		    res.numeric[ri, "errflag"] = skato.res$errflag
		  }
		} else {
		  res.numeric[ri,"p"] <- res.numeric[ri,"pmin"] <- 1
		  res.numeric[ri,"rho"] <- 0
      res.numeric[ri, "errflag"] <- 0
		}			
		res.numeric[ri,"cmaf"] = sum(maf,na.rm=TRUE)
		res.numeric[ri,"nsnps"] = sum(maf!= 0, na.rm =T)
		res.numeric[ri,"nmiss"] = sum(n.miss, na.rm =T)
		if(verbose){
			pb.i <- pb.i+1
			utils::setTxtProgressBar(pb, pb.i)
		}
	}
	if(verbose) close(pb)
	return(cbind(res.strings,res.numeric))
}


skatO_getp <- function(U,V, R, w, rho,method = "davies", gene=NULL){
	##Input:
	#U: score vector (length p)
	#R: p x p weight matrix for skat
	#w: burden weights
	#rho: vector of rhos in [0,1] 
	#method: method for calculating Normal quadratic form distribution
	#gene: The name of the region - used for error reporting
	
	##Output: a list with elemeents
	#minp: the minimum p-value
	#actualp: the actual p-value
	#rho: the value of rho which gave the minp
	#ps: the whole vector of p-values
	#errflag: 0 if no problem, 1 if quantile issue, 2 if integration issue
	
	satterthwaite <- function(a, df) {
        if (any(df > 1)) {
        	a <- rep(a, df)
        }
        tr <- mean(a)
        tr2 <- mean(a^2)/(tr^2)
        list(scale = tr * tr2, df = length(a)/tr2)
    }
	
	errflag = 0
	Q.skat <- crossprod(R%*%U) # SKAT 
	Q.burden <- (t(w)%*%U)^2 # burden
	Qs <- (1-rho)*Q.skat + rho*Q.burden
			
	lambdas <- ps <- NULL
	ps <- numeric(length(rho))
	for(i in 1:length(rho)){
		PC <- eigen((1-rho[i])*crossprod(R)+ rho[i]*outer(w,w),symmetric=TRUE)
		v.sqrt <- with(PC,{ values[values < 0] <- 0; (vectors)%*%diag(sqrt(values))%*%t(vectors) })
		lam <- eigen( zapsmall(v.sqrt%*%V%*%v.sqrt),only.values=TRUE,symmetric=TRUE)$values
		lam <- lam[lam != 0]
			
		lambdas <- c(lambdas, list( lam ))
		tmpP <- pchisqsum2(Qs[i],lambda=lambdas[[i]],method=method, acc=1e-7)
		if(tmpP$errflag != 0){
			errflag <- 1
			ps[i] <- pchisqsum2(Qs[i],lambda=lambdas[[i]],method="liu")$p
		} else {
			ps[i] <- tmpP$p	
		}
	}
			
	minp <- min(ps)
	Ts <- numeric(length(rho))
			
	for(i in 1:length(rho)){					
		sat <- satterthwaite(lambdas[[i]],rep(1,length(lambdas[[i]])))
		upper <- stats::qchisq(minp/20,df=sat$df,lower.tail=FALSE)*sat$scale		
		tmpT <- try(stats::uniroot(function(x){pchisqsum2(x,lambda=lambdas[[i]],method=method,acc=1e-5)$p- minp }, interval=c(1e-10,upper))$root, silent = TRUE)
		if(class(tmpT) == "try-error"){
			#warning(paste0("Problem finding quantiles in gene ", gene, ", p-value may not be accurate"))
			Ts[i] <- Qs[i]
			errflag <- 2
		} else {
			Ts[i] <- tmpT
		}
	}
	
	v11 <- R%*%V%*%R
	v12 <- R%*%V%*%w	
	v22 <- as.numeric(t(w)%*%V%*%w)
	V.cond <- v11 - outer( v12, v12 )/v22
	lambda.cond <- eigen(V.cond,only.values=TRUE,symmetric=TRUE)$values
			
	EDec <- eigen(V.cond,symmetric=TRUE)
	D <- zapsmall(diag(EDec$values))
	diag(D)[zapsmall(diag(D)) > 0] <- 1/sqrt(diag(D)[zapsmall(diag(D)) > 0])
	diag(D)[diag(D) <= 0 ] <- 0
	#meanvec <- t(EDec$vectors)%*%D%*%(EDec$vectors)%*%(v12)/c(v22) 
	meanvec <- as.numeric(D%*%t(EDec$vectors)%*%(v12)/c(v22))
			
	Fcond <- function(x,method){
		pp <- qmax <- numeric(length(x))
		for(i in 1:length(x)){
				qmax[i] <- min( ( (Ts[rho !=1 ] - rho[rho != 1]*x[i])/(1-rho)[rho !=1]) )
			if(any(x[i] > Ts[rho == 1]) ){ 
				pp[i] <- 1
			} else {
				p.tmp <- pchisqsum2(qmax[i], lambda=lambda.cond, delta = meanvec^2*x[i], method = method, acc=min(minp,1e-5) )
				if(p.tmp$errflag != 0) stop("Error in integration! using Liu p-value")
				pp[i] = p.tmp$p
			}
		}			
		return(pp)
	}
			
	if(any(lambda.cond > 0)){
		integrand <- function(x){stats::dchisq(x,1)*Fcond(x*v22,method=method)}
		integral <- try(stats::integrate(Vectorize(integrand),lower=0,upper=Inf, subdivisions = 200L, rel.tol=min(minp/100,1e-4)), silent = TRUE)
		if (class(integral) == "try-error" ) {
       		integrand <- function(x){stats::dchisq(x,1)*Fcond(x*v22,method="liu")}
       		integral <- stats::integrate(Vectorize(integrand),lower=0,upper=Inf)
       		errflag <- 3
       	} else {
       		if(integral$message != "OK") errflag <- 2
       	}
    	actualp <- integral[1]$value
    } else {
    	#cat(".")
    	actualp = minp
    }	
	return(list("actualp"= actualp, "minp" = minp, "rho" = rho[which.min(ps)], "ps" = ps, "errflag" = errflag))
}


