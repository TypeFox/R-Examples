###############################################################################
# R user interface to summary test of a multivarivate linear model
# Author: Yi Wang (yi dot wang at computer dot org)
# 05-Jan-2010
###############################################################################

summary.manyglm <- function(object, resamp="pit.trap", test="wald", p.uni="none", nBoot=1000, cor.type=object$cor.type, show.cor = FALSE, show.est=FALSE, show.residuals=FALSE, symbolic.cor = FALSE, show.time=FALSE, show.warning=FALSE,... ) 
{
    tol = object$tol
    allargs <- match.call(expand.dots = FALSE)
    dots <- allargs$...
    if ("rep.seed" %in% names(dots)) rep.seed <- dots$rep.seed
    else rep.seed <- FALSE
    if ("ld.perm" %in% names(dots)) ld.perm <- dots$ld.perm
    else ld.perm <- FALSE
    if ("bootID" %in% names(dots)) bootID <- dots$bootID
    else bootID <- NULL

    if (show.time==FALSE) st=0
    else st=1

    if (show.warning==TRUE) warn=1
    else warn=0

    if (cor.type!="I" & test=="LR") {
       warning("The likelihood ratio test can only be used if correlation matrix of the abundances is is assumed to be the Identity matrix. The Wald Test will be used.")
       test <- "Wald"
    }
    if (any(class(object)=="manylm")) {
        if ( test == "LR" ) 
           return(summary.manylm(object, resamp=resamp, test="LR", p.uni=p.uni, nBoot=nBoot, cor.type=cor.type, show.cor=show.cor, show.est=show.est, show.residuals=show.residuals, symbolic.cor=symbolic.cor, tol=tol, ld.perm=ld.perm, bootID=bootID, ... ))
	else {   
	   warning("For an manylm object, only the likelihood ratio test and F test are supported. So the test option is changed to `'F''. ")
           return(summary.manylm(object, resamp=resamp, test="F", p.uni=p.uni, nBoot=nBoot, cor.type=cor.type, show.cor=show.cor, show.est=show.est, show.residuals=show.residuals, symbolic.cor=symbolic.cor, tol=tol, ld.perm=ld.perm, bootID=bootID, ... ))
	}   
    }
    else if (!any(class(object)=="manyglm"))
       stop("The function 'summary.manyglm' can only be used for a manyglm or manylm object.")
    
    nRows = nrow(object$y)
    nVars = ncol(object$y)
    nParam = ncol(object$x)
    Y <- matrix(as.integer(object$y), nrow=nRows, ncol=nVars)
    X <- as.matrix(object$x, "numeric")

    w <- object$weights
    if (is.null(w)) w  <- rep(1, times=nRows)
    else {     
        if (!is.numeric(w))  stop("'weights' must be a numeric vector")
        if (any(w < 0)) stop("negative 'weights' not allowed")
    }

    # the following values need to be converted to integer types  


	  if (object$family == "poisson") { familynum <- 1; linkfun = 0 } 
	  else if (object$family == "negative.binomial") { familynum <- 2; linkfun = 0 }
	  else if (object$family == "binomial(link=logit)") { familynum <- 3; linkfun = 0 }
	  else if (object$family == "binomial(link=cloglog)") { familynum <- 3; linkfun = 1}
	  else stop("'family' not recognised. See ?manyglm for currently available options.") 
	
  
    if(object$theta.method == "ML") methodnum <- 0
    else if (object$theta.method == "Chi2") methodnum <- 1 
    else if (object$theta.method == "PHI") methodnum <- 2 
    else stop("'method' not defined. Choose one of 'PHI', 'ML', 'Chi2' for an manyglm object") 
    if (substr(resamp,1,1)=="c") resampnum <- 0  #case
    # To exclude case resampling
    #if (resamp=="case") stop("Sorry, case resampling is not yet available.")
    else if (substr(resamp,1,4)=="resi") resampnum <- 1  # residual
    else if (resamp=="score") resampnum <- 2  # score
    else if (substr(resamp,1,4) =="perm") resampnum <- 3 # permuation
#    else if (substr(resamp,1,1) =="f") resampnum <- 4 # free permuation
    else if (substr(resamp,1,4) ==  "mont") resampnum <- 5 # montecarlo 
    else if (substr(resamp,1,3) ==  "pit") resampnum <- 8 # PIT residual bootstrap 
    else stop("'resamp' not defined. Choose one of 'case', 'resid', 'score', 'perm.resid', 'montecarlo', 'pit.trap'")

    # allows case and parametric bootstrap only for binomial regression
    if (familynum == 3 && ( (resampnum !=5) && (resampnum!=8) ) ) {
       warning("'montecarlo' or 'pit.trap' should be used for binomial regression. Resampling option is changed to 'pit.trap'.")
       resamp <- "pit.trap"
       resampnum <- 8 
    }
  
    if (substr(test,1,1) == "w") testnum <- 2 # wald
    else if (substr(test,1,1) == "s") testnum <- 3 #score
    else if (substr(test,1,1) == "L") testnum <- 4 #LR
    else stop("'test'not defined. Choose one of 'wald', 'score', 'LR' for an manyglm object.") 
 
    if (cor.type == "R") corrnum <- 0
    else if (cor.type == "I") corrnum <- 1
    else if (cor.type == "shrink") corrnum <- 2
    else stop("'cor.type' not defined. Choose one of 'I', 'R', 'shrink'")

    if (substr(p.uni,1,1) == "n"){
       pu <- 0
       calc.pj <- FALSE
    } else if(substr(p.uni,1,1) == "u"){
       pu <- 1
       calc.pj <- TRUE 
    } else if(substr(p.uni,1,1)=="a"){
       pu <- 2
       calc.pj <- TRUE
    } else if(substr(p.uni,1,1) == "s"){
       pu <- 3
       calc.pj <- TRUE
    } else
       stop("'p.uni' not defined. Choose one of 'single', 'adjusted', 'unadjusted', 'none'.")

    if (ld.perm && is.null(bootID)) {
       warning("bootID not supplied. Calc bootID on the fly (default)...")
       ld.perm <- FALSE
    }

    if (!is.null(bootID)) {
       ld.perm <- TRUE
       nBoot <- dim(bootID)[1]
       if (max(bootID)>nRows) {
          bootID <- NULL
          cat(paste("Invalid bootID -- sample id larger than no. of observations. Generate bootID matrix on the fly.","\n"))
       }
       else {
          if (resamp == "score") {
             cat(paste("Using <double> bootID from input for score resampling.","\n"))
          }
          else { # all other methods resample the matrix index 
             if (is.integer(bootID)) {
                 cat(paste("Using <int> bootID matrix from input.","\n"))
                 if (max(bootID)==nRows) # to fit the format in C i.e. (0, nObs-1)
                     bootID <- matrix(as.integer(bootID-1), nrow=nBoot, ncol=nRows)
             }
             else {
                 bootID <- NULL
                 cat(paste("Invalid bootID -- sample id for methods other than 'score' resampling should be integer numbers up to the no. of observations. Generate bootID matrix on the fly.","\n"))
            }
         }
      }
    }



    if (corrnum == 2 | resampnum==5 ) {
       # get the shrinkage estimates
       shrink.param <- c(rep(NA,nParam+2))
       tX <- matrix(1, nRows, 1)
       if ( object$cor.type=="shrink") shrink.param[1] <- object$shrink.param
       else shrink.param[1] <- ridgeParamEst(dat=Y, X=X, weights=w, only.ridge=TRUE, tol=tol)$ridgeParameter
       objH0 <- manyglm(Y~tX, family=object$family, cor.type="shrink")       
       shrink.param[2] <- objH0$shrink.param
       objH0 <- manyglm(Y~0+X[,-1],family=object$family,cor.type="shrink")
       shrink.param[3] <- objH0$shrink.param
       for (i in 2:nParam) {
           objH0 <- manyglm(Y~X[,-i], family=object$family, cor.type="shrink")
           shrink.param[i+2] <- objH0$shrink.param
       }
    }
    else if (corrnum == 0) shrink.param <- c(rep(1,nParam+2))
    else if (corrnum == 1) shrink.param <- c(rep(0,nParam+2))
    
    modelParam <- list(tol=tol, regression=familynum, link=linkfun, estimation=methodnum, stablizer=0, n=object$K, maxiter=object$maxiter, maxiter2=object$maxiter2, warning=warn)
    # note that nboot excludes the original dataset
    testParams <- list(tol=tol, nboot=nBoot, cor_type=corrnum, test_type=testnum, resamp=resampnum, reprand=rep.seed, punit=pu, showtime=st, warning=warn)
    if(is.null(object$offset)) O <- matrix(0, nrow=nRows, ncol=nVars)
    else O <- as.matrix(object$offset)
    ######## Call Summary Rcpp #########
#    val <- .Call("RtoGlmSmry", modelParam, testParams, Y, X, O, bootID, shrink.param, PACKAGE="mvabund")
    val <- RtoGlmSmry(modelParam, testParams, Y, X, O, bootID, shrink.param)

    ######## Collect Summary Values ######## 
    smry  <- list()
    rankX <- object$rank

    # residual statistics  (from the previous R codes)
    r <- as.matrix(object$residuals)
    rss <- t(r)%*% r
    rdf <- nRows - rankX   # residual rdf
    resvar <- rss/rdf 
    genVar <- det(resvar)     

    if (is.null(object$terms))
        stop("invalid 'lm' object:  no 'terms'")
    Qr <- object$qr
    if (is.null(Qr)) Qr <- qr(object$x)

    p1 <- 1:rankX
    R  <- chol2inv(Qr$qr[p1, p1, drop = FALSE])   # = inv X'X   pxp matrix 
    genVarB <- genVar*(c(R)[1+0:(rankX-1)*(rankX+1)])
    est <- object$coefficients[Qr$pivot[p1], , drop=FALSE]
    est <- cbind(genVarB, est)
    dimnames(est)[[2]][1] <-  c("Gen. Var")

    # residual correlations
    if (show.cor)  {
       nrowR <- nrow(R)
       correlation <- matrix(NA,nrow=nrowR*nrow(resvar),ncol= nrowR*nrow(resvar))
       se <- c(outer(sqrt(c(R)[1+0:(rankX-1)*(rankX+1)]),sqrt(c(resvar)[1+0:(nVars-1)*(nVars+1)] ) ))
       for (i in 1:nrow(resvar))
       for (j in 1:nrow(resvar))
           correlation[((1:nrowR)+nrowR*(i-1)),((1:nrowR) + nrowR*(j-1))]<-(R*resvar[i,j])
       correlation <- correlation / outer(se, se)
       corrnames <- rep(paste("b_", 1:rankX, sep=""), times=nVars)
       corrnames <- paste(corrnames, rep(1:nVars, each=rankX), sep="")
       colnames(correlation) <- rownames(correlation) <- corrnames 
    }
    else correlation <- NULL

    dimnames(R) <- dimnames(X)[c(2,2)]

    # test statistics
    if (!is.null(test)) {    
       testname <- paste(test,"value")
       pname <- paste("Pr(>",test,")", sep="")
    } else {
       testname <- "no test"
       pname    <- ""
    }
    responseNames <- colnames(object$y)
    explanNames <- colnames(X)

    # significance
    significance <- cbind(val$signific, val$Psignific)
    dimnames(significance) <- list(NULL, c(testname, pname)) 
    xnames <- colnames(X) 
    dimnames(significance)[[1]] <- xnames
    
    # unitvariate tests
    if (calc.pj & !is.null(test)){
       unit_signic <- t(val$unitsign)
       unit_signic.p <- t(val$Punitsign)
       rownames(unit_signic) <- rownames(unit_signic.p) <- c(responseNames)
       colnames(unit_signic) <- colnames(unit_signic.p) <- c(explanNames)
    }
    else {
       unit_signic <- NULL
       unit_signic.p <- NULL
    }

    # overall statistics 
    df.int  <- if (attr(object$terms, "intercept")) 1 else 0
    rankX <- object$rank
    overaltest <- c(val$multstat, val$Pmultstat, rankX - df.int, rdf)
    names(overaltest) <- c(testname, pname, "num df", "den df")
    if (calc.pj) {
        unit_overaltest <- cbind(val$unitmult, val$Punitmult)
        dimnames(unit_overaltest) <- list(responseNames, c(testname, pname)) 
    } 
    else unit_overaltest <- NULL

    # display flags
    smry$show.est <- show.est
    smry$show.cor <- show.cor
    smry$show.residuals <- show.residuals
    smry$symbolic.cor <- symbolic.cor

    # resampling flags
    smry$p.uni <- p.uni
    smry$test  <- test
    smry$resamp <- resamp
    smry$cor.type <- cor.type

    # parameter values
    smry$shrink.param <- shrink.param
    smry$nBoot <- nBoot 
    smry$n.bootsdone <- val$nSamp
    smry$n.iter.sign <- val$nBoot - val$nSamp
    smry$aliased <- c(is.na(object$coefficients)[,1])

    # residual stats
    smry$residuals <- object$residuals 
    smry$df <- c(rankX, rdf, NCOL(Qr$qr))
    smry$genVar <- genVar 
    smry$est <- est     # coefficient estimates, i.e. Beta
    smry$est.stderr <- object$stderr.coefficients
    smry$cov.unscaled <- R
    smry$correlation <- correlation
    # test stats
    smry$coefficients <- significance
    smry$uni.test <- unit_signic
    smry$uni.p <- unit_signic.p
    smry$statistic <- overaltest
    smry$statistic.j <- unit_overaltest

    class(smry) <- "summary.manyglm"
    return(smry)
}
setGeneric("summary")
setMethod("summary", "manyglm", summary.manyglm)
