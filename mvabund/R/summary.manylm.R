###############################################################################
# R user interface to summary test of a multivarivate linear model
# Author: Yi Wang (yi dot wang at computer dot org)
# 05-Jan-2010
###############################################################################

summary.manylm <- function(object, nBoot=1000,resamp="residual", test="F", cor.type=object$cor.type, shrink.param=NULL, p.uni="none", studentize=TRUE, R2="h", show.cor = FALSE, show.est=FALSE, show.residuals=FALSE, symbolic.cor = FALSE, tol=1.0e-6, ... ) 
{
    allargs <- match.call(expand.dots = FALSE)
    dots <- allargs$...
    if ("rep.seed" %in% names(dots)) rep.seed <- dots$rep.seed
    else rep.seed <- FALSE
    if ("ld.perm" %in% names(dots)) ld.perm <- dots$ld.perm
    else ld.perm <- FALSE
    if ("bootID" %in% names(dots)) bootID <- dots$bootID
    else bootID <- NULL

    if(!any(class(object)=="manylm"))
       stop("The function 'summary.manylm' can only be used for a manylm object.")

    ######## Construct Summary Input Arguments ########
    nRows = nrow(object$y)
    nVars = ncol(object$y)
    nParam = ncol(object$x)
    Y <- as.matrix(object$y, nrow=nRows, ncol=nVars)
    X <- as.matrix(object$x, "numeric")
    if(substr(p.uni,1,1) == "n"){
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

    # the following values need to be converted to integer types 
    if (cor.type == "I") corr <- 1
    else if (cor.type == "R") corr <- 0
    else if (cor.type == "shrink") corr <- 2
    else stop("No such correlation type.") 
    
    # To exclude case resampling
    #if (resamp=="case") stop("Sorry, case resampling is not yet available.")
    if (resamp=="case") resam <- 0
    else if (resamp == "residual") resam <- 1
    else if (resamp == "score") resam <- 2
    else if (resamp == "perm.resid") resam <- 3
    else stop("No such resampling method.") 
 
    if (test=="LR") testype <- 0
    else if (test == "F") testype <- 1 
    else stop("No such test method.") 
 
    # estimate ridge parameter if cor.type is not "shrink" when fitting the model
    w <- object$weights
    if (is.null(w)){
       ### Fit the multivariate LM.
       w   <- rep(1, times=nRows)
    }
   else {
       if (!is.numeric(w))  stop("'weights' must be a numeric vector")
       if (any(w < 0)) stop("negative 'weights' not allowed")
   }

   if (cor.type=="shrink" & is.null(shrink.param)) {
       shrink.param <- ridgeParamEst(dat=Y, X=X, weights=w, only.ridge=TRUE, doPlot=FALSE, tol=tol)$ridgeParameter
       # to simplify later computation
       if (shrink.param == 0) cor.type <- "I"
       if (shrink.param == 1) cor.type <- "R"
       if (abs(shrink.param)>1)
          stop("the absolute 'shrink.param' should be between 0 and 1")
   }
   else if (cor.type == "I") 
      shrink.param <- 1 
   else if (cor.type == "R")
      shrink.param <- 0 


    if (!is.null(bootID)) {
       nBoot<-dim(bootID)[1]
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


    if (studentize) st <- 1
    else st <- 0
    if (substr(R2,1,1) == "h") {
       rsq <- 0
       R2name <- c("Hooper's R-squared")
    }
    else if (substr(R2, 1, 1) == "v") {
       rsq <- 1
       R2name <- c("Vector's R-squared")
    }
    else 
       stop ("No such R2 method.")   

    # construct for param list      
    params <- list(tol=tol, nboot=nBoot, cor_type=corr, shrink_param=shrink.param, test_type=testype, resamp=resam, reprand=rep.seed, studentize=st, punit=pu, rsquare=rsq)

    ######## Call Summary Rcpp #########
    val <- RtoSmryCpp(params, Y, X, bootID)

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
    rownames(R) <- rownames(X)[p1]
    colnames(R) <- colnames(X)[p1]

    # test statistics
    if(!is.null(test)) {    
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
    overaltest <- c(val$multstat, val$Pmultstat, numdf = rankX - df.int, dendf = rdf)
    names(overaltest) <- c(testname, pname, "num df", "den df")
    if (calc.pj) {
        unit_overaltest <- cbind(val$unitmult, val$Punitmult) 
        dimnames(unit_overaltest) <- list(responseNames, c(testname, pname)) 
    } 
    else unit_overaltest <- NULL

    # resampling flags
    smry$p.uni <- p.uni
    smry$test  <- test
    smry$cor.type <- cor.type
    smry$resamp <- resamp
    # display flags
    smry$show.est <- show.est
    smry$show.cor <- show.cor
    smry$show.residuals <- show.residuals
    smry$symbolic.cor <- symbolic.cor
    # parameter values
    smry$shrink.param <- shrink.param
    smry$nBoot <- nBoot 
    smry$n.bootsdone <- val$nSamp
    smry$n.iter.sign <- nBoot - val$nSamp
    smry$aliased <- c(is.na(object$coef)[,1])
    # residual stats
    smry$residuals <- r
    smry$df <- c(rankX, rdf, NCOL(Qr$qr))
    smry$genVar <- genVar 
    smry$est <- est 
    smry$cov.unscaled <- R
    smry$correlation <- correlation
    # test stats
    smry$coefficients <- significance
    smry$uni.test <- unit_signic
    smry$uni.p <- unit_signic.p
    smry$R2 <- R2name
    smry$r.squared <- c(val$R2)
    smry$statistic <- overaltest
    smry$statistic.j <- unit_overaltest

    class(smry) <- "summary.manylm"
    return(smry)
}
setGeneric("summary")
setMethod("summary", "manylm", summary.manylm)
