###############################################################################
# R user interface to anova test for comparing multivariate linear models 
# Author: Yi Wang (yi dot wang at computer dot org) and David Warton
# Last modified: 24-Mar-2015
###############################################################################

anova.manyglm <- function(object, ..., resamp="pit.trap", test="LR", p.uni="none", nBoot=999, cor.type=object$cor.type, block = NULL, show.time="total", show.warning=FALSE, rep.seed=FALSE, bootID=NULL)
{
    if (cor.type!="I" & test=="LR") {
        warning("The likelihood ratio test can only be used if correlation matrix of the abundances is is assumed to be the Identity matrix. The Wald Test will be used.")
        test <- "wald"
    }
   
    if (show.time=="none") st=0
    else if (show.time=="all") st=1
    else st=2

    if (show.warning==TRUE) warn=1
    else warn=0

    if (any(class(object) == "manylm")) {
        if ( test == "LR" ) 
        return(anova.manylm(object, ..., resamp=resamp, test="LR", p.uni=p.uni, nBoot=nBoot, cor.type=cor.type, shrink.param=object$shrink.param, bootID=bootID))
        else {
        warning("For an manylm object, only the likelihood ratio test and F test are supported. So the test option is changed to `'F''. ")
        return(anova.manylm(object, resamp=resamp, test="F", p.uni=p.uni, nBoot=nBoot, cor.type=cor.type, bootID=bootID, ... ))
    }
    }   
    else if (!any(class(object)=="manyglm"))
        stop("The function 'anova.manyglm' can only be used for a manyglm or manylm object.")

    #check if any non manylm object in ...
    objects <- list(object, ...)    
    dots <- list(...)
    ndots <- length(dots)
    if (ndots>0) {
       which <- rep(TRUE, ndots+1)
       for (i in 1:ndots) {
           if (!any(class(dots[[i]])=="manyglm")){
              objectname <- names(dots[i])
              warning(paste(objectname, "is not a manyglm object nor a valid argument- removed from input, default value is used instead."))
              which[i+1] <- FALSE
           }
       }
       objects <- objects[which]
    }

    nModels = length(objects)
    nRows <- nrow(object$y)
    nVars <- ncol(object$y)
    nParam <- ncol(object$x)
    dimnam.a <- dimnames(object$y)[[2]]
    if (is.null(dimnam.a)) dimnam.a <- paste("abund", 1:nVars)

    Y <- matrix(as.integer(object$y), nrow=nRows, ncol=nVars) 
    if (is.null(Y)) {
    #      mu.eta <- object$family$mu.eta
        eta <- object$linear.predictor
        Y <- object$fitted.values + object$residuals * log(eta)      
     }

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

    if (object$theta.method == "ML") methodnum <- 0
    else if (object$theta.method == "Chi2") methodnum <- 1 
    else if (object$theta.method == "PHI") methodnum <- 2

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
    if ( (familynum==3)&&(resampnum!=5)&&(resampnum!=8) ) { 
    #    if (override==T)    
    #       warning("'montecarlo' or 'pit.trap' should be used for binomial regression. Only proceeding because override=T, but I'd rather not...")
    #    else
    #    {
        warning("'montecarlo' or 'pit.trap' should be used for binomial regression.")
    #       resamp <- "pit.trap"
    #       resampnum <- 8
    #    }
    }
    
    if (substr(test,1,1) == "w") testnum <- 2 # wald
    else if (substr(test,1,1) == "s") testnum <- 3 #score
    else if (substr(test,1,1) == "L") testnum <- 4 #LR
    else stop("'test'not defined. Choose one of 'wald', 'score', 'LR' for an manyglm object.")  

    if (resampnum==0 && testnum!=2) # case resampling and score/LR test
       warning("case resampling with score and LR tests is under development. try case resampling with wald test.")

    if (cor.type == "R") {
        corrnum <- 0
        if ( nVars > nRows ) # p>N 
           warning("number of variables is greater than number of parameters so R cannot be estimated reliably -- suggest using cor.type='shrink'.")
    }        
    else if (cor.type == "I") corrnum <- 1
    else if (cor.type == "shrink") corrnum <- 2
    else stop("'cor.type' not defined. Choose one of 'I', 'R', 'shrink'")  

    if(substr(p.uni,1,1) == "n"){
       pu <- 0
       calc.pj <- adjust.pj <- FALSE
    } else if(substr(p.uni,1,1) == "u"){
       pu <- 1
       calc.pj <- TRUE
       adjust.pj <- FALSE
    } else if(substr(p.uni,1,1)=="a"){
       pu <- 2
       calc.pj <- adjust.pj <- TRUE
    } else
       stop("'p.uni' not defined. Choose one of 'adjusted', 'unadjusted', 'none'.")

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

#DW additions
    if(is.null(block)==FALSE)
    {
      tb=table(block)
      nLevels = length(tb)
      if(any(tb!=nRows/nLevels))
      {   
        print(tb) 
        stop("Sorry, block needs to be a balanced factor - same number of rows for each level")
      }
      else
      {
        blockIDs = vector("list",nLevels)
        for(i.level in 1:nLevels)
          blockIDs[[i.level]] = which(block==names(tb)[i.level])
        unlistIDs = unlist(blockIDs) #needed to match each resampled observation with its correct location
      }
      #then each iteration...
      if(is.null(bootID)) #generate a bootID matrix if required
        samp = matrix(sample(nLevels,nLevels*nBoot,replace=TRUE),ncol=nLevels)
      else
        samp=bootID
      bootID = matrix(NA,nBoot,nRows)
      for(iBoot in 1:nBoot)
        bootID[iBoot,unlistIDs] = unlist(blockIDs[samp[iBoot,]]) #unlistIDs is needed to make sure each unlisted blockID ends up in the right place
      bootID = bootID-1 #to fit the format in C, 0 to nRows.
      cat(paste("Using block resampling...","\n"))
#      print(bootID)
    }

    # construct for param list     
    modelParam <- list(tol=object$tol, regression=familynum, link=linkfun, maxiter=object$maxiter, maxiter2=object$maxiter2, warning=warn, estimation=methodnum, stablizer=FALSE, n=object$K)
    # note that nboot excludes the original data set
    testParams <- list(tol=object$tol, nboot=nBoot, cor_type=corrnum, test_type=testnum, resamp=resampnum, reprand=rep.seed, punit=pu, showtime=st, warning=warn)
    if(is.null(object$offset)) O <- matrix(0, nrow=nRows, ncol=nVars)
    else O <- as.matrix(object$offset)
    # ANOVA
    if (nModels==1)
    {
       # test the significance of each model terms
       X <- object$x
       varseq <- object$assign
       resdev <- resdf <- NULL
       tl <- attr(object$terms, "term.labels")
       # if intercept is included
       if (attr(object$terms,"intercept")==0)
       {
          minterm = 1
          nterms = max(1, varseq)
       }
       else
       {
          minterm = 0
          nterms <- max(0, varseq)+1
          tl <- c("(Intercept)", tl)
       }
       if ( nParam==1 )
           stop("An intercept model is comparing to itself. Stopped")

       XvarIn <- matrix(ncol=nParam, nrow=nterms, 1)
       for ( i in 0:(nterms-2))
       { # exclude object itself
           XvarIn[nterms-i, varseq>i+minterm] <- 0 # in reversed order
           ncoef <- nParam-length(varseq[varseq>i+minterm])
           resdf <- c(resdf, nRows-ncoef)
       }

       resdf <- c(resdf, object$df.residual)
       # get the shrinkage estimates
       tX <- matrix(1, nrow=nRows, ncol=1)
       if (corrnum==2 | resampnum==5){ # shrinkage or montecarlo bootstrap
#          shrink.param <- c(rep(NA, nterms))
          # use a single shrinkage parameter for all models
          if (object$cor.type == "shrink") {
#                 shrink.param[1] <- object$shrink.param
              shrink.param <- rep(object$shrink.param,nterms)
          }
      else  {
#              shrink.param[1] <- ridgeParamEst(dat=object$residuals, X=tX, 
#                        only.ridge=TRUE)$ridgeParam      
  
              lambda <- ridgeParamEst(dat=object$residuals, X=tX, 
                         only.ridge=TRUE)$ridgeParam          
              shrink.param <- rep(lambda, nterms)
          }
#          for ( i in 0:(nterms-2)){ # exclude object itself
#              fit <- .Call("RtoGlm", modelParam, Y, X[,varseq<=i+minterm,drop=FALSE], 
#                 PACKAGE="mvabund")
#              shrink.param[nterms-i] <- ridgeParamEst(dat=fit$residuals, 
#                      X=tX, only.ridge=TRUE)$ridgeParam # in reversed order
#           }
       }   
       else if (corrnum == 0) shrink.param <- c(rep(1, nterms))
       else if (corrnum == 1) shrink.param <- c(rep(0, nterms))
#       resdev <- c(resdev, object$deviance) 
       nModels <- nterms
       ord <- (nterms-1):1
       topnote <- paste("Model:", deparse(object$call))
    }   
    else {
        targs <- match.call(expand.dots = FALSE)
        if (targs[[1]] == "example" || any(class(object) == "traitglm"))
            modelnamelist <- paste("Model", format(1:nModels))
        else    
            modelnamelist <- as.character(c(targs[[2]], targs[[3]]))

        resdf   <- as.numeric(sapply(objects, function(x) x$df.residual))
        ####### check input arguments #######
        # check the order of models, so that each model is tested against the next smaller one 
        ord <- order(resdf, decreasing=TRUE) 
        objects <- objects[ord]
        resdf <- resdf[ord]
        modelnamelist <- modelnamelist[ord]

        # get the shrinkage estimates
        if (corrnum == 2 | resampnum == 5) { # shrinkage or parametric bootstrap
        shrink.param <- c(rep(NA,nModels))
            tX <- matrix(1, nrow=nRows, ncol=1)
        for ( i in 1:nModels ) {
            if (objects[[i]]$cor.type == "shrink") 
                    shrink.param[i] <- objects[[i]]$shrink.param
            else shrink.param[i] <- ridgeParamEst(dat=objects[[i]]$residuals, X=tX, only.ridge=TRUE)$ridgeParam 
        }
    }
        else if (corrnum == 0) shrink.param <- c(rep(1,nModels))
        else if (corrnum == 1) shrink.param <- c(rep(0,nModels))

        # Test if models are nested, construct the full matrix and XvarIn 
#        XNull <- as.matrix(objects[[1]]$x, "numeric")
#        ind <- matrix(ncol=1, nrow=nModels)
#        for ( i in 2:nModels ) {
#            XAlt  <- as.matrix(objects[[i]]$x, "numeric")
#            Xarg  <- cbind(XAlt, XNull)
#            tmp <- qr(Xarg)
#            Xplus <- qr(XAlt)
#            if ( tmp$rank == Xplus$rank ) {
#               Beta <- qr.coef(Xplus, XNull)  # equivalent to (XAlt\XNull) in matlab 
#               # The following gets the left null space of beta, ie.LT=null(t(beta));
#               # note that LT is an orthogonal complement of Beta, and [Beta, LT] together forms the orthogonal basis that span the column space of XAlt
#               # For some reason, it must be null(beta) instead of null(t(beta)) in R to get the same answer in matlab.
#               tmp <- qr(Beta)
#               set <- if(tmp$rank == 0) 1:ncol(Beta) else  - (1:tmp$rank)
#               LT <- qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
#               # to get the dimension of Xnull
#               ind[nModels+2-i, 1] <- dim(XNull)[2]
#               XNull <- cbind(XNull, XAlt%*%LT)
#            } 
#        X <- XNull
        X <- as.matrix(objects[[nModels]]$x, "numeric") # XAlt
        nParam <- dim(X)[2] 
        XvarIn <- matrix(ncol=nParam, nrow=nModels, as.integer(0))  
        Xnames <- list()   # formula of each model
        for ( i in 1:nModels ) {
            nx = dim(as.matrix(objects[[i]]$x, "numeric"))[2]
            XvarIn[nModels+1-i, 1:nx] <- as.integer(1) 
        }
#browser()
        Xnames <- lapply(objects, function(x) paste(deparse(formula(x), 
                         width.cutoff=500), collapse = "\n")) 
        topnote <- paste(modelnamelist, ": ", Xnames, sep = "", collapse = "\n")
        tl <- modelnamelist 
        if (tl[1]==tl[2]) {
	    warning(paste("Two identical models. Second model's name changed to ", tl[2], "_2", sep=""))
	    tl[2] <- paste(tl[2], "_2", sep="")
	}
        ord <- (nModels-1):1
    }

    ######## call resampTest Rcpp #########   
#        val <- .Call("RtoGlmAnova", modelParam, testParams, Y, X, O,
 #                XvarIn, bootID, shrink.param, PACKAGE="mvabund")

    val <- RtoGlmAnova(modelParam, testParams, Y, X, O, XvarIn, bootID, shrink.param)

    # prepare output summary
    table <- data.frame(resdf, c(NA, val$dfDiff[ord]), 
                 c(NA, val$multstat[ord]), c(NA, val$Pmultstat[ord])) 
    uni.p <- matrix(ncol=nVars,nrow=nModels) 
    uni.test <- matrix(ncol=nVars, nrow=nModels)
    uni.p[2:nModels, ] <- val$Pstatj[ord,]   
    uni.test[2:nModels, ] <- val$statj[ord,]

    anova <- list()
    # Supplied arguments
    anova$family <- object$family
    anova$p.uni <- p.uni
    anova$test  <- if (test=="LR") "Dev" else test
    anova$cor.type <- cor.type
    anova$resamp <- if (resamp=="montecarlo") "parametric" else resamp
    anova$nBoot <- nBoot 
    # estimated parameters
    anova$shrink.param <- shrink.param
    anova$n.bootsdone <- val$nSamp
    # test statistics
    anova$table <- table 
    anova$uni.p <- uni.p
    anova$uni.test <- uni.test

    ########### formal displays #########
    # Title and model formulas
    title <- if (test=="LR") "Analysis of Deviance Table\n"
             else "Analysis of Variance Table\n" 
    attr(anova$table, "heading") <- c(title, topnote) 
    attr(anova$table, "title") <- "\nMultivariate test:\n" 
    # make multivariate table 
    if (!is.null(test)) {
       testname <- anova$test
       pname    <- paste("Pr(>",anova$test,")", sep="")
    } else {
       testname <- "no test"
       pname    <- ""
    }
    
    dimnames(anova$table) <- list(tl, c("Res.Df", "Df.diff", testname, pname))
   
    # make several univariate tables, if required 
    attr(anova$uni.test, "title") <- attr(anova$uni.p, "title") <- "Univariate Tests:"
    dimnames(anova$uni.p) <- dimnames(anova$uni.test) <- list(tl, dimnam.a)

    if(p.uni=="none") #hack fix because univariate P's were all coming out as 1/B for "none" (!?)
      anova$uni.p[is.numeric(anova$uni.p)] <- NA

    anova$block = block
    class(anova) <- "anova.manyglm"
    return(anova)
}
