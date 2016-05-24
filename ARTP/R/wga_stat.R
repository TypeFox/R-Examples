# History: May 02 2008 Add getWaldTest
#          May 14 2008 Make getWald test more general
#          May 30 2008 Add getMAF function
#          Jun 03 2008 Use matchNames function in getWaldTest.
#          Jun 04 2008 Add getSummary function
#          Jun 10 2008 Add setUpSummary function
#          Jun 16 2008 Rename getLoglike to getLoglike.glm
#          Jun 18 2008 Redo getWaldTest
#          Jun 30 2008 Let getWaldTest call setUpSummary
#          Jul 11 2008 Add getPermutation function
#          Jul 15 2008 In getWaldTest, use the try function
#          Jul 17 2008 Add getGenoCounts
#          Jul 18 2008 Remove missing names from getWaldTest
#                      Remove inc.Beta.se option
#          Jul 21 2008 Generalize likelihood ratio to lists
#                      Add getXBeta function
#          Jul 25 2008 Add getDesignMatrix function
#          Aug 05 2008 Add effects functions
#          Aug 08 2008 Compute standard errors for effects
#                      Add getLinearComb.var
#          Aug 12 2008 Change getDesignMatrix
#          Aug 13 2008 Add getModelData
#          Aug 13 2008 Use getModelData in callGLM
#          Aug 14 2008 Change to getPermutation function
#          Aug 21 2008 Change to effects.init function
#          Aug 25 2008 Add effects function
#          Aug 29 2008 Add waldTest.main
#          Aug 29 2008 Add getSummary.main
#          Sep 29 2008 Add getORfromLOR and getCI
#          Oct 01 2008 Fix bug in computing LR p-value
#          Oct 22 2008 Fix bug in getDesignMatrix
#          Oct 23 2008 Add pvalue.normal function
#          Nov 18 2008 Add unadjustedGLM.counts
#          Dec 05 2008 Catch errors in unadjustedGLM.counts
#          Jan 15 2009 Add GC.adj.pvalues and inflationFactor functions
#          Feb 03 2009 Change in getMAF
#          Mar 26 2009 Add swap2cols.cov function
#          Apr 04 2009 Fix bug in getDesignMatrix
#          Apr 27 2009 Modify effects.init return list
#          Apr 28 2009 Change in effects.init for the first row of a
#                      stratified effects table.
#          Jun 04 2009 Add heterTest function
#          Jun 12 2009 Add freqCounts.var function
#                      Add standardize.z function
#          Jun 19 2009 Fix bug in getDesignMatrix for interaction
#                      matrices of 1 column, and setting the colnames
#          Jul 15 2009 Add inflation factor argument to GC.Adj.Pvalues function
#          Jul 20 2009 Update freqCounts.var for left endpoint = right endpoint
#          Jul 28 2009 Update freqCounts.var for left endpoint = right endpoint
#          Jul 30 2009 Update freqCounts.var to include data frame
#          Aug 03 2009 Change output of effects.init
#          Oct 06 2009 Change levels in freqCounts.var for leftEndClosed
#          Oct 09 2009 Update freqCounts.var to pass in labels
#          Oct 23 2009 Add function to check the convergence of an object
#          Oct 26 2009 Add functions for likelihood ratio test, Wald test
#          Oct 28 2009 Add function dsgnMat
#          Dec 28 2009 Let a colon be the seperator in snp.effects for the interaction
#          Dec 29 2009 Add option removeInt to dsgnMat
#          Mar 03 2010 Update getWaldTest, getEstCov for snp.matched class
#          Mar 03 2010 Update snp.effects for the snp.matched class
#          Mar 13 2010 Add method option to getSummary, getWaldTest
#                      Compute stratified effects in snp.effects
#                      Add function to print snp.effects object
#          Mar 15 2010 Add method option to snp.effects
#          Mar 18 2010 Add function getGenoStats
#          Aug 18 2010 Add base1.name option to effects.init
#          Sep 23 2010 Update inflationFactor function
#          Oct 12 2010 Update waldTest.main 
#          Nov 09 2010 Update waldTest.main to return NA if chisq test < 0
#          Nov 29 2010 Add option to getXBeta
#                      Add getExtremeSubs function
#          Dec 23 2010 Add getOR.CI function
#          Jan 12 2011 Remove extended option in getDesingMatrix
#          Nov 01 2011 Add getPermutation.strata function
#          Feb 22 2012 Add myrmvnorm function

# Function to return point estimates and covariance matrix from an object
getEstCov <- function(fit) {

  clss <- class(fit)
  if (any(clss == "snp.logistic")) {
    methods <- c("UML", "CML", "EB")
    ret <- list(methods=methods)
    for (method in methods) {
      temp <- fit[[method, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[method]] <- list(estimates=temp$parms, cov=temp$cov)
      } 
    }
    return(ret)
  } else if (any(clss == "glm")) {
    parms <- fit$coefficients
    fit   <- summary(fit)
    cov   <- fit$cov.scaled
  } else if (any(clss == "coxph")) {
    parms         <- fit$coefficients
    cov           <- fit$var
    vnames        <- names(parms)
    rownames(cov) <- vnames
    colnames(cov) <- vnames
  } else if (any(clss == "vglm")) {
    parms <- fit@coefficients
    fit   <- summary(fit)
    cov   <- fit@cov.unscaled
  } else if (any(clss == "snp.matched")) {
    methods <- c("CLR", "CCL", "HCL")
    ret <- list(methods=methods)
    for (method in methods) {
      temp <- fit[[method, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[method]] <- list(estimates=temp$parms, cov=temp$cov)
      } 
    }
    return(ret)
  } else {
    parms <- fit$parms
    if (is.null(parms)) parms <- fit$coefficients 
    cov <- fit$cov
    if (is.null(cov)) cov <- fit$cov.scaled 
  }

  list(estimates=parms, cov=cov)

} # END: getEstCov

# Function to get the log-likelihood from a glm object
getLoglike.glm <- function(model) {
  # model

  # AIC = -2*loglike + 2*(# of parms)
  (2*model$rank - model$aic)/2 

} # END: getLoglike.glm

# Main function for likelihood ratio test
likelihoodRatio.main <- function(ll1, rank1, ll2, rank2) {

  # ll1, ll2   log-likelihood values

  df     <- abs(rank1 - rank2)
  test   <- 2*abs(ll1 - ll2)
  if (!df) {
    pvalue <- 1
  } else {
    pvalue <- pchisq(test, df, lower.tail=FALSE)
  }
  list(test=test, df=df, pvalue=pvalue)

} # END: likelihoodRatio.main

# Function to return the log-likelihood and rank of an object
loglikeAndRank <- function(fit) {

  clss <- class(fit)
  if (any(clss == "glm")) {
    rank <- fit$rank
    ll   <- (2*rank - fit$aic)/2
  } else if (any(clss == "coxph")) {
    rank <- sum(!is.na(fit$coefficients))
    ll   <- max(fit$loglik)
  } else if (any(clss == "vglm")) {
    rank <- fit@rank
    ll   <- fit@criterion$loglikelihood
  } else {
    rank <- fit$rank
    ll   <- fit$loglike
  } 

  list(loglike=ll, rank=rank)

} # END: loglikeAndRank

# Function to do a likelihood ratio test
likelihoodRatio <- function(model1, model2) {
  # model1    Return object from glm, lm, snp.logistic, coxph, vglm or list
  #           with names "loglike" and "rank"
  # model2

  l1 <- loglikeAndRank(model1)
  l2 <- loglikeAndRank(model2)

  ret <- likelihoodRatio.main(l1$loglike, l1$rank, l2$loglike, l2$rank) 

  ret

} # END: likelihoodRatio

# Function to compute the Wald test (2 - sided) 
waldTest.main <- function(parms, cov, parmNames) {

  # parms      Parameter vector
  # cov        Covariance matrix
  # parmNames  Character or numeric vector of parameters to test

  df     <- length(parmNames)
  nrcov  <- nrow(cov)
  vnames <- names(parms)

  if (is.numeric(parmNames)) {
    temp <- parmNames %in% 1:nrcov
    vpos <- parmNames[temp]
    np   <- length(vpos)
    if (!np) return(list(test=NA, df=0, pvalue=NA))

    # Update parms and cov
    parms <- parms[vpos]
    cov   <- cov[vpos, vpos] 

  } else {
    # Remove missing names in parmNames
    vnames <- vnames[vnames %in% parmNames]
     
    # Update the parameter vector (the name is kept if length(vnames) = 1)
    parms <- parms[vnames]

    # Remove missing values
    temp   <- !is.na(parms)
    parms  <- parms[temp]
    vnames <- vnames[temp]

    # Check for error
    np <- length(parms)
    if (!np) return(list(test=NA, df=0, pvalue=NA))
    
    # Update cov
    cov <- cov[vnames, vnames]
  }

  if (np == 1) {
    test   <- parms/sqrt(cov)
    pvalue <- 2*pnorm(abs(test), lower.tail=FALSE) 
    return(list(test=test, df=np, pvalue=pvalue))
  } 

  # See if matrix is invertible
  temp <- try(solve(cov), silent=TRUE)
  if (class(temp) == "try-error") {
    return(list(test=NA, df=np, pvalue=NA))
  } 

  # Get the test statistic
  dim(parms) <- c(np, 1)
  test       <- t(parms) %*% temp %*% parms
  dim(test)  <- NULL
  if (test >= 0) {
    pvalue   <- pchisq(test, df=np, lower.tail=FALSE) 
  } else {
    pvalue   <- NA
  }

  list(test=test, df=np, pvalue=pvalue)

} # END: waldTest.main

# Function to compute the Wald test (2 - sided)
getWaldTest <- function(fit, parmNames, method=NULL) {

  # fit        Return object from glm, list with names "coefficients"
  #            and "cov.scaled", return object from snp.logistic or snp.matched.
  # parmNames  Character or numeric vector of parameters to test
  # method     

  estcov <- getEstCov(fit)
  methods <- estcov[["methods", exact=TRUE]]
  if (!is.null(methods)) {
    ret <- list()
    for (m in methods) {
      temp <- estcov[[m, exact=TRUE]]
      if (!is.null(temp)) ret[[m]] <- waldTest.main(temp$estimates, temp$cov, parmNames)
    } 
    if (!is.null(method)) ret <- ret[[method, exact=TRUE]]
    return(ret)
  } 

  return(waldTest.main(estcov$estimates, estcov$cov, parmNames))    
  
} # END: getWaldTest

# Function to compute the threshold p-value
threshold.trunc.prod <- function(pvals, threshold=0.05) {

  # pvals      Vector or matrix of p-values   
  # threshold  The default is 0.05

  dim(pvals) <- NULL

  # Get the pvalues less than threshold
  pvals <- pvals[((pvals < threshold) & (!is.na(pvals)))]

  ret <- exp(sum(log(pvals)))

  ret

} # END: rankTruncate.pvals

# Function to compute rank truncated p-value
rank.trunc.prod <- function(pvals, k=10) {

  # pvals   Vector or matrix of p-values   
  # k       Maximum number of p-values to use

  # Get the sorted p-values
  pvals <- sort(pvals)
  n     <- min(length(pvals), k)
  ret   <- exp(sum(log(pvals[1:n])))

  ret

} # END: rankTruncate.pvals

# Function to compute score test for logistic reg
score.logReg <- function(fit, mat) {

  # fit      Return object from glm with x=TRUE and y=TRUE in the call
  # mat      A single factor, matrix, or data frame of variables to test

  # See if mat is a factor
  if (is.factor(mat)) mat <- data.frame(mat)
  if (is.data.frame(mat)) {
    mat <- as.matrix(createDummy(mat)$data)
  }

  # Get the number of columns of mat
  df <- ncol(mat)
  if (is.null(df)) df <- 1

  temp <- exp(fit$linear.predictors)
  p    <- temp/(1 + temp)  

  # Add the new vector to x
  x  <- cbind(fit$x, mat)
  
  # Get the gradient
  U      <- colSums(matrixMultVec(x, fit$y-p, by=2)) 
  n      <- length(U)
  dim(U) <- c(n, 1)

  # Free memory
  rm(fit, mat)
  temp <- gc(verbose=FALSE)

  # Let p = p*(1-p)
  p <- p*(1 - p)

  # Get the negative Hessian
  hess <- matrix(data=NA, nrow=n, ncol=n)
  for (i in 1:n) {
    temp      <- p*x[, i]
    hess[i, ] <- colSums(matrixMultVec(x, temp, by=2)) 
    hess[, i] <- hess[i, ]
  }

  # Invert the hessian
  hess <- try(solve(hess), silent=TRUE)
  if (class(hess) == "try-error") {
    warning("Singular hessian matrix")
    return(list(test=NA, df=df, pvalue=NA))
  } 

  # Compute the chi-squared test statistic
  test <- t(U) %*% hess %*% U
  dim(test) <- NULL

  pvalue <- pchisq(test, df=df, lower.tail=FALSE)

  list(test=test, df=df, pvalue=pvalue)

} # score.logReg

# Function to compute minor allele frequency
getMAF <- function(genotype, sub.vec=NULL, controls=0) {

  # genotype       Vector of genotypes coded as 0, 1, 2, NA
  #                No default
  # sub.vec        NULL or case/control vector for only using
  #                a subset of the genotypes.
  #                The default is NULL
  # controls       Vector of values describing the controls in 
  #                sub.vec.
  #                The default is 0.

  temp <- !is.na(genotype)
  if (!is.null(sub.vec)) temp <- temp & (sub.vec %in% controls)
  genotype <- genotype[temp]
  ng <- length(genotype)
  if (!ng) return(NA)

  freq <- (sum(genotype==1) + 2*sum(genotype==2))/(2*ng)
  freq

} # END: getMAF 

# Function to return summary information for parameters
getSummary.main <- function(parms, cov, sided=2) {

  # parms   Vector of parameters
  # cov     Covariance matrix 
  # sided   1 or 2 

  if (sided != 1) sided <- 2
  cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")

  n    <- length(parms)
  ret  <- matrix(data=NA, nrow=n, ncol=4)
  pnames <- names(parms)
  rownames(ret) <- pnames
  colnames(ret) <- cols
  ret[, 1] <- parms
  
  cols <- colnames(cov)
  cov  <- sqrt(diag(cov))
  names(cov) <- cols
  
  # Get the correct order
  if (is.null(pnames)) pnames <- 1:n
  cov <- cov[pnames]
  ret[, 2] <- cov

  ret[, 3] <- parms/cov 
  ret[, 4] <- sided*pnorm(abs(ret[, 3]), lower.tail=FALSE)
  ret

} # END: getSummary.main

# Function to return summary information for parameters
getSummary <- function(fit, sided=2, method=NULL) {

  # fit     Return object from glm, snp.logistic, or a list
  #         with names "parms" and "cov".
  # sided   1 or 2 

  clss <- class(fit)

  # snp.logistic
  if (any(clss %in% "snp.logistic")) {
    if (is.null(method)) {
      methods <- c("UML", "CML", "EB")
    } else {
      methods <- toupper(method)
    }
    ret <- list()
    for (m in methods) {
      temp <- fit[[m, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[m]] <- getSummary.main(temp$parms, temp$cov, sided=sided)
      } 
    }
    return(ret)
  } 

  # snp.matched
  if (any(clss %in% "snp.matched")) {
    if (is.null(method)) {
      methods <- c("CLR", "CCL", "HCL")
    } else {
      methods <- toupper(method)
    }
    ret <- list()
    for (m in methods) {
      temp <- fit[[m, exact=TRUE]]
      if (!is.null(temp)) {
        ret[[m]] <- getSummary.main(temp$parms, temp$cov, sided=sided)
      } 
    }
    return(ret)
  } 

  # GLM
  if ("glm" %in% clss) fit <- summary(fit)
  if (class(fit) == "summary.glm") {
    cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
    if (sided != 1) sided <- 2

    fit$coefficients[, 4] <- 
       sided*pnorm(abs(fit$coefficients[, 3]), lower.tail=FALSE)
    colnames(fit$coefficients) <- cols
    return(fit$coefficients)
  }

  # List
  parms <- fit$parms
  if (is.null(parms)) {
    parms <- fit$coefficients
    if (is.null(parms)) return(NULL)
  }
  cov <- fit$cov
  if (is.null(cov)) {
    cov <- fit$cov.scaled
    if (is.null(cov)) return(NULL)
  }
  return(getSummary.main(parms, cov, sided=sided))      

} # END: getSummary

# Function to take a parameter vector and covariance matrix and
#  output an nx2 matrix that can be used with getWaldTest
setUpSummary <- function(parms, cov) {

  cnames         <- colnames(cov)
  nc             <- ncol(cov)
  if (is.null(nc)) nc <- 1
  temp           <- matrix(data=NA, nrow=nc, ncol=2)
  colnames(temp) <- c("Estimate", "Std. Error")
  rownames(temp) <- cnames
  temp[, 2]      <- sqrt(diag(cov))
 
  # Match the parameter names (there could be NAs in the vector of
  #   point estimates)
  if ((!is.null(cnames)) && (!is.null(names(parms)))) {
    parms <- parms[cnames]
  }
  temp[, 1] <- parms
 
  temp <- list(coefficients=temp, cov.scaled=cov)
  temp
 
} # END: setUpSummary

# Function to call for permutations with a stratification variable
getPermutation.strata <- function(vec, start=NULL, stop=NULL) {

  # vec           Vector of the strata variable (example family ids)
  #               This vector MUST be sorted
  # start         Starting indices for each unique value of vec
  # stop          Stopping indices for each unique value of vec

  if ((is.null(start)) || (is.null(stop))) {
    levels  <- table(vec)
    nlevels <- length(levels)
    start   <- integer(nlevels)
    stop    <- integer(nlevels)
    a       <- 1
    for (i in 1:nlevels) {
      start[i] <- a
      b        <- a + levels[i] - 1
      stop[i]  <- b
      a        <- b + 1 
    } 
  } else {
    nlevels <- length(start)
  }

  ret <- integer(length(vec))
  for (i in 1:nlevels) {
    ids <- start[i]:stop[i]
    ret[ids] <- sample(ids, replace=FALSE)
  }

  list(perm=ret, nlevels=nlevels, start=start, stop=stop)

} # END: getPermutation.strata

# Function to call for permutations
getPermutation <- function(fit0, nsub, perm.method=1) {

  # fit0           Base model fit
  # nsub
  # perm.method    1-3

  if (perm.method == 3) {
    # For gaussian family
    # Permute the residuals
    errors <- sample(fit0$residuals)

    # Add to linear predictor
    response <- fit0$linear.predictors + errors

    perm <- 1:nsub
  } else if (perm.method == 2) {
    # For binomial family
    perm     <- 1:nsub
    response <- rbinom(nsub, 1, fit0$fitted.values)    
  } else {
    perm     <- sample(1:nsub)
    response <- fit0$y
  }

  list(response=response, perm=perm)

} # END: getPermutation

# Function to call glm
callGLM <- function(y, X.main=NULL, X.int=NULL, int.vec=NULL,
                    family="binomial", prefix="SNP_", retX=TRUE,
                    retY=TRUE, inc.int.vec=1, int.vec.base=0) {

  # y           Response vector
  # X.main      Matrix of main effects (without intercept and int.vec)
  # X.int       Matrix for interactions
  # int.vec     Interaction vector or factor 
  # family
  # inc.int.vec 0 or 1 to include the interacting vector in the model

  temp <- getModelData(y, int.vec, X.main=X.main, X.int=X.int, 
                    prefix=prefix, inc.snp=inc.int.vec, 
                   snp.base=int.vec.base)

  y <- temp$y
  X <- temp$design

  if (is.null(prefix)) {
    fit <- glm(y ~ X-1, family=family,
            model=FALSE, x=retX, y=retY) 
  } else {
    fit <- glm(y ~ .-1, family=family, data=data.frame(X),
            model=FALSE, x=retX, y=retY)
  }

  fit

} # END: callGLM

# Function to return genotype counts
getGenoCounts <- function(snp, exclude=c(NA, NaN), check=1) {

  ret <- table(snp, exclude=exclude)
  if ((check) && (length(ret) < 3)) {
    temp        <- rep.int(0, times=3)
    names(temp) <- c("0", "1", "2")
    nm          <- names(ret)
    temp[nm]    <- ret
    ret         <- temp
  }
  ret

} # END: getGenoCounts

# Function to compute XBeta (linear.predictors) by matching the names
getXBeta <- function(X, beta, drop=NULL) {

  X <- as.matrix(X)
  
  if (!is.null(drop)) {
    temp <- !(names(beta) %in% drop)
    beta <- beta[temp]
  }
  cnames  <- intersect(colnames(X), names(beta))
  print("Variables used:")
  print(cnames)
  X       <- removeOrKeepCols(X, cnames, which=1)
  if (!is.numeric(X)) {
    dimX   <- dim(X)
    temp   <- colnames(X)
    X      <- as.numeric(X)
    dim(X) <- dimX
    colnames(X) <- temp
  }

  beta    <- beta[cnames]
  b2      <- beta
  dim(b2) <- c(length(beta), 1)
  ret     <- X %*% b2
  list(X=X, beta=beta, XBeta=ret)

} # END: getXBeta

# Function to return the model data for callGLM
getModelData <- function(y, snp, X.main=NULL, X.int=NULL, 
                    prefix=NULL, inc.snp=1, snp.base=0) {

  pflag <- !is.null(prefix)
  if (pflag) cnames <- colnames(X.main)

  # Append y and intercept to X.main
  X.main <- cbind(y, 1, X.main)
  if (pflag) colnames(X.main) <- c("y", "Intercept", cnames)

  mat <- getDesignMatrix(snp, X.main=X.main, X.int=X.int,
                   inc.snp=inc.snp, X.hasInt=1, prefix=prefix,
                   snp.base=snp.base)
 
  # Remove the response from mat
  y   <- mat[, 1]
  mat <- removeOrKeepCols(mat, 1, which=-1) 
  
  list(y=y, design=mat)

} # END: getModelData

# Function to return a design matrix.
# Missing values are automatically removed
getDesignMatrix <- function(snp, X.main=NULL, X.int=NULL,
                   inc.snp=1, X.hasInt=0, prefix=NULL,
                   snp.base=0) {
  
  # snp         Vector or factor. If a factor, see snp.base
  # X.main      Matrix for main effects
  # X.int       Matrix for interactions with snp
  # inc.snp     0 or 1 to include snp in the model
  # X.hasInt    0 or 1 if X.main has an intercept column
  # prefix      NULL or snp prefix for variable names
  #             If set to NULL, the default variable names
  #             from model.matrix will be kept.
  # snp.base    Baseline category for snp that will be left
  #             out of the returned matrix

  # Input matrices should have column names !!!
  # If column names for X.int contain a colon, then there will be a problem
  # Watch for X.main = NULL

  intFlag <- !is.null(X.int)
  if (intFlag) {
    intNames   <- colnames(X.int)
    intIds     <- grep(":", intNames)
    intIdsFlag <- length(intIds)
    if (intIdsFlag) {
      stop("ERROR: X.int column names cannot contain a colon (:)")
      intNames <- gsub(":", ".", intNames)
      colnames(X.int) <- intNames
    } 
  }
  pFlag   <- !is.null(prefix)
  snpFlag <- !is.null(snp)
  if (snpFlag) {
    facFlag <- is.factor(snp)
  } else {
    facFlag <- 0
  } 
  if (!snpFlag) {
    inc.snp <- 0
    intFlag <- 0
  }
 
  # An intercept will always be included
  if (!X.hasInt) X.main <- cbind(rep(1, times=length(snp)), X.main)

  if (intFlag) {
    #colnames(X.int) <- NULL
    mat <- model.matrix(~X.main + snp*X.int - 1 - X.int)
  } else {
    if (snpFlag) {
      mat <- model.matrix(~X.main + snp - 1)
    } else {
      mat <- model.matrix(~X.main - 1)
    }
  }

  # Set the column names
  if (pFlag) {
    assign <- attributes(mat)$assign
    max    <- max(assign)
    cnames <- colnames(mat)

    # Main effects
    ids  <- assign == 1
    temp <- cnames[ids]
    # Start from pos 7 (X.main is 6 chars)
    cnames[ids] <- substring(temp, 7)

    # Intercept column
    cnames[1] <- "Intercept"

    # SNP
    if (max > 1) {
      ids  <- assign == 2
      temp <- cnames[ids]
      # Start from pos 4 (snp is 3 chars)
      cnames[ids] <- paste(prefix, substring(temp, 4), sep="")

      # Interactions
        if (max > 2) {
          if (facFlag) {
            string <- "_"
          } else {
            string <- ""
          }

          ids  <- assign == 3
          temp <- cnames[ids]
          nint <- length(intNames)

          # Remove the string "snp"
          temp <- substring(temp, 4)
          temp <- unlist(strsplit(temp, ":", fixed=TRUE))
          n    <- length(temp)
  
          # temp is a character vector, every 2 consecutive elements
          # was from the same list  
          even <- seq(from=2, to=n, by=2)
          temp[even] <- substring(temp[even], 6)
          odd <- seq(from=1, to=n-1, by=2)
          temp[odd] <- paste(prefix, temp[odd], sep="")        
          if (nint == 1) {
            temp[even] <- intNames
          } 
          cnames[ids] <- paste(temp[odd], string, temp[even], sep="")

        } # END: if (max > 2)

    } # END: if (max > 1)

    colnames(mat) <- cnames

  } # END: if (pFlag)

  # Remove columns if needed
  if (inc.snp) {
    if (facFlag) {

      # Make sure baseline column is not in the model
      if (pFlag) {
        var <- paste(prefix, snp.base, sep="")
        if (intFlag) {
          var <- c(var, paste(prefix, snp.base, "_", intNames, sep=""))
        }
      } else {
        var <- paste("snp", snp.base, sep="")
        if (intFlag) {
          var <- c(var, paste(var, ":X.int", intNames, sep=""))
        }
      }
      temp <- var %in% colnames(mat)
      if (any(temp)) mat <- removeOrKeepCols(mat, var[temp], which=-1)
    }
  } else {
    # Get the snp column numbers
    ids      <- attributes(mat)$assign == 2
    if (any(ids)) {
      snp.cols <- (1:length(ids))[ids]
      mat <- removeOrKeepCols(mat, snp.cols, which=-1)
    } 
  }

  mat

} # END: getDesignMatrix

# Function to compute an effects table and standard errors
effects.init <- function(parms, cov, var1, var2, levels1, levels2, 
                     base1=0, base2=0, int.var=NULL, effects=1,
                     sep1="_", base2.name="baseline", base1.name="baseline") {

  # parms       Vector of parameter estimates
  # var1        Vector of parameter names for one of the variables. If the
  #             length of this vector is > 1, then it is assumed that var1
  #             is a categorical variable.
  # var2        
  # levels1     Levels for var1 (snp)
  # levels2     Levels for var2 Can be NULL for a categorical variable
  # base1
  # base2
  # int.var     For no interaction, set to NULL. Otherwise a var1 x var2 
  #             matrix of interaction parameter names. int.var can also
  #             be a vector if length(var1) = 1 or length(var2) = 1.
  #             The order of this matrix must match the order of var1 and
  #             var2.
  #             The default is NULL.
  # effects     1 or 2  1 = joint, 2 = stratified
  # sep1        String to separate var1 with its levels
  #             The default is "_".

  int.flag <- !is.null(int.var)
  nv1      <- length(var1)
  nv2      <- length(var2)
  contv1   <- nv1 == 1  
  contv2   <- nv2 == 1
  joint    <- effects == 1

  # Check the number of interaction variables
  if (int.flag) {
    if (length(int.var) != nv1*nv2) {
      stop("ERROR with int.var")
    }
  }

  # Get the levels of the continuous variables, otherwise we assume
  #  the values are 0-1. 
  if (contv1) {
    nlev1   <- length(levels1)
    p1      <- rep(parms[var1], times=nlev1)
    names1  <- paste(var1, levels1, sep=sep1)
    v1      <- rep(var1, times=nlev1)
  } else {
    levels1 <- c(0, rep.int(1, times=nv1))
    nlev1   <- length(levels1)
    p1      <- c(0, parms[var1])
    base1   <- 0
    names1  <- c(base1.name, var1)
    v1      <- c(var1[1], var1)
  }
  if (contv2) {
    nlev2   <- length(levels2)
    p2      <- rep(parms[var2], times=nlev2)
    names2  <- paste(var2, levels2, sep="_")
    v2      <- rep(var2, times=nlev2)
  } else {
    levels2 <- c(0, rep.int(1, times=nv2))
    nlev2   <- length(levels2)
    p2      <- c(0, parms[var2])
    base2   <- 0
    names2  <- c(base2.name, var2) 
    v2      <- c(var2[1], var2)
  }

  # Initialize the matrix of interaction parms
  pint <- matrix(data=0, nrow=nlev1, ncol=nlev2)

  # Get the interaction parm
  if (int.flag) {
    # Remove dimension of int.var if <= 1 categorical var
    if (sum(contv1 + contv2) != 0) dim(int.var) <- NULL 

    if (contv1 && contv2) {
      pint[,] <- parms[int.var]
      int.var <- matrix(data=int.var, nrow=nlev1, ncol=nlev2)
    } else {
      if (!contv1 && contv2) {
        temp <- c(0, parms[int.var])
        for (i in 1:nlev2) pint[, i] <- temp
        int.var <- rep(c(int.var[1], int.var), times=nlev2)
        dim(int.var) <- c(nlev1, nlev2)
      } else if (contv1 && !contv2) {
        temp <- c(0, parms[int.var])
        for (i in 1:nlev1) pint[i, ] <- temp
        int.var <- rep(c(int.var[1], int.var), each=nlev1)
        dim(int.var) <- c(nlev1, nlev2)
      } else {
        # Both are categorical
        pint[1, ] <- 0
        for (i in 2:nlev1) {
          pint[i, ] <- c(0, parms[int.var[i-1,]])
        }
       
        # Add baselines to int.var
        temp <- int.var
        int.var <- matrix(data="", nrow=nlev1, ncol=nlev2)
        int.var[2:nlev1, 2:nlev2] <- temp

        # We can assign anything to these 
        int.var[, 1] <- temp[1,1]
        int.var[1, ] <- temp[1,1]
      }
    }
  } 

  # Get the baseline value(s)
  if (joint) {
    base <- base1*p1[1] + base2*p2[1] + base1*base2*pint[1, 1]
    base <- exp(base)
    base <- rep(base, times=nlev2)
  } else {
    base <- rep(NA, times=nlev2)
    for (j in 1:nlev2) {
      temp <- base1*p1[1] + levels2[j]*p2[j] + base1*levels2[j]*pint[1, j]
      base[j] <- exp(temp) 
    }
  }

  # Initialize matrix of effects 
  eff <- matrix(data=NA, nrow=nlev1, ncol=nlev2)
  colnames(eff) <- names2
  rownames(eff) <- names1
  # Loop over the levels
  for (i in 1:nlev1) {
    for (j in 1:nlev2) {
      temp <- levels1[i]*p1[i] + levels2[j]*p2[j] + 
              levels1[i]*levels2[j]*pint[i, j]
      temp <- exp(temp)
      #if ((joint) || (levels1[i] != base1)) temp <- temp/base[j]    
      #eff[i, j] <- temp
      eff[i, j] <- temp/base[j]
    }
  }

  # Initialize matrix for standard errors
  se <- matrix(data=NA, nrow=nlev1, ncol=nlev2)
  colnames(se) <- names2
  rownames(se) <- names1

  # Get base levels
  sebase1 <- rep(base1, times=nlev1)
  if (joint) {
    sebase2 <- rep(base2, times=nlev2)
  } else {
    sebase2 <- levels2
  }

  # Loop over the levels
  for (i in 1:nlev1) {
    for (j in 1:nlev2) {
      b1 <- sebase1[i]
      b2 <- sebase2[j]

      # Get vector of variables and coefficients
      vars <- c(v1[i], v2[j])
      coef <- c(levels1[i] - b1, levels2[j] - b2)

      if (int.flag) {
        vars <- c(vars, int.var[i, j])
        coef <- c(coef, levels1[i]*levels2[j] - b1*b2)
      }

      # SE
      se[i, j] <- sqrt(getLinearComb.var(vars, cov, coef=coef))
    }
  }
 
  logEffects <- log(eff)
  lower95    <- exp(logEffects - 1.96*se)
  upper95    <- exp(logEffects + 1.96*se)  

  list(effects=eff, lower95=lower95, upper95=upper95, 
       logEffects=logEffects, logEffects.se=se)

} # END: effects.init

# Function to compute an effects table and standard errors from
#  the snp.logistic output
snp.effects <- function(fit, var, var.levels=c(0,1), method=NULL) {

  # fit      Output from snp.logistic or snp.matched
  # var      Name of variable to get effects with snp
  # var.levels   (For continuous var) First level is assumed to be
  #          the baseline level 
  #          The default is NULL.

  if (length(var) != 1) stop("Only 1 variable can be specified")

  # Determine the input object
  temp <- class(fit)
  if (temp == "snp.logistic") {
    which   <- 1
    snp     <- fit$model.info$snpName
    methods <- c("UML", "CML", "EB")
    cnames  <- colnames(fit$UML$cov)
    if (is.null(fit$UML)) return(NULL)
  } else if (temp == "snp.matched") {
    which   <- 2
    snp     <- fit$model.info$snp.vars
    methods <- c("CLR", "CCL", "HCL")
    for (m in methods) {
      temp <- fit[[m, exact=TRUE]]
      if (!is.null(temp)) {
        cnames <- colnames(temp$cov)
        break
      }
    }
  } else {
    stop("fit must be of class snp.logistic or snp.matched")
  }

  if (!is.null(method)) {
    temp <- methods %in% method
    methods <- methods[temp]
    if (!length(methods)) stop("Incorrect method")
  }

  levels <- var.levels
  sep1   <- "_"
  nsnp   <- length(snp)

  if (!(var %in% fit$model.info$main.vars)) {
    stop("var must be a main effect variable")
  }
  
  if (var %in% fit$model.info$factors) {
    
    facFlag    <- 1
    levels     <- levels(fit$model.info$data[, var])
    temp2      <- paste(var, "_", levels, sep="") 
  
    # Get the variables in the model fit
    temp       <- temp2 %in% cnames
    var2       <- temp2[temp]
    base2.name <- temp2[!temp] 
    temp       <- length(base2.name)
    if (!temp) base2.name <- "baseline"
    if (temp > 1) base2.name <- base2.name[1]
  } else {
    facFlag    <- 0
    var2       <- var
    base2.name <- NULL
  }

  levels1 <- 0:2
  if (is.null(levels)) {
    if (is.factor(fit$model.info$data[, var])) {
      base2 <- 0
    } else {
      stop("levels must be specified for a continuous var")
    }
  } else {
    base2 <- levels[1]
  }
  
  # Initialize the return list
  ret   <- list()

  intFlag <- 0
  int.var <- NULL
  if (var %in% fit$model.info$int.vars) intFlag <- 1

  for (var1 in snp) {
    if (intFlag) int.var <- paste(var1, ":", var2, sep="")
  
    for (method in methods) {
      temp2 <- fit[[method]]
      if (is.null(temp2)) next
      tlist <- list()
      eff1 <- effects.init(temp2$parms, temp2$cov, var1, var2, 
                 levels1, levels, base1=0, base2=base2,
                 int.var=int.var, effects=1, sep1=sep1,
                 base2.name=base2.name)
      eff2 <- effects.init(temp2$parms, temp2$cov, var1, var2, 
                 levels1, levels, base1=0, base2=base2,
                 int.var=int.var, effects=2, sep1=sep1,
                 base2.name=base2.name)
      eff3 <- effects.init(temp2$parms, temp2$cov, var2, var1, 
                 levels, levels1, base1=base2, base2=0,
                 int.var=int.var, effects=2, sep1=sep1,
                 base2.name="0")

      # Set attributes  
      attr(eff1, "var1")    <- var1 
      attr(eff1, "var2")    <- var2 
      attr(eff1, "levels1") <- levels1 
      attr(eff1, "levels2") <- levels 
      attr(eff2, "var1")    <- var1
      attr(eff2, "var2")    <- var2
      attr(eff2, "levels1") <- levels1
      attr(eff2, "levels2") <- levels
      attr(eff3, "var1")    <- var2
      attr(eff3, "var2")    <- var1
      attr(eff3, "levels1") <- levels
      attr(eff3, "levels2") <- levels1

      temp <- list(JointEffects=eff1, StratEffects=eff2, StratEffects.2=eff3)
      #temp <- list(JointEffects=eff1, StratEffects=eff2)

      class(temp) <- "snp.effects.method" 

      if (nsnp == 1) {
        ret[[method]] <- temp
      } else {
        tlist[[method]] <- temp
      }      
    }
    if (nsnp != 1) ret[[var1]] <- tlist
  }
  class(ret) <- "snp.effects"

  ret

} # END: snp.effects

# Function to compute the variance of a linear combination of parms
getLinearComb.var <- function(vars, cov, coef=NULL) {

  # vars       Vector variable names or indices in cov
  # cov        Covariance matrix 
  # coef       Coefficients for vars. The order must be the same
  #            The default is all coefficients are 1
 
  n <- length(vars)
  if (is.null(coef)) coef <- rep(1, times=n)
  if (n != length(coef)) stop("ERROR with parms and/or coef")

  sum <- 0
  # Sum up the variances
  for (i in 1:n) {
    sum <- sum + coef[i]*coef[i]*cov[vars[i], vars[i]]
  }
  if (n == 1) return(sum)

  # Sum up the covariances
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      sum <- sum + 2*coef[i]*coef[j]*cov[vars[i], vars[j]]
    }
  }

  sum

} # END: getLinearComb.var

# Function to compute odds ratios and standard errors from
# log-odds ratios and se. 
getORfromLOR <- function(data, lor.var, lor.se.var,
                         or="OR", or.se="OR.SE") {

  data[, or] <- exp(data[, lor.var])

  # Change the standard errors
  temp <- data[, lor.se.var]
  data[, or.se] <- temp*data[, or]

  data

} # END: getORfromLOR

# Function to compute confidence intervals
getCI <- function(data, var, var.se, se=1, 
                  lower="LOWER", upper="UPPER") {
 
  zcrit <- qnorm(0.05/2, lower.tail=FALSE)
  if (se) {
    temp  <- zcrit*data[, var.se]
  } else {
    temp  <- zcrit*sqrt(data[, var.se])
  }

  data[, lower] <- data[, var] - temp
  data[, upper] <- data[, var] + temp 

  data

} # END: getCI

# Function to compute normal pvalues
pvalue.normal <- function(test, sided=2) {
  
  sided*pnorm(test, lower.tail=FALSE)

} # END: pvalue.normal

# Function to perform an unadjusted analysis based on genotype
#  frequency counts for cases and controls
unadjustedGLM.counts <- function(file.list, op=NULL) {

  #################################################################
  # file.list      List of type file.list with additional fields:
  #  caseCounts    Variables for the genotype frequencies 0, 1, 2
  #                among the cases.
  #                No default
  #  controlCounts Variables for the genotype frequencies 0, 1, 2
  #                among the controls.
  #                No default.
  #  covar         Covariate in the model.
  #                The default is c(0, 1, 2)
  #  caseOrder     Order for caseCounts in terms of covar
  #                The default is c(1, 2, 3)
  #  controlOrder  Order for controlCounts in terms of covar
  #                The default is c(1, 2, 3)
  #  caseSep       The default is "/"
  #  controlSep    The default is "/"
  #################################################################
  # op            List with names:
  #  outfile      The default is NULL
  #  copyVars     Variables to copy to the output data set
  #################################################################

  file.list <- default.list(file.list,
       c("file", "file.type", "header", "delimiter", "caseCounts",
         "controlCounts", "covar", "caseOrder", "controlOrder",
         "caseSep", "controlSep"),
       list("ERROR", 3, 1, "\t", "ERROR", "ERROR", c(0,1,2),
            1:3, 1:3, "/", "/"),
       error=c(1,0,0,0,1,1,0,0,0,0,0)
       )
  
  covar  <- file.list$covar
  nn     <- length(covar)
  v1     <- file.list$caseCounts
  v0     <- file.list$controlCounts
  n1     <- length(v1)
  n0     <- length(v0)
  order1 <- file.list$caseOrder
  order0 <- file.list$controlOrder
  if (n1 == nn) {
    v1    <- v1[order1]
    flag1 <- 1
  } else {
    flag1 <- 0
    sep1  <- file.list$caseSep
  }
  if (n0 == nn) {
    v0    <- v0[order0]
    flag0 <- 1
  } else {
    flag0 <- 0
    sep0  <- file.list$controlSep
  }

  # Read in the data
  x <- loadData(file.list$file, file.list)

  vars <- getListName(op, "copyVars")
  if (!is.null(vars)) {
    x <- removeOrKeepCols(x, c(vars, v1, v0), which=1)
  } else {
    x <- removeOrKeepCols(x, c(v1, v0), which=1)
  }
  x  <- unfactor.all(x)
  nr <- nrow(x)

  # Add new variables
  newVars <- c("beta", "se", "test", "pvalue")
  for (var in newVars) x[, var] <- NA

  mat <- matrix(data=NA, nrow=nn, ncol=2)
  for (i in 1:nr) {
    mat[] <- NA

    # Case counts
    if (flag1) {
      mat[, 1] <- as.numeric(unlist(x[i, v1]))
    } else {
      mat[, 1] <- as.numeric(getVecFromStr(x[i, v1], delimiter=sep1))
    }
    
    # Control counts
    if (flag0) {
      mat[, 2] <- as.numeric(unlist(x[i, v0]))
    } else {
      mat[, 2] <- as.numeric(getVecFromStr(x[i, v0], delimiter=sep0))
    }
    temp <- try(glm(mat~covar, family=binomial), silent=TRUE)
    if (class(temp)[1] != "try-error") {
      temp <- summary(temp)$coefficients
      if (nrow(temp) == 2) x[i, newVars] <- temp[2,]
    } 
  }

  temp <- getListName(op, "outfile")
  if (!is.null(temp)) {
    write.table(x, file=temp, sep="\t", row.names=FALSE, quote=FALSE)
  }

  x

} # END: unadjustedGLM.counts

# Function to compute the inflation factor
inflationFactor <- function(tests, squared=0, df=1) {

  # tests    Vector of Z-test statistics or squared Z-test
  #          statistics, or p-values
  # squared  0 or 1 if the tests are already squared

  i1  <- qchisq(0.5, df=df)
  if (!squared) tests <- tests*tests
  i2  <- median(tests, na.rm=TRUE)
  ret <- i2/i1
  ret

} # END: inflationFactor

# Function to compute genomic control adjusted p-values
GC.adj.pvalues <- function(tests, pvals=NULL, ifac=NULL) {

  # tests    Vector of Z-test statistics
  # ifac     NULL or the inflation factor to use
  #          The default is NULL

  test2 <- tests*tests
  if (is.null(pvals)) pvals <- 2*pnorm(abs(tests), lower.tail=FALSE)
  if (is.null(ifac)) ifac <- inflationFactor(test2, squared=1)
  ret  <- pchisq(test2/ifac, df=1, lower.tail=FALSE)
  ret

} # END: GC.adj.pvalues

# Function to swap 2 columns of a covariance matrix
swap2cols.cov <- function(mat, col1, col2, errorCheck=0) {

  # mat         Covariance matrix
  # col1        Column name or number
  # col2        Column name or number
  # errorCheck  0 or 1. If set to 1, then the matrix mat must
  #             have both row and column names for the error
  #             check to be done.
  #             The default is 0.
 
  if (col1 == col2) return(mat)

  nr <- nrow(mat)
  nc <- ncol(mat)
  if (nr != nc) stop("ERROR in swap2cols.cov: mat is not a square matrix")

  # Get row and column names
  cnames <- colnames(mat)
  rnames <- rownames(mat)
  cflag  <- !is.null(cnames)
  rflag  <- !is.null(rnames)

  # Get column numbers if variable names are passed in
  cflag1 <- is.character(col1)
  cflag2 <- is.character(col2)
  if (cflag1 || cflag2) {
    if (!cflag) stop("ERROR in swap2col.cov: mat must have column names")
    col1 <- (1:nr)[temp == col1]
    col2 <- (1:nr)[temp == col2]
  }

  # Let col1 be the smaller column
  temp <- col1
  if (col2 < col1) {
    col1 <- col2
    col2 <- temp
  }

  if (errorCheck && cflag && rflag) {
    mat0 <- mat
  } else {
    errorCheck <- 0
  }

  # Save col1 
  save <- mat[, col1]

  # Change columns col1 and col2
  if (col1 > 1) {
    vec <- 1:(col1-1) 
    mat[vec, col1] <- mat[vec, col2]
    mat[vec, col2] <- save[vec]
  }
  mat[col1, col1] <- mat[col2, col2]
  if (col1+1 < col2) {
    vec <- (col1+1):(col2-1)
    mat[vec, col1] <- mat[vec, col2]
    mat[vec, col2] <- save[vec]
  }
  mat[col2, col1] <- mat[col1, col2]
  mat[col2, col2] <- save[col1]
  if (col2 < nr) {
    vec <- (col2+1):nr
    mat[vec, col1] <- mat[vec, col2]
    mat[vec, col2] <- save[vec]
  }

  # Change rows col1 and col2
  vec  <- c(col1, col2)
  temp <- (1:nr)[-vec] 
  for (i in vec) {
    for (j in temp) mat[i, j] <- mat[j, i]
  }

  # Change row/col names
  if (cflag) {
    temp          <- cnames[col1]
    cnames[col1]  <- cnames[col2]
    cnames[col2]  <- temp
    colnames(mat) <- cnames
  }
  if (rflag) {
    temp          <- rnames[col1]
    rnames[col1]  <- rnames[col2]
    rnames[col2]  <- temp
    rownames(mat) <- rnames
  }

  # Error check
  if (errorCheck) {
    for (i in rnames) {
      for (j in cnames) {
        if (mat[i, j] != mat0[i, j]) {
          stop("ERROR in swap2cols: with the error check")
        }
      }
    }
  } 

  mat

} # END: swap2cols.cov

# Function to perform test for heterogeneity (logistic regression only)
heterTest <- function(data, X.vars, group.var, snp.var, op=NULL) {

  # data       Data frame
  # X.vars     Variables to be adjusted for
  # group.var  Variable that defines the groups 
  # snp.var    Name of the SNP variable
  # op         List with names:
  #   print    0 or 1 to print model summaries
  #            The default is 0
  #   levels   The levels of group.var that are to be used.
  #            The default is NULL so that all levels of group.var
  #            will be used.

  op <- default.list(op, c("print"), list(0))
  print <- op$print

  levels <- getListName(op, "levels")
  if (is.null(levels)) levels <- unique(data[, group.var])
  nlevels <- length(levels)
  
  # Variables in the model
  X.vars <- unique(c(X.vars, snp.var))

  nTest <- 0
  minP  <- 9999
  for (i in 1:(nlevels - 1)) {
    for (j in (i+1):nlevels) {
      # Get the subgroups to be included in the model
      vec  <- c(levels[i], levels[j])
      temp <- data[, group.var] %in% vec
      temp[is.na(temp)] <- FALSE   
      dat2 <- data[temp, ]

      # Define the design matrix
      X <- dat2[, X.vars]
      
      # Define the response. Set one group to 1
      y <- rep.int(0, times=nrow(X))
      temp <- dat2[, group.var] == vec[1]
      y[temp] <- 1

      nTest <- nTest + 1
      fit <- try(glm(y ~ ., data=X, family="binomial"), silent=TRUE)
      if ("try-error" %in% class(fit)) next
      if (!fit$converged) next

      s <- summary(fit)
      if (print) {
        print(vec)
        print(table(dat2[, group.var]))
        print(s)
      }
      temp <- s$coefficients
      pval <- temp[snp.var, 4]
      minP <- min(minP, pval)    
    }
  }  

  minP <- min(1, minP*nTest)

  minP

} # heterTest

# Function to compute frequency counts for a vector of intervals.
# The returned object will be a data frame of frequency counts for
# the partitioned intervals or if data is passed in, then the returned
# object will be data with a new column.
freqCounts.var <- function(vec, intervals, leftEndClosed=1, data=NULL, 
                           newVar="newVar", newVarCats=NULL) {
 
  # vec              Numeric vector
  #                  No default
  # intervals        Numeric vector of intervals. Ex: c(0, 10, 20, 50, 200)
  #                  No default
  # leftEndClosed    0 or 1 Set to 1 for intervals [a , b).
  #                  Note: if a=b, then the interval is [a, a] and the next
  #                  interval will be (a, b) 
  # data             Data frame to be returned with the new variable newVar
  #                  on it with the disjoint categories.
  #                  The default is NULL
  # newVar           The name of the new variable if data is passed in
  #                  The default is "newVar".  
  # newVarCats       Vector of categories for newVar
  #                  The default is NULL

  intervals <- sort(intervals)
  intervals <- c(-Inf, intervals, Inf)
  n         <- length(intervals) - 1
  ret <- data.frame(rep.int(0, times=n+1))
  colnames(ret) <- "FREQ"
  if (leftEndClosed) {
    left  <- "["
    right <- ")"
    lop   <- ">="
    rop   <- "<"
  } else {
    left  <- "("
    right <- "]"
    lop   <- ">"
    rop   <- "<="
  }

  dFlag <- !is.null(data)
  if (dFlag) {
    data[, newVar] <- "MISSING"
    # Check for integers
    temp <- (vec == as.integer(vec))
    temp[is.na(temp)] <- TRUE
    if (all(temp)) {
      intFlag <- 1
    } else {
      intFlag <- 0
    }
    catFlag <- !is.null(newVarCats)
  } 

  rnames <- NULL
  flag   <- 0
  for (i in 1:n) {
    a <- intervals[i]
    b <- intervals[i+1]
    
    if (a == b) {
      text  <- paste("(vec == ", a, ")", sep="") 
      rtemp <- paste("[", a, ", ", b, "]", sep="")
      dstr  <- as.character(a)
      flag <- 1
    } else {
      if (flag) {
        # Previous had a == b, so left interval should be open to obtain disjoint sets
        text  <- paste("(vec", ">", a, ") & (vec", rop, b, ")", sep="")
        rtemp <- paste("(", a, ", ", b, right, sep="")
      } else {
        text  <- paste("(vec", lop, a, ") & (vec", rop, b, ")", sep="")
        rtemp <- paste(left, a, ", ", b, right, sep="")
      } 
      flag <- 0
      if (dFlag) {
        if (a == -Inf) {
          if (leftEndClosed) {
            dstr <- paste("lt", b, sep="")
          } else {
            dstr <- paste("lteq", b, sep="")
          }
        } else if (b == Inf) {
          if ((intFlag) & (!leftEndClosed)) a <- a + 1
          dstr <- paste(a, "plus", sep="")
        } else {
          if (intFlag) {
            if (!leftEndClosed) {
              a <- a + 1
            } else {
              b <- b - 1
            }
          }
          dstr <- paste(a, "to", b, sep="")
        }
      }
    }
    # Get the logical vector
    temp <- eval(parse(text=text))
    temp[is.na(temp)] <- FALSE
    ret[i, "FREQ"] <- sum(temp, na.rm=TRUE)
    rnames <- c(rnames, rtemp)

    if ( (dFlag) & (any(temp)) ){ 
      if (catFlag) dstr <- newVarCats[i]
      data[temp, newVar] <- dstr
    }
  }
  
  # Count the number of missing
  rnames <- c(rnames, "NA")
  ret[n+1, "FREQ"] <- sum(is.na(vec))
  rownames(ret) <- rnames

  if (dFlag) {
    print(ret)
    return(data)
  }
  ret

} # freqCounts.var

# Function to standardize a continuous vector
standardize.z <- function(vec) {

  mu <- mean(vec, na.rm=TRUE)
  se <- sqrt(var(vec, na.rm=TRUE))
  ret <- (vec - mu)/se
  ret 

} # END: standardize.z

# Function to create a design matrix
dsgnMat <- function(data, vars, facVars, removeInt=1) {

  # data        Data frame
  # vars        Character vector of variable names or a formula
  # facVars     Character vector of factor names

  if (is.null(vars)) return(list(designMatrix=NULL, newVars=NULL)) 

  # See if vars is a character string containing a formula
  if ((length(vars) == 1) && (substr(vars, 1, 1) == "~")) {
    vars <- as.formula(vars)
  }

  # Determine if vars is a formula
  if ("formula" %in% class(vars)) {
    # Get the design matrix
    design <- model.matrix(vars, data=data)

    # Remove the intercept, if needed
    newVars <- colnames(design)
    if (removeInt) {
      if (newVars[1] == "(Intercept)") {
        design  <- removeOrKeepCols(design, 1, which=-1)
        newVars <- newVars[-1]
      }
    }

    return(list(designMatrix=design, newVars=newVars))    
  }

  design  <- removeOrKeepCols(data, vars, which=1)
  newVars <- NULL
  if (!is.null(facVars)) {
    temp <- vars %in% facVars
    if (any(temp)) {
      temp    <- vars[temp]
      temp    <- createDummy(design, vars=temp)
      design  <- temp$data
      newVars <- temp$newVars
    }
  } 
  design <- as.matrix(design)

  # Check for constant variables
  design <- checkForConstantVar(design, msg=1)$data

  if (!removeInt) {
    # Add intercept
    cnames <- colnames(design)
    design <- cbind(1, design)
    colnames(design) <- c("Intercept", cnames)
  }

  # Make sure matrix is numeric
  d <- dim(design)
  cnames <- colnames(design)
  design <- as.numeric(design)
  dim(design) <- d
  colnames(design) <- cnames

  list(designMatrix=design, newVars=newVars)

} # END: dsgnMat

# Function to return genotype stats
getGenoStats <- function(vec, MAF=1, freqCounts=1) {

  # Vec must be numeric coded as 0-1-2 or NA

  if (MAF) {
    MAF2 <- getMAF(vec)
  } else {
    MAF2 <- NULL
  }
  if (freqCounts) {
    counts <- getGenoCounts(vec)
  } else {
    counts <- NULL
  }
  n.miss   <- sum(is.na(vec))
  len      <- length(vec)
  missRate <- n.miss/len
  n        <- len - n.miss

  list(MAF=MAF2, n.miss=n.miss, freqCounts=counts, missRate=missRate, n=n)

} # END: getGenoStats

# Function to return the extreme subjects based on a score
getExtremeSubs <- function(id, cc, score, n=500, study=NULL) {

  # Get cases with lowest score
  temp   <- cc == 1
  id1    <- id[temp]
  temp   <- sort(score[temp], index.return=TRUE)$ix
  id1    <- id1[temp]
  case   <- id1[1:n]

  # Get controls with highest score
  temp   <- cc == 0
  id1    <- id[temp]
  temp   <- sort(score[temp], decreasing=TRUE, index.return=TRUE)$ix
  id1    <- id1[temp]
  cntl   <- id1[1:n]

  if (!is.null(study)) {
    # Save info for the chosen subjects
    temp <- id %in% c(case, cntl)
    s2   <- study[temp]
    cc2  <- cc[temp]

    # Let each study have an equal number of cases and controls
    ustudy <- unique(study)
    nstudy <- length(ustudy)

    # Remove chosen subjects
    temp  <- !(id %in% c(case, cntl))
    id    <- id[temp]
    cc    <- cc[temp]
    score <- score[temp]
    study <- study[temp]   
    for (i in 1:nstudy) {
      # Get the study counts for the chosen subjects
      ncase <- sum((s2 == ustudy[i]) & (cc2 == 1))
      ncntl <- sum((s2 == ustudy[i]) & (cc2 == 0))

      if (ncase < ncntl) {
        value <- 1
        dec   <- FALSE
      } else if (ncase > ncntl) {
        value <- 0
        dec   <- TRUE
      } else {
        next
      }
      # The number of subjects we need
      m <- abs(ncase - ncntl)

      temp   <- (cc == value) & (study == ustudy[i])
      id1    <- id[temp]
      temp   <- sort(score[temp], decreasing=dec, index.return=TRUE)$ix
      id1    <- id1[temp]
      subs   <- id1[1:m]

      if (value == 1) {
        case <- c(case, subs)
      } else {
        cntl <- c(cntl, subs)
      }
    }
  }

  list(case=case, control=cntl) 

} # END: getExtremeSubs

# Function to add the OR and confidence intervel onto a data frame or matrix
getOR.CI <- function(x, op=NULL) {

  op <- default.list(op, 
         c("beta.var", "se.var", "OR.name", "CI.name", "alpha", "digits"),
           list("Beta", "SE", "OR", "OR.CI", 0.05, 4))

  z <- qnorm(1 - (op$alpha)/2)

  cn   <- colnames(x)
  beta <- as.numeric(x[, op$beta.var])
  se   <- as.numeric(x[, op$se.var])
  or   <- exp(beta)
  l    <- exp(beta - z*se)
  l    <- round(l, digits=op$digits)
  u    <- exp(beta + z*se)
  u    <- round(u, digits=op$digits)
  ci   <- paste("(", l, ", ", u, ")", sep="") 
  x    <- cbind(x, or, ci) 
  colnames(x) <- c(cn, op$OR.name, op$CI.name)

  x

} # END: getOR.CI

# Function to generete multi-variate random normal vectors
myrmvnorm <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
    method = c("eigen", "svd", "chol"), saveObj=NULL) {

  # Taken from rmvnorm function in the mvtnorm package. It has been modified to
  #   return and input an object that is used for generating the random vectors, so 
  #   that a singular value decomposition or cholesky decomposition does not
  #   have to be redone.
  # n
  # mean
  # sigma
  # method
  # saveObj   For efficiency in calling the function more than once with 
  #           the same input data

    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
        warning("sigma is numerically not symmetric")
    }
    method <- match.arg(method)
    if (is.null(saveObj)) {
      if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
            t(ev$vectors)
      }
      else if (method == "svd") {
        sigsvd <- svd(sigma)
        if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
      }
      else if (method == "chol") {
        retval <- chol(sigma, pivot = TRUE)
        o <- order(attr(retval, "pivot"))
        retval <- retval[, o]
      }
      saveObj <- retval

    } # END: if (is.null(retval)) 
    
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% saveObj
    retval <- sweep(retval, 2, mean, "+")
    colnames(retval) <- names(mean)

    list(randomVectors=retval, saveObj=saveObj)

} # END: myrmvnorm





