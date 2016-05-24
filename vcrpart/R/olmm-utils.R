##' -------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2015-08-21
##'
##' Description:
##' Utility functions for the olmm function (see olmm.R). Some
##' functions are experimental and not listed in the namespace.
##'
##' References and dependencies:
##' statmod:         http://cran.r-project.org/web/packages/statmod/index.html
##' ucminf:          http://cran.r-project.org/web/packages/ucminf/index.html
##'
##' Contents:
##'
##' Functions for olmm:
##' olmm_expandQP:        expand grid for numerical integration
##' olmm_fn
##' olmm_gn
##' olmm_optim_setup:     set up algorithm function
##' olmm_optim_warnings:
##' olmm_coefShortLabs:   short labels for coefficient names
##' olmm_check_mm:        check and modify model matrix
##' olmm_start:           set initial values
##' olmm_mergeMm:         merge the predictor-variable and predictor-invariant
##'                       model matrices
##'                     
##' gcauchy:              derivate of dcauchy
##' glogis:               derivate of dlogis
##' gnorm:                derivate of dnorm
##'
##' Functions for estfun.olmm:
##' olmm_scoreVar:        computes the variance of the observation scores.
##' olmm_scoreCovWin:     computes the intra-subject covariance of
##'                       the observation scores.
##' olmm_scoreCovBet:     computes the between-subject covariance
##'                       of the observation scores.
##' olmm_f_decormat:      the equation to be optimized to zero:
##'                       the difference between the adjusted
##'                       intra-subject covariance and the adjusted
##'                       between-subject covariance of likelihood scores.
##' olmm_g_decormat:      the derivation of olmm_fDecorrelate.
##' olmm_rename:
##' olmm_decormat:        computes the transformation matrix for
##'                       removing the intra-subject
##'                       correlation of ML scores.
##'
##' Modifications:
##' 2014-09-08: partial substitution of 'rep' by 'rep.int'
##' 2014-09-07: updated descriptions for undocumented functions
##' 
##' To do:
##' - replace ranefChol fac initial value with covariance matrix
##' - add multiple family arguments
##' - make olmm_refit_MC working (compare version 0.1-5)
##' -------------------------------------------------------- #

##' -------------------------------------------------------- #
##' Expand quadrature points to the dimension of random
##' effects.
##'
##' @param x nodes (or weights) for one dimension.
##' @param q number of random coefficients.
##'
##' @return A grid matrix of nodes (or weights) of dimension q.
##' -------------------------------------------------------- #

olmm_expandQP <- function(x, q) {

  if (length(q)^q * q > 2^31 - 1)
    stop("number of quadrature weights is too large. Decrease 'nAGQ'")
  rval <- matrix(0, length(x)^q, q)
  for (i in 1:q) rval[, i] <- x[rep(1:length(x), each = length(x)^(i - 1L), length.out = length(x)^q)]
  return(rval)
}


olmm_fn <- function(par, restricted) {
  if (!exists("object")) object <- new(Class = "olmm")
  parNew <- object$coefficients
  parNew[!restricted] <- par[!restricted]
  .Call("olmm_update_marg", object, as.numeric(parNew), PACKAGE = "vcrpart")
  return(-object$logLik)
}


olmm_gr <- function(par, restricted) {
  if (!exists("object")) object <- structure(list(), class = "olmm")
  scoreNew <- object$score
  scoreNew[restricted] <- c(0.0)
  return(-scoreNew)
}


##' -------------------------------------------------------- #
##' Processes the algorithm argument.
##'
##' @param x       argument 'optim' of 'olmm' call
##' @param numGrad argument 'numGrad' of 'olmm' call
##' @param env     environment of the optimization 
##'
##' @return A prepared list to be assigned to the eval command
##'    invoking the optimization (using the eval function).
##' -------------------------------------------------------- #

olmm_optim_setup <- function(x, env = parent.frame()) {  

  numGrad <- x$numGrad
  if (x$fit == "optim" && is.null(x$method)) x$method <- "Nelder-Mead"
  if (x$fit == "optim" && !x$method %in% c("BFGS", "CG", "L-BFGS-B"))
    numGrad <- TRUE
  x <- x[!names(x) %in% c("start", "restricted")]
  rval <- list(par = NULL, # first 4 arguments must be exactly in that order!
               fn = olmm_fn,
               gr = if (numGrad) NULL else olmm_gr,
               restricted = NULL,
               fit = x$fit)
  rval <- append(rval, x[intersect(names(formals(rval$fit)), names(x))])
  if (rval$fit == "nlminb") names(rval)[1:3] <- names(formals(nlminb)[1:3])
  
  ## set environment for objective function and gradient
  environment(rval[[2]]) <- env
  if (!numGrad) environment(rval[[3]]) <- env
  
  return(rval)
}


olmm_optim_warnings <- function(output, FUN) {
  
  if (FUN == "optim") {
    switch(as.character(output$convergence),
           "1" = warning("Stopped by small step (xtol)."),
           "20" = warning("Inadmissible initial parameters."),
           "21" = warning("Inadmissible immediate parameters."),
           "10" = warning("Degeneracy of the Nelder-Mead simplex."),
           "51" = warning("Warning from 'L-BFGS-B'."),
           "52" = warning("Error from 'L-BFGS-B'."),
           NULL)
  }
  
  if (FUN == "nlminb") {
    if (output$convergence != 0) warning(output$message)
  }
  
  if (FUN == "ucminf") {
    
    switch(as.character(output$convergence),
           "2" =
           if (output$info["maxgradient"] > 1e-3)
           warning("Stopped by small step (xtol)."),
           "3" = 
           warning("Stopped by function evaluation limit (maxeval)."),
           "4" =
           if (output$info["maxgradient"] > 1e-3)
           warning("Stopped by zero step from line search."),
           "-2" =
           warning("Computation did not start: length(par) = 0."),
           "-4" =
           warning("Computation did not start: stepmax is too small."),
           "-5" =
           warning("Computation did not start: grtol or xtol <= 0"),
           "-6" =
           warning("Computation did not start: maxeval <= 0"),NULL)
  }
}


##' -------------------------------------------------------- #
##' Proposes abbreviations for names of coefficients.
##'
##' @param object an olmm object
##'
##' @return Vector of character strings with abbreviations
##'    for names of coefficients.
##' -------------------------------------------------------- #

olmm_coefShortLabs <- function(object) {

  ## abbreviations for predictor-variable fixed effects
  abbCe <- names(fixef(object, "predictor-variable"))
  abbCe <- gsub("Eta", "", abbCe, fixed = TRUE)
  abbCe <- gsub("(Intercept)", "I", abbCe, fixed = TRUE)
  abbCe <- abbreviate(abbCe)

  ## abbreviations for predictor-invariant fixed effects
  abbGe <- names(fixef(object, "predictor-invariant"))
  abbGe <- abbreviate(abbGe)

  ## abbreviations for Choleski factors
  q <- object$dims["q"]
  abbReCF <- paste("reCF", 1:(q * (q + 1L) / 2L), sep = "")
  
  return(c(abbCe, abbGe, abbReCF))
}


##' -------------------------------------------------------- #
##' Merge a model matrix of predictor variable effects with a
##' model matrix of predictor invariant effects
##'
##' @param x a model matrix for predictor-variable effects.
##' @param y a model matrix for predictor-invariant effects. 
##'
##' @return A model matrix.
##' -------------------------------------------------------- #

olmm_merge_mm <- function(x, y, deleteIntY = TRUE) {

  if (ncol(y) == 0L) deleteIntY <- FALSE
  if (deleteIntY & "(Intercept)" %in% colnames(y)) {
    rval <- cbind(x, y[, -1L, drop = FALSE])
    attr(rval, "assign") <- c(attr(x, "assign"), attr(y, "assign")[-1L])
    attr(rval, "merge") <- rep(c(1L, 2L), c(ncol(x), ncol(y) - 1L))
  } else {
    rval <- cbind(x, y)
    attr(rval, "assign") <- c(attr(x, "assign"), attr(y, "assign"))
    attr(rval, "merge") <- rep(c(1L, 2L), c(ncol(x), ncol(y)))
  }
  attr(rval, "contrasts") <-
    append(attr(x, "contrasts"), attr(y, "contrasts"))
  attr(rval, "contrasts") <-
    attr(rval, "contrasts")[!duplicated(names(attr(rval, "contrasts")))]
  rownames(rval) <- rownames(x)
  return(rval)
}


##' -------------------------------------------------------- #
##' Check and modify model matrix.
##'
##' @param x a model matrix.
##'
##' @return A model matrix.
##' -------------------------------------------------------- #

olmm_check_mm <- function(x) {
  
  qr.x <- qr(x, LAPACK = FALSE)
  rank <- qr.x$rank
  
  ## automated drops
  if (rank < ncol(x)) {
      
    warning("design matrix appears to be column rank-deficient ",
            "so dropping some coefs")
    
    dropterms <- function(x, keep) {
      rval <- x[, keep, drop = FALSE]
      attr(rval, "assign") <- attr(x, "assign")[keep]
      attr(rval, "merge") <- attr(x, "merge")[keep]
      attr(rval, "contrasts") <- attr(x, "contrasts")
      if ("orig.colnames" %in% names(attributes(x))) {
        attr(rval, "orig.colnames") <- attr(x, "orig.colnames")
      } else {
        attr(rval, "orig.colnames") <- colnames(x)
      }
      return(rval)
    }
   
    subs <- qr.x$pivot[1:qr.x$rank]
    x <- dropterms(x, subs)
    
    if (rank != qr(x, LAPACK = FALSE)$rank)
      stop("determination of full column rank design matrix failed")
  } else {

    attr(x, "orig.colnames") <- colnames(x)
  } 
  storage.mode(x) <- "double"
  return(x)
}


##' -------------------------------------------------------- #
##' Evaluate user-defined initial values for the parameters.
##'
##' @param start      named vector with user-defined initial
##'    values
##' @param dims       dimension slot of the olmm object
##' @param parNames   list with elements 'fixef' and
##'    'ranefCholFac' including the parameter names
##' @param X          model matrix for fixed effects
##' @param W          model matrix for random effects
##' @param eta        matrix for storing the linear predictor
##' @param ranefElMat matrix for transforming the ranefCholFac
##'    vectors to matrices
##'
##' @return A list with elements 'fixef', 'ranefCholFac' and
##'    'coefficients' that are used as initial values.
##' -------------------------------------------------------- #

olmm_start <- function(start, dims, parNames, X, W, eta, ranefElMat) {

  ## checks
  if (!is.null(start)) {
    if (!is.numeric(start) || is.null(names(start)))
      stop("'start' argument should be a named numeric vector")
    if (!all(names(start) %in% unlist(parNames)))
      start <- start[names(start) %in% unlist(parNames)]
      }
  
  ## set fixed effects
  
  ## default values
  intDef <- switch(as.character(dims["family"]),
                   "1" = qlogis(ppoints(dims["nEta"])),
                   rep.int(0.0, dims["nEta"]))
  fixef <- c(matrix(c(rep.int(intDef, dims["pInt"]), rep.int(0.0, dims["nEta"] * (dims["pCe"] - dims["pInt"]))), dims["pCe"], dims["nEta"], byrow = TRUE), rep.int(0.0, dims["pGe"]))
  names(fixef) <- parNames$fixef
  
  ## overwrite with 'start'
  subs <- intersect(names(fixef), names(start))
  fixef[subs] <- start[subs]
  
  ## make a fixed effects matrix
  fixef <- rbind(matrix(as.numeric(fixef[1:(dims["pCe"] * dims["nEta"])]), dims["pCe"], dims["nEta"], byrow = FALSE), if (dims["pGe"] > 0) matrix(rep(as.numeric(fixef[(dims["pCe"] * dims["nEta"] + 1):dims["p"]]), each = dims["nEta"]), dims["pGe"], dims["nEta"], byrow = TRUE) else NULL)
  rownames(fixef) <- colnames(X)
  colnames(fixef) <-  colnames(eta)
  
  if (dims["family"] == 3L && dims["pCe"] > 0L) {
    
    ## transform adjacent category parameters into baseline category
    ## parameters
    T <- 1 * lower.tri(diag(dims["nEta"]), diag = TRUE)
    fixef[1:dims["pCe"],] <- fixef[1:dims["pCe"],] %*% T
  }
  
  ## set random effect variance components
  
  ## default values
  ranefCholFac <- c(ranefElMat %*%  c(diag(dims["q"])))
  names(ranefCholFac) <- parNames$ranefCholFac
  
  ## overwrite with 'start'
  subs <- intersect(names(ranefCholFac), names(start))
  ranefCholFac[subs] <- start[subs]
  
  ## make a matrix
  ranefCholFac <-
    matrix(t(ranefElMat) %*% ranefCholFac, dims["q"], dims["q"])
  ## set row- and column names
  tmp <- c((if (dims["qCe"] > 0) paste("Eta", rep(seq(1, dims["nEta"], 1), each = dims["qCe"]), ":", rep.int(colnames(W)[attr(W, "merge") == 1], dims["nEta"]), sep = "") else NULL), (if (dims["qGe"] > 0) colnames(W)[attr(W, "merge") == 2] else NULL))
  rownames(ranefCholFac) <- tmp
  colnames(ranefCholFac) <- tmp
  
  if (dims["family"] == 3 && dims["qCe"] > 0) {
    
    ## transform adjacent category parameters into baseline category
    ## parameters
    for (i in 1:dims["qCe"]) {
      subs <- seq(i, dims["qCe"] * dims["nEta"], dims["qCe"])
      for (j in 1:(dims["nEta"] - 1)) {
        ranefCholFac[subs[j], ] <-
          colSums(ranefCholFac[subs[j:dims["nEta"]], ]) 
      }
    }
    ranefCholFac <- t(chol(ranefCholFac %*% t(ranefCholFac)))
  }
  
  ## coefficients
  coefficients <- numeric()
  if (dims["pCe"] > 0L) coefficients <- fixef[1:dims["pCe"],]
  if (dims["pGe"] > 0) coefficients <- c(coefficients,fixef[(dims["pCe"] + 1):dims["pEta"], 1])
  coefficients <- c(coefficients, c(ranefElMat %*% c(ranefCholFac)))
  names(coefficients) <- unlist(parNames)

  return(list(fixef = fixef,
              ranefCholFac = ranefCholFac,
              coefficients = coefficients))
}


##' -------------------------------------------------------- #
##' Computes variance of scores.
##'
##' @param scores  the scores of the model
##' @param subject the subject vector
##' -------------------------------------------------------- #

olmm_scoreVar <- function(scores, subject) {
  
  return(crossprod(scores) / nrow(scores))
}


##' -------------------------------------------------------- #
##' Computes within-subject covariance of scores.
##'
##' @param scores  a numeric matrix of scores.
##' @param subject a factor vector that assigns entries in
##'    'scores' to the grouping factor.
##' -------------------------------------------------------- #

olmm_scoreCovWin <- function(scores, subject) {
  
  Ni <- table(subject)
  return((crossprod(apply(scores, 2, tapply, subject, sum)) -
          crossprod(scores)) / sum(Ni * (Ni - 1)))
}


##' -------------------------------------------------------- #
##' Computes between-subject covariance of scores.
##'
##' @param scores  a numeric matrix of scores.
##' @param subject a factor vector that assigns entries in
##'    'scores' to the grouping factor.
##' -------------------------------------------------------- #

olmm_scoreCovBet <- function(scores, subject) {
  
  Ni <- table(subject)
  return(-crossprod(apply(scores, 2, tapply, subject, sum)) /
         (nrow(scores)^2 - nrow(scores) - sum(Ni * (Ni - 1))))
}


##' -------------------------------------------------------- #
##' Difference between adjusted within- and adjusted between-
##' subject covariance.
##' -------------------------------------------------------- #

olmm_f_decormat <- function(T, Tindex, sVar, sCovWin, sCovBet, Nmax) {

  adjScoreCovWin <-
    sVar %*% t(T) + T %*% sVar + (Nmax - 2) * T %*% sVar %*% t(T) +
    sCovWin + (Nmax - 2) * sCovWin %*% t(T) + (Nmax - 2) * T %*% t(sCovWin) +
      ((Nmax - 1)^2 - (Nmax - 2)) * T %*% sCovWin %*% t(T)

  adjScoreCovBet <-
    sCovBet +
        (Nmax - 1) * sCovBet %*% t(T) + (Nmax - 1) * T %*% t(sCovBet) +
        (Nmax - 1)^2 * T %*% sCovBet %*% t(T)
  
  rval <- adjScoreCovWin - adjScoreCovBet
  subs <- which(!duplicated(c(Tindex)) & c(Tindex) != 0)
  return(c(rval[subs]))
}

##' -------------------------------------------------------- #
##' Gradient for the olmm_f_decormat() function above
##' -------------------------------------------------------- #

olmm_g_decormat <- function(T, Tindex, sVar, sCovWin, sCovBet, Nmax) {

  k <- ncol(T)
  nPar <- max(Tindex)
  subs <- which(!duplicated(c(Tindex)) & c(Tindex) != 0)
  rval <- matrix(0, nPar, nPar)
  ## each iteration evaluates the gradient for one element in T
  ## regarding all elements in the objective function
  for (i in 1:nPar) {
    ind <- 1L * (Tindex == i)
    val <- T[which(ind == 1L)][1]

    gAdjScoreCovWin <-
      sVar %*% t(ind) + ind %*% sVar + 2 * val * (Nmax - 2) * ind %*% sVar %*% t(ind) +
        (Nmax - 2) * sCovWin %*% t(ind) + (Nmax - 2) * ind %*% sCovWin +
          2 * val * ((Nmax - 1)^2 - (Nmax - 2)) * ind %*% sCovWin %*% t(ind)

    gAdjScoreCovBet <-
      (Nmax - 1) * sCovBet %*% t(ind)  + (Nmax - 1) * ind %*% sCovBet +
           2 * val * (Nmax - 1)^2 * ind %*% sCovBet %*% t(ind)
    
    rval[, i] <- (gAdjScoreCovWin + gAdjScoreCovBet)[subs]
  }
  return(rval)
}


##' -------------------------------------------------------- #
##' Modifies the 'Eta[1-9]+' labels of coefficients of
##' 'olmm' objects.
##'
##' @param x       a character vector, a named numeric vector
##'    or a named matrix.
##' @param levels  the response levels.
##' @param family  an object of class 'family.olmm'
##' @param etalab  character string. Allowed are "int", "char" or
##'    "eta". 
##'
##' @return An object of the same class as 'x' but with new
##'    labels for category-specific coefficients.
##' -------------------------------------------------------- #

olmm_rename <- function(x, levels, family, etalab = c("int", "char", "eta")) {

  etalab <- match.arg(etalab)

  if(etalab != "eta") {

    if (etalab == "int") seq_along(levels)
    rename <- function(names) {
      names <- strsplit(names, ":") 
      names <- lapply(names, function(name) {
        if (substr(name[1L], 1, 3) == "Eta") {
          cat <- as.numeric(substr(name[1L], 4, nchar(name[1])))
          name[1L] <-
            switch(family$family,
                   cumulative = paste(levels[cat:(cat + 1L)], collapse = "|"),
                   adjacent = paste(levels[cat:(cat + 1L)], collapse = "|"),
                   baseline = paste(c(levels[cat], rev(levels)[1L]), collapse = "|"),
                   name[1L])
        }
        name <- paste(name, collapse = ":")
        return(name)
      })
      return(unlist(names))
    }

    if (is.character(x)) {
      x <- rename(x)
    } else if (is.matrix(x)) {
      if (!is.null(rownames(x))) {
        rownames(x) <- rename(rownames(x))
      } else {
        if (nrow(x) > 0L) rownames(x) <- rep.int("", nrow(x))
      }
      if (!is.null(colnames(x))) colnames(x) <- rename(colnames(x))
    } else if (is.numeric(x)) {
      if (!is.null(names(x))) names(x) <- rename(names(x))
    }
  }
  
  return(x)  
}


##' -------------------------------------------------------- #
##' Computes the T matrix for the pre-decorrelation transformation
##' of scores.
##'
##' @param scores  a matrix of scores.
##' @param subject a factor, the grouping factor for entries in
##'    'scores'.
##' @param parm    a subset of coefficients for which the matrix
##'    is to be computed.
##' @param control a list with parameters produced by predecor_control.
##'
##' @return A list with the frame and model matrices of the simulated
##'    predictors
##' -------------------------------------------------------- #

olmm_decormat <- function(scores, subject, control = predecor_control()) {
  
  stopifnot(inherits(control, "predecor_control"))
  Nmax <- max(table(subject))
  
  ## estimate variances and covariances
  sVar <- olmm_scoreVar(scores, subject)
  sCovWin <- olmm_scoreCovWin(scores, subject)
  sCovBet <- olmm_scoreCovBet(scores, subject)
  
  ## reduce to coefficient subset if intended
  k <- ncol(scores)
  
  ## set initial values
  T <- matrix(0, k, k, dimnames = list(rownames(sVar), colnames(sVar)))
  Tindex <- matrix(0, k, k, dimnames = list(rownames(sVar), colnames(sVar)))
  subs <- if (control$symmetric) lower.tri(T, TRUE) else matrix(TRUE, k, k)  
  if (control$minsize > 0L) # omit off-diagonal terms which are zero
    subs[crossprod(apply(abs(scores) > 0, 2, tapply, subject, sum)) <= control$minsize & crossprod(abs(scores) > 0) <= control$minsize] <- FALSE  
  nPar <- sum(subs)
  Tindex[subs] <- 1:nPar
  if (control$symmetric)
    Tindex[upper.tri(Tindex, FALSE)] <- t(Tindex)[upper.tri(Tindex, FALSE)]
  par <- rep.int(0, nPar)
  fEval <- rep.int(Inf, nPar)
  nit <- 0; error <- FALSE; eps <- 2 * control$reltol;
  
  ## optimize by Newton's algorithm
  while (!error && nit < control$maxit && eps >= control$reltol) {
    nit <- nit + 1
    fEval <- olmm_f_decormat(T, Tindex, sVar, sCovWin, sCovBet, Nmax)
    gEval <- olmm_g_decormat(T, Tindex, sVar, sCovWin, sCovBet, Nmax)
    par <- try(solve(gEval, -fEval) + par, silent = TRUE)
    eps <- max(abs(fEval / sCovBet[subs]))
    if (class(par) != "try-error") {
      T[] <- c(0, par)[Tindex + 1]
      if (control$verbose & nit > 1)
        cat("\nnit =", nit,
            "max|f| =", format(max(abs(fEval)), digits = 3, scientific = TRUE),
            "max|diff/f| =", format(eps, digits = 3, scientific = TRUE))
    }
    if (class(par) == "try-error" || eps > 1e100 || is.nan(eps))
      error <- TRUE
  } 
  
  ## check convergence
  if (nit >= control$maxit | error) {
    mess <- paste("optimization not converged: nit =", nit,
                  "reltol =", format(eps, digits = 3, scientific = TRUE), "\n")
    if (control$verbose) { cat("\n"); cat(mess); }
    if (!control$silent) warning(mess)
  } else {
    if (control$verbose)
      cat("\noptimization converged: nit =", nit,
          "max|diff/f| =", format(eps, digits = 3, scientific = TRUE), "\n")
  }
  
  attr(T, "conv") <- as.integer((nit != control$maxit) & !error)
  attr(T, "nit") <- nit
  attr(T, "eps") <- eps
  rownames(T) <- colnames(T) <- colnames(scores)
  
  ## return transformation matrix
  return(T)
}
