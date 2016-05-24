## =============================================================================
## Utility functions
## =============================================================================

## -----------------------------------------------------------------------------
## combinations of n elements from a vector p (length (p) > = n)
## -----------------------------------------------------------------------------
combin <- function(n, v) {
  if (n == 1)
    matrix(data = v, ncol = 1)
  else if (n >= length(v))
    matrix(data = v, nrow = 1)
  else
    rbind(cbind(v[1], combin(n-1, v[-1])), combin(n, v[-1]))
}

## -----------------------------------------------------------------------------
## Function to generate collinearities
## for a given set of parameter combinations (cc)
## -----------------------------------------------------------------------------
# new version dd 16/04/2014 (Jeremy David Silver)

collFun <- function (cc, normSens, npar, iNa) {
  n <- ncol(cc) ## calculate n once
  Collset <- matrix(0,ncol = npar + 2, nrow = nrow(cc)) ## allocate memory rather than use rbind
  for (i in 1:nrow(cc)) {
    ii <- cc[i, ]
    S <- normSens[, ii]
    Nident <- (iNa != 0 & iNa %in% ii)
    if (Nident) {
      id <- Inf
    } else {
## only calculate eigenvalues, use symmetry, make the cross product more efficient
      id <- 1/sqrt(max(0,tail(eigen(crossprod(S),only.values=TRUE,symmetric=TRUE)$value,1)))
    }
    psub <- rep(0, npar)
    psub[ii] <- 1
    Collset[i,] <- c(psub, n, id) ## insert into pre-allocated memory
  }
  return(Collset)
}

# Replaced 16/04/2014
collFun_OLD <- function(cc, normSens, npar, iNa) {
  Collset <- NULL

  for (i in 1:nrow(cc)) {  # for each parameter combination
    ii    <- cc[i,]          # vector with parameters in set
    S     <- normSens[,ii]   # select relevant columns of sensitivity functions

    Nident <- (iNa != 0 & iNa %in% ii)  # Is it nonidentifiable?.

    if (Nident) {
        id <- Inf
    } else {
        id  <- 1 / sqrt(min(eigen(t(S) %*% S)$value)) # collinearity
    }

    # output: psub= whether par is in (1) set or not (0), n= number of
    # parameters, id = the collinearity value
    psub     <- rep(0, npar)
    psub[ii] <- 1
    n        <- ncol(cc)
    Collset  <- rbind(Collset, c(psub, n, id))
  }
  return(Collset)
}

## =============================================================================
## Main function: Collinearity indices
## =============================================================================

collin <- function(sensfun, parset = NULL, N = NULL, which = NULL, maxcomb = 5000) {

  #-----------------
  # 1. check input
  if (is.null(colnames(sensfun)))
    colnames(sensfun) <- 1:ncol(sensfun)

  #-----------------
  # 2. If observed *variables* are specified ..
  if (!is.null(which)) {
    nx  <- attr(sensfun, "nx")
    var <- attr(sensfun, "var")
    TYP <- attr(sensfun, "Type")

    if (! is.numeric(which)) { #.. by name
      ln <- length(which)
      Select <- which (var %in% which)
      if(length(Select) != ln)
        stop("not all variables in 'which' are in 'sensfun'")
    } else {                   #.. by value
       Select <- which
       if (max(Select) > nx)
         stop("index in 'which' too large")
    }
    ii <- NULL

    if (TYP == 1)
      for (i in Select)  ii <- c(ii, ((i-1)*nx):(i*nx))
    else
      for (i in Select)  ii <- c(ii, (nx[i]+1):nx[i+1])
    # select relevant part of sensitivity functions
    sensfun <- sensfun[ii,]
  }

  # cleanup
  if (colnames(sensfun)[1]=="x" && colnames(sensfun)[2] == "var")
    Sens <- sensfun[,-(1:2)]
  else
    Sens <- sensfun

  npar <- ncol(Sens)
  L2   <- sqrt(colSums(Sens*Sens))

  iNa <- 0
  # Check for non-identifiable parameters
  if (any(L2 == 0 ) ) {
    iNa <- which(L2 == 0)
    warning (paste("Sensitivity of parameter", colnames(Sens)[iNa], "is 0! "))
  }

  # check if work requested not too large...
# TOGGLED OFF 16/04/2014
#  if (npar > 14 & is.null(parset) & is.null(N))
#    warning ("will reduce collinearity estimates: too many combinations")

  # normalise sensitivity functions
  normSens <- t(t(Sens) / L2)


  Collin <- NULL      # Will contain the collinearity coefficients

  if (is.null(parset)) {   # parameter combination not given, but the number of params
    pset   <- 1:npar

    if (is.null(N))        # All parameter combination (from 2 to total number)
      nset <- 2:npar
    else                   # only given number of parameters
      nset <- N

    for (n in nset) {
      numcomb <- choose(npar, n) # number of combinations requested
      if (numcomb < maxcomb) {
        cc  <- combin(n, pset)   # all combinations of n pars from pset.
        # collinearity of the parameter sets in cc
        Collin <- rbind(Collin, collFun(cc, normSens, npar, iNa))
      } else if (!is.null(N)) 
        stop ("too many parameter combinations")
      else {
       warning ("will not produce collinearity estimates for n = ", n)        
       warning ("number of combinations = ", numcomb, "maxcomb = ", maxcomb)        
      }      
    }
  } else {                 # parameter combination specified ..
    if (! is.vector(parset))
      stop("'parset' should be a vector")

    if (is.character(parset)) {  #.. if by name: find indices
      ln <- length(parset)
      pnames<-colnames(Sens)
      parset <- which(pnames %in% parset)
      if (length(parset) != ln)
        stop ("Not all parameters in parset known")
    }

    parset <- matrix(data = parset, nrow = 1)
    Collin <- rbind(Collin, collFun(parset, normSens, npar, iNa))
  }

  Collin <- as.data.frame(Collin)

  class(Collin) <- c("collin", "data.frame")
  names(Collin) <- c(colnames(Sens), "N", "collinearity")

  return(Collin)
}

## -----------------------------------------------------------------------------
## S3 methods of collin
## -----------------------------------------------------------------------------

plot.collin <- function(x, ...) {

  dots <- list(...)
  dots$ylab <- if(is.null(dots$ylab)) "Collinearity index" else dots$ylab
  dots$xlab <- if(is.null(dots$xlab)) "Number of parameters" else dots$xlab
  dots$main <- if(is.null(dots$main)) "Collinearity" else dots$main

  nc <- ncol(x)
  do.call("stripchart", c(alist(x[,nc] ~ x[,nc-1], method="stack",
    vertical = TRUE), dots))
}


## -----------------------------------------------------------------------------

print.collin <- function (x, ...)
  print(format(as.data.frame(unclass(x)), digits = 2, scientific = FALSE, ...))
