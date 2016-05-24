##' Compute the standard errors for the Bayes factors estimates.
##'
##' Uses the batch means method to compute the standard errors for
##' Bayes factors.
##' @title Computation of standard errors for Bayes factors
##' @param pargrid A data frame with components "linkp", "phi", "omg",
##' "kappa". Each row gives a combination of the parameters to compute
##' the new standard errors.
##' @param runs A list with outputs from the function
##' \code{\link{mcsglmm}} or \code{\link{mcstrga}}.
##' @param nbatches An integer scalar or vector of the same length as
##' runs indicating the number of batches to create for computing the
##' variance using the samples of the first stage.
##' @param bfsize1 A scalar or vector of the same length as
##' \code{runs} with all integer values or all values in (0, 1]. How
##' many samples (or what proportion of the sample) to use for
##' estimating the Bayes factors at the first stage. The remaining
##' sample will be used for estimating the standard errors in the
##' second stage. Setting it to 1 will perform only the first stage.
##' @param method Which method to use to calculate the Bayes factors:
##' Reverse logistic or Meng-Wong.
##' @param reference Which model goes in the denominator.
##' @param transf Whether to use the transformed sample mu for the
##' computations. Otherwise it uses z.
##' @param bmmcse.size Size for computing the batch means Monte-Carlo
##' variance using the samples of the second stage.
##' @return A list with components
##' \itemize{
##' \item \code{pargrid} The inputted pargrid augmented with the computed standard
##' errors.
##' \item \code{bfEstimate} The estimates of the Bayes factors
##' \item \code{bfSigma} The covariance matrix for the Bayes factors
##' estimates.
##' }
##' @references Roy, V., Tan, A. and Flegal, J. (2015). Estimating
##' standard errors for importance sampling estimators with multiple
##' Markov chains. Technical report, Iowa State University.
##' \url{http://lib.dr.iastate.edu/stat_las_preprints/34}
##' @export 
bfse <- function(pargrid, runs, nbatches, bfsize1 = 0.80,
                 method = c("RL", "MW"), reference = 1, transf = FALSE, 
                 bmmcse.size = "sqroot")
{
  ## Check input
  if (any(sapply(runs, class) != "geomcmc")) {
    stop ("Input runs is not a list with elements of class geomcmc")
  }
  nruns <- length(runs)
  nbatches <- rep(as.integer(nbatches), length.out = length(runs))
  if (any(nbatches <= 1)) stop ("Argument nbatches must >= 2")

  out <- list() # Function output

  ## Compute Bayes factors
  bf1call <- bf1skel(runs, bfsize1, method, reference, transf)
  Nout1 <- bf1call$N1
  Ntot1 <- sum(Nout1)
  logbf <- bf1call$logbf
  dvec <- exp(logbf)
  runsi <- rep(seq_len(nruns), Nout1)
  zeta <- -logbf + log(Nout1/Ntot1)
  logYYnum <- bf1call$logLik1 + matrix(zeta, Ntot1, nruns, TRUE)
  logYYdff <- logYYnum - apply(logYYnum, 1, max)
  logYY <- logYYdff - log(rowSums(exp(logYYdff)))
  YY <- exp(logYY)               # r is columns
  ## WW <- groupColMeans(YY, runsi) # r is columns

  ## Split to batches
  batchl <- mapply(function(n, e) {ne <- as.integer(n/e); ner <- n%%e;
                                   rep(c(ne + 1L, ne), c(ner, e - ner))}, 
                   Nout1, nbatches, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  batchi <- lapply(batchl, function(b) rep(seq_along(b), b))
  YYbar <- mapply(groupColMeans, split(data.frame(YY), rep(1:nruns, Nout1)),
                  batchi, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  WW <- t(sapply(1:nruns, function(i) colSums(YYbar[[i]]*batchl[[i]])/Nout1[i]))
  
  ## Compute matrix Omega
  Omglist <- sapply(1:nruns, function(l) {
    YW <- YYbar[[l]] - matrix(WW[l, ], nbatches[l], nruns, TRUE)
    rowSums(apply(YW*sqrt(batchl[[l]]), 1, tcrossprod))/(nbatches[l]-1)*
      (Nout1[l]/Ntot1)}) # k^2 rows
  Omg <- matrix(rowSums(Omglist), nruns, nruns)

  ## Compute matrix Beta
  Bet <- diag(colMeans(YY*(1-YY)), nruns, nruns)
  Bet[lower.tri(Bet)] <- -colMeans(YY[, rep(1:(nruns-1), (nruns-1):1)] *
                    YY[, unlist(lapply(2:nruns, function(i) seq.int(i, nruns)))])
  Bet[upper.tri(Bet)] <- t(Bet)[upper.tri(Bet)]
  ## Compute M-P inverse
  BMP <- chol2inv(chol(Bet + 1/nruns)) - 1/nruns

  ## Compute Sigma
  BD <- (BMP[, reference] - BMP[, -reference, drop = FALSE]) *
    matrix(dvec[-reference], nruns, nruns - 1L, TRUE)
  Sigma <- t(BD) %*% Omg %*% BD

  ## Check for Ntot2 == 0
  Nout2 <- bf1call$N2
  Ntot2 <- sum(Nout2)
  if (Ntot2 == 0) {
    pargrid$SE <- NA
    out$pargrid <- pargrid
    out$bfEstimate <- dvec
    out$bfSigma <- Sigma
    return (out)
  }

  ## Extract the grid for the new parameters
  family <- bf1call$family
  corrfcn <- bf1call$corrfcn
  pargrid <- .check_pargrid(pargrid, family, corrfcn)
  phi <- pargrid$phi
  omg <- pargrid$omg
  kappa <- pargrid$kappa
  nu <- pargrid$nu
  kg <- NROW(pargrid)

  ## Compute log-posterior for the new grid
  if (transf) {
    froutine <- "llikfcnmu"
  } else {
    froutine <- "llikfcnz"
  }
  sample <- bf1call$sample2
  y <- bf1call$response
  l <- bf1call$weights
  F <- bf1call$modelmatrix
  dm <- bf1call$distmat
  betm0 <- bf1call$betm0
  betQ0 <- bf1call$betQ0
  ssqdf <- bf1call$ssqdf
  ssqsc <- bf1call$ssqsc
  tsqdf <- bf1call$tsqdf
  tsqsc <- bf1call$tsqsc
  dispersion <- bf1call$dispersion
  n <- NROW(F)
  p <- NCOL(F)
  icf <- match(corrfcn, c("matern", "spherical", "powerexponential"))
  ifam <- match(family, eval(formals(mcsglmm)$family), 0L)
  lglk <- matrix(0, Ntot2, kg)
  tsq <- if (ifam == 0) tsqsc else dispersion 
  fcall <- .Fortran(froutine,
                    lglk = lglk,
                    as.double(phi), as.double(omg), as.double(nu),
                    as.double(kappa),
                    as.double(sample), as.integer(Ntot2), as.double(y),
                    as.double(l), as.double(F), as.double(dm),
                    as.double(betm0), as.double(betQ0), as.double(ssqdf),
                    as.double(ssqsc), max(tsqdf, 0), as.double(tsq),
                    as.integer(icf), as.integer(n), as.integer(p),
                    as.integer(kg), as.integer(ifam))
  lglk <- fcall$lglk

  ## Compute vector V
  llik2 <- bf1call$logLik2
  mxll2 <- apply(llik2, 1, max)
  llik2mm <- llik2 - mxll2
  logVVden <- log(rowSums(exp(llik2mm + matrix(log(Nout2) - logbf,
                      Ntot2, nruns, TRUE)))) + mxll2
  logVVN <- lglk - logVVden   # kg columns
  VV <- Ntot2*exp(logVVN)     # kg columns

  ## 2nd term in the variance
  secalc <- tapply(1:Ntot2, rep(1:nruns, Nout2),
                   function(jj) bmmcse(VV[jj, , drop = FALSE],
                                       size = bmmcse.size),
                   simplify = FALSE)
  semat <- matrix(unlist(secalc), nruns, kg, TRUE)
  VT2 <- colSums(Nout2*semat^2)/Ntot2

  ## Compute vector C
  logVhN <- llik2[, -reference, drop = FALSE] - logVVden
  logCterms <- logVhN[, rep(1:(nruns-1), each = kg), drop = FALSE] +
    logVVN[, rep(1:kg, nruns - 1), drop = FALSE]
  logCmx <- apply(logCterms, 2, max)
  logCCN <- matrix(logCmx + log(colSums(exp(logCterms -
                   matrix(logCmx, Ntot2, (nruns - 1)*kg, TRUE)))),
                   nruns - 1, kg)
  CC <- exp(logCCN)*(Nout2[-reference]/dvec[-reference]^2)

  ## 1st term in the variance
  VT1 <- (Ntot2/Ntot1)*colSums(c(Sigma) *
          CC[rep(1:(nruns-1), each = nruns-1), , drop = FALSE] *
          CC[rep(1:(nruns-1), nruns-1), , drop = FALSE])

  ## Return
  pargrid$SE <- sqrt((VT1 + VT2)/Ntot2)
  out$pargrid <- pargrid
  out$bfEstimate <- dvec
  out$bfSigma <- Sigma/Ntot1
  out$VT1 <- VT1/Ntot2
  out$VT2 <- VT2/Ntot2
  return (out)
}


groupColMeans <- function (x, i) {
  ## Compute the mean at each column of x by subsetting its rows
  ## according to the factor i. Returns a matrix with as many columns
  ## as x and as many rows as the number of levels in i.
  x <- as.matrix(x)
  rx <- nrow(x)
  cx <- ncol(x)
  m <- tapply(seq_len(rx), i, function(jj) colMeans(x[jj, , drop = FALSE]),
              simplify = FALSE)
  matrix(unlist(m), ncol = cx, byrow = TRUE)
}

bmmcse <- function(x, size)
{
  ## Purpose: Compute batch-means standard error of univariate MC samples
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x   A matrix of MC samples. The SE corresponding to each
  ## column of x will be returned.
  ## size   The size of each batch.
  ## ----------------------------------------------------------------------
  x <- as.matrix(x)
  n <- nrow(x)
  if (missing(size)) {
    size <- "sqroot"
  } else if (length(size) != 1) stop ("Argument size must have length 1.")
  if (is.character(size)) {
    i <- pmatch(size, c("sqroot", "sqrt", "cuberoot", "cubert"), nomatch = 0L,
                duplicates.ok = FALSE)
    if (i == 1 | i == 2) {
      size <- as.integer(sqrt(n))
    } else if (i == 3 | i == 4) {
      size <- as.integer(n^(1/3))
    } else stop ("Character size must be one of sqroot or cuberoot.")
  } else {
    size <- as.numeric(size)
    if (size <= 0) {
      stop ("Argument size must be positive.")
    } else if (size <= 1) {
      size <- as.integer(size*n)
    } else if (size > n) stop ("Argument size is larger than sample size.")
  }
  ## Determine the number of batches and the size of each batch
  nbatches <- as.integer(n/size)
  if (nbatches <= 1L) stop ("Number of batches too small. Reduce size.")
  nrem <- n%%nbatches
  batchl <- rep(c(size + 1L, size), c(nrem, nbatches - nrem))
  batchf <- rep(seq_len(nbatches), batchl)
  ## Compute batch mean and overall mean
  bm <- groupColMeans(x, batchf)
  m <- colSums(bm*batchl)/n
  ## Compute the batch variance
  bmmm <- sweep(bm, 2, m, check.margin = FALSE)
  bv <- colSums(batchl*bmmm^2)/(nbatches - 1)
  sqrt(bv/n)
}
