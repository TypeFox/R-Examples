## Skew-symmetric association model (van der Heijden & Mooijaart, 1995)

HMSkew <- function(row, col, inst=NULL) {
  list(predictors = list(substitute(row), substitute(col), substitute(row), substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("%s * %s - %s * %s",
                   predLabels[3], predLabels[2], predLabels[1], predLabels[4])
       },
       call = as.expression(match.call()),
       common = c(1, 1, 2, 2)
       )
   }
class(HMSkew) <- "nonlin"

hmskew <- function(tab, nd.symm=NA, diagonal=FALSE,
                   weighting=c("marginal", "uniform", "none"),  rowsup=NULL, colsup=NULL,
                   se=c("none", "jackknife", "bootstrap"), nreplicates=100, ncpus=getOption("boot.ncpus"),
                   family=poisson, weights=NULL, start=NULL, etastart=NULL, tolerance=1e-8, iterMax=5000,
                   trace=FALSE, verbose=TRUE, ...) {
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(nrow(tab) != ncol(tab))
      stop("tab must be a square table for asymmetry models")

  if(!all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for asymmetry models")

  if(!is.na(nd.symm) && nd.symm < 0)
      stop("nd.symm must be NA, zero or positive")

  if(!is.na(nd.symm) && nd.symm/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric association cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!is.null(rowsup) && !is.matrix(rowsup))
      stop("'rowsup' must be a matrix")

  if(!is.null(colsup) && !is.matrix(colsup))
      stop("'colsup' must be a matrix")

  if(!is.null(rowsup) && ncol(rowsup) != ncol(tab))
      stop("'rowsup' must have one column for each column in 'tab'")

  if(!is.null(colsup) && nrow(colsup) != nrow(tab))
      stop("'colsup' must have one row for each row in 'tab'")

  if(length(dim(tab)) > 2)
      tab <- margin.table(tab, 1:2)

  tab <- prepareTable(tab, TRUE, rowsup, colsup)
  vars <- names(dimnames(tab))


  if(diagonal && !is.na(nd.symm))
      diagstr <- sprintf(" + Diag(%s, %s)", vars[1], vars[2])
  else
      diagstr <- ""

  basef <- sprintf("Freq ~ %s + %s%s", vars[1], vars[2], diagstr)

  if(is.na(nd.symm))
      basef2 <- sprintf("%s + Symm(%s, %s)",
                        basef, vars[1], vars[2])
  else if(nd.symm > 0)
      basef2 <- sprintf("%s + instances(MultHomog(%s, %s), %i)",
                        basef, vars[1], vars[2], nd.symm)
  else
      basef2 <- basef

  if(!is.null(start) && is.na(start)) {
      cat("Running base model to find starting values...\n")

      # Setting tolerance to a value below 1e-6 can lead to convergence issues with large tables
      args <- list(formula=as.formula(basef),
                   data=tab, family=family, weights=weights,
                   tolerance=1e-3, iterMax=iterMax, verbose=verbose, trace=trace)

      base <- do.call("gnm", c(args, list(...)))

      args$method <- "coefNames"
      args$formula <- as.formula(basef2)
      npar <- length(do.call("gnm", args))

      # Use base model without symmetric interaction as this seems to give better results
      # when quasi-symmetry is specified.
      # Using NA for all linear parameters usually works better than their values in the base model
      start <- c(rep(NA, npar), residEVD(base, 1, skew=TRUE))

      if(is.null(etastart))
          etastart <- as.numeric(predict(base))

      cat("Running real model...\n")
  }

  f <- sprintf("%s + HMSkew(%s, %s)", basef2, vars[1], vars[2])

  args <- list(formula=as.formula(f), data=tab,
               family=family, weights=weights, start=start, etastart=etastart,
               tolerance=tolerance, iterMax=iterMax, verbose=verbose, trace=trace)

  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  if(!is.na(nd.symm) && nd.symm > 0) {
      model$assoc <- assoc.rc.symm(model, weighting=weighting, rowsup=rowsup, colsup=colsup)
      class(model) <- c("hmskew", "rc.symm", "rc", "assocmod", class(model))
  }
  else {
      class(model) <- c("hmskew", "assocmod", class(model))
  }

  model$assoc.hmskew <- assoc.hmskew(model, weighting=weighting, rowsup=rowsup, colsup=colsup)

  model$call.gnm <- model$call
  model$call <- match.call()


  if(se %in% c("jackknife", "bootstrap")) {
      assoc1 <- if(!is.na(nd.symm) && nd.symm > 0) assoc.rc.symm else assoc.hmskew
      assoc2 <- if(!is.na(nd.symm) && nd.symm > 0) assoc.hmskew else NULL

      jb <- jackboot(se, ncpus, nreplicates, tab, model, assoc1, assoc2,
                     weighting,  rowsup=rowsup, colsup=colsup,
                     family, weights, verbose, trace, start, etastart, ...)

      if(!is.na(nd.symm) && nd.symm > 0) {
          model$assoc$covtype <- se
          model$assoc$covmat <- jb$covmat1
          model$assoc$adj.covmats <- jb$adj.covmats1
          model$assoc$jack.results <- jb$jack.results1
          model$assoc$boot.results <- jb$boot.results1

          model$assoc.hmskew$covtype <- se
          model$assoc.hmskew$covmat <- jb$covmat2
          model$assoc.hmskew$adj.covmats <- jb$adj.covmats2
          model$assoc.hmskew$jack.results <- jb$jack.results2
          model$assoc.hmskew$boot.results <- jb$boot.results2
      }
      else {
          model$assoc.hmskew$covtype <- se
          model$assoc.hmskew$covmat <- jb$covmat
          model$assoc.hmskew$adj.covmats <- jb$adj.covmats
          model$assoc.hmskew$jack.results <- jb$jack.results
          model$assoc.hmskew$boot.results <- jb$boot.results
      }
  }
  else {
      if(!is.na(nd.symm) && nd.symm > 0) {
          model$assoc$covtype <- se
          model$assoc$covmat <- numeric(0)
          model$assoc$adj.covmats <- numeric(0)
          model$assoc$boot.results <- numeric(0)
          model$assoc$jack.results <- numeric(0)
      }

      model$assoc.hmskew$covtype <- se
      model$assoc.hmskew$covmat <- numeric(0)
      model$assoc.hmskew$adj.covmats <- numeric(0)
      model$assoc.hmskew$boot.results <- numeric(0)
      model$assoc.hmskew$jack.results <- numeric(0)
  }

  model
}

assoc.hmskew <- function(model, weighting=c("marginal", "uniform", "none"),
                         rowsup=NULL, colsup=NULL, ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  tab <- prepareTable(model$data, TRUE)
  vars <- names(dimnames(tab))

  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(apply(tab, 1, sum, na.rm=TRUE) + apply(tab, 2, sum, na.rm=TRUE))
  else if(weighting == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))


  mu <- parameters(model)[pickCoef(model, sprintf("HMSkew.*(\\Q%s\\E)",
                                            paste(rownames(tab), collapse="\\E|\\Q")))]
  mu1 <- mu[1:nrow(tab)]
  mu2 <- mu[-(1:nrow(tab))]

  if(!(length(mu1) == nrow(tab) && length(mu2) == nrow(tab)))
      stop("No dimensions found. Are you sure this is a van der Heijden & Mooijaart skewness model?")


  sc <- cbind(mu1, mu2)

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], dim(tab)[3], nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  ## Normalize, cf. Clogg & Shihadeh (1994), eq. 5.3 et 5.4 (p. 83)
  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")
  # Technique proposed in van der Heijden & Mooijaart (1995), Appendix C
  # Weighted SVD taken from Goodman (1991), Appendix 4
  lambda <- sc[,2] %o% sc[,1] - sc[,1] %o% sc[,2]
  lambda0 <- sqrt(p %o% p) * lambda # Eq. A.4.3
  sv <- svd(lambda0)
  sc[,1:2] <- diag(1/sqrt(p)) %*% sv$u[,1:2] # Eq. A.4.7
  phi <- sv$d[1:2]

  # Since both dimensions share the same singular value, we cannot distinguish them
  # and their order is random. The sign of the association reconstructed from scores
  # must be the same as the original one: we use one cell, since all signs are the same.
  # This operation affects the results, while the next one is merely a display convention.
  sc[, 1] <- sc[, 1] * sign((sc[1, 2] * sc[2, 1] - sc[1, 1] * sc[2, 2])/lambda[1, 2])

  # Use the convention that we want the null score to be for the first category
  # on the second dimension, and a positive score on the first dimension
  # (like original article does with *last* category).
  # The SVD always sets to 0 the score of the first category on one dimension.
  # These two operations do not change the actual association.
  if(which.min(abs(sc[1,])) == 1)
      sc <- sc %*% matrix(c(0, 1, -1, 0), 2, 2)

  if(sc[1, 1] < 0)
      sc <- -sc

  # The reference category is not really at 0, and this makes the display ugly
  sc[abs(sc) < sqrt(.Machine$double.eps)] <- 0

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(sc)[3] <- 1
  colnames(sc) <- colnames(phi) <- paste("Dim", 1:2, sep="")
  rownames(sc) <- rownames(tab)

  if(length(dg) > 0) {
      # Diag() sorts coefficients alphabetically!
      dg[,order(rownames(tab))] <- dg

      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")

      if(nrow(dg) > 1)
          rownames(dg) <- dimnames(tab)[[3]]
      else
          rownames(dg) <- "All levels"
  }

  if(length(dim(tab)) == 3) {
      row.weights <- apply(tab, c(1, 3), sum, na.rm=TRUE)
      col.weights <- apply(tab, c(2, 3), sum, na.rm=TRUE)
  }
  else {
      row.weights <- as.matrix(apply(tab, 1, sum, na.rm=TRUE))
      col.weights <- as.matrix(apply(tab, 2, sum, na.rm=TRUE))
  }

  obj <- list(phi = phi, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)

  ## Supplementary rows/columns
  if(!is.null(rowsup) || !is.null(colsup)) {
      sup <- sup.scores.rc(model, tab, obj, rowsup, colsup, symmetry="skew-symmetric", str="HMSkew")

      obj$row <- sup$row
      obj$col <- sup$col
      obj$row.weights <- sup$row.weights
      obj$col.weights <- sup$col.weights
      obj$row.sup <- sup$row.sup
      obj$col.sup <- sup$col.sup
  }

  class(obj) <- c("assoc.hmskew", "assoc.symm", "assoc")
  obj
}

assoc <- function(model, weighting, rowsup, colsup, ...) UseMethod("assoc", model)
