## RC(M) models with skew-symmetric association (Yamaguchi, 1990)

YRCSkew <- function(row, col, rowinf, rowsup, inst=NULL) {
  list(predictors = list(substitute(row), substitute(col), substitute(rowinf), substitute(rowsup)),
       term = function(predLabels, varLabels) {
           sprintf("%s * %s * (%s - %s) - %s * %s * (%s - %s)",
                   predLabels[3], predLabels[1], predLabels[2], predLabels[1],
                   predLabels[4], predLabels[2], predLabels[1], predLabels[2])
       },
       call = as.expression(match.call()),
       common = c(1, 1, 2, 2)
       )
  }
class(YRCSkew) <- "nonlin"

yrcskew <- function(tab, nd.symm=NA, nd.skew=1, diagonal=FALSE,
                    weighting=c("marginal", "uniform", "none"), se=c("none", "jackknife", "bootstrap"),
                    nreplicates=100, ncpus=getOption("boot.ncpus"),
                    family=poisson, weights=NULL, start=NA, etastart=NULL, tolerance=1e-8, iterMax=15000,
                    trace=FALSE, verbose=TRUE, ...) {
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(nrow(tab) != ncol(tab))
      stop("tab must be a square table for asymmetry models")

  if(!all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for symmetric model")

  if(!is.na(nd.symm) && nd.symm <= 0)
      stop("nd.symm must be NA or strictly positive")

  if(is.na(nd.skew))
      stop("nd.skew must be strictly positive")

  if(nd.skew > 1)
      warning("nd.skew > 1 is untested. You are on your own!")

  if(!is.na(nd.symm) && nd.symm/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric association cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(length(dim(tab)) > 2)
      tab <- margin.table(tab, 1:2)

  tab <- prepareTable(tab, TRUE)
  vars <- names(dimnames(tab))


  if(diagonal && !is.na(nd.symm))
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""

  if(is.na(nd.symm))
      basef <- sprintf("Freq ~ %s + %s %s+ Symm(%s, %s)",
                       vars[1], vars[2], diagstr, vars[1], vars[2])
  else
      basef <- sprintf("Freq ~ %s + %s %s+ instances(MultHomog(%s, %s), %i)",
                       vars[1], vars[2], diagstr, vars[1], vars[2], nd.symm)


  if(!is.null(start) && is.na(start)) {
      # Without good starting values, estimation can fail when running start-up iterations
      # with large tables
      cat("Running base model to find starting values...\n")

      args <- list(formula=as.formula(basef),
                   data=tab, weights=weights, family=family,
                   tolerance=1e-3, iterMax=iterMax, verbose=verbose, trace=trace)

      base <- do.call("gnm", c(args, list(...)))

      start <- c(parameters(base), rep(NA, nd.skew * (nrow(tab) + 2)))

      if(is.null(etastart))
          etastart <- as.numeric(predict(base))

      cat("Running real model...\n")
  }

  f <- sprintf("%s + instances(YRCSkew(%s, %s, factor(ifelse(as.numeric(%s) < as.numeric(%s), 1, 0)), factor(ifelse(as.numeric(%s) > as.numeric(%s), 1, 0))), %s)",
               basef,
               vars[1], vars[2],
               vars[1], vars[2], vars[1], vars[2], nd.skew)

  args <- list(formula=eval(as.formula(f)), data=tab,
               constrain="YRCSkew\\(.*\\)0$",
               family=family, weights=weights, start=start, etastart=etastart,
               tolerance=tolerance, iterMax=iterMax, verbose=verbose, trace=trace)

  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  class(model) <- c("yrcskew", "rc.symm", "rc", "assocmod", class(model))

  model$call.gnm <- model$call
  model$call <- match.call()

  if(!is.na(nd.symm))
      model$assoc <- assoc.rc.symm(model, weighting=weighting)
  else
      model$assoc <- list()

  model$assoc.yrcskew <- assoc.yrcskew(model, weighting=weighting)


  if(se %in% c("jackknife", "bootstrap")) {
      assoc1 <- if(is.na(nd.symm)) assoc.yrcskew else assoc.rc.symm
      assoc2 <- if(is.na(nd.symm)) NULL else assoc.yrcskew

      jb <- jackboot(se, ncpus, nreplicates, tab, model, assoc1, assoc2,
                     weighting, NULL, NULL, family, weights,
                     verbose, trace, start, etastart, ...)

      if(!is.na(nd.symm)) {
          model$assoc$covtype <- se
          model$assoc$covmat <- jb$covmat1
          model$assoc$adj.covmats <- jb$adj.covmats1
          model$assoc$jack.results <- jb$jack.results1
          model$assoc$boot.results <- jb$boot.results1

          model$assoc.yrcskew$covtype <- se
          model$assoc.yrcskew$covmat <- jb$covmat2
          model$assoc.yrcskew$adj.covmats <- jb$adj.covmats2
          model$assoc.yrcskew$jack.results <- jb$jack.results2
          model$assoc.yrcskew$boot.results <- jb$boot.results2
      }
      else {
          model$assoc.yrcskew$covtype <- se
          model$assoc.yrcskew$covmat <- jb$covmat
          model$assoc.yrcskew$adj.covmats <- jb$adj.covmats
          model$assoc.yrcskew$jack.results <- jb$jack.results
          model$assoc.yrcskew$boot.results <- jb$boot.results
      }
  }
  else {
      if(!is.na(nd.symm)) {
          model$assoc$covtype <- se
          model$assoc$covmat <- numeric(0)
          model$assoc$adj.covmats <- numeric(0)
          model$assoc$boot.results <- numeric(0)
          model$assoc$jack.results <- numeric(0)
      }

      model$assoc.yrcskew$covtype <- se
      model$assoc.yrcskew$covmat <- numeric(0)
      model$assoc.yrcskew$adj.covmats <- numeric(0)
      model$assoc.yrcskew$boot.results <- numeric(0)
      model$assoc.yrcskew$jack.results <- numeric(0)
  }

  model
}

## Model with separate skew-symmetric association (RC_SK)
assoc.yrcskew <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
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


  ## Get skew association coefficients
  sc <- matrix(NA, nrow(tab), 0)

  nd <- 0
  while(TRUE) {
      musk <- parameters(model)[pickCoef(model, sprintf("YRCSkew\\(.*inst = %i\\)(\\Q%s\\E)",
                                                  nd+1, paste(rownames(tab), collapse="\\E|\\Q")))]

      if(length(musk) != nrow(tab))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, musk)
  }

  if(nd <= 0) {
      musk <- parameters(model)[pickCoef(model, sprintf("YRCSkew.*\\)(\\Q%s\\E)$",
                                                  paste(rownames(tab), collapse="\\E|\\Q")))]

      if(length(musk) == nrow(tab)) {
          nd <- 1
          sc <- cbind(sc, musk)
      }
      else {
          stop("No skew dimensions found. Are you sure this is a Yamaguchi RC_SK model?")
      }
  }

  skew <- parameters(model)[pickCoef(model, "YRCSkew.*\\)1$")]
  if(length(skew) != nd)
      stop("Skew coefficients not found. Are you sure this is a Yamaguchi RC_SK model?")

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], ncol=nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd) {
    lambda <- lambda + abs(skew[i]) * sc[,i] %o% sc[,i]
  }
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  sv <- svd(lambda0)
  sc[] <- diag(1/sqrt(p)) %*% sv$u[,1:nd] # Eq. A.4.7
  # Preserve the sign of skew
  phisk <- sign(skew) * sv$d[1:nd]


  ## Prepare objects
  phisk <- rbind(c(phisk))
  dim(sc)[3] <- 1
  names(phisk) <- NULL
  colnames(sc) <- colnames(phisk) <- paste("Dim", 1:nd, sep="")
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

  row.weights <- as.matrix(apply(tab, 1, sum, na.rm=TRUE))
  col.weights <- as.matrix(apply(tab, 2, sum, na.rm=TRUE))

  obj <- list(phi = phisk, row = sc, col = sc,  diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)

  class(obj) <- c("assoc.yrcskew", "assoc.symm", "assoc")
  obj
}


## Model with skew-symmetric association homogeneous to symmetric association (HM_(S+SK))
# assoc.yrcskew.homog <- function(model, weighting=c("marginal", "unit"), ...) {
#   if(!"gnm" %in% class(model)) stop("model must be a gnm object")
# 
#   # gnm doesn't include coefficients for NA row/columns, so get rid of them too
#   tab <- as.table(model$data[!is.na(rownames(model$data)), !is.na(colnames(model$data))])
# 
#   weighting <- match.arg(weighting)
# 
#   # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
#   if(weighting == "marginal")
#       p <- prop.table(apply(tab, 1, sum, na.rm=TRUE) + apply(tab, 2, sum, na.rm=TRUE))
#   else
#       p <- rep(1/nrow(tab), nrow(tab))
# 
#   ## Get symmetric association coefficients
#   sc <- matrix(NA, nrow(tab), 0)
# 
#   nd <- 0
#   while(TRUE) {
#       mu <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(.*inst = %s\\)", nd+1))]
# 
#       if(!(length(mu) == nrow(tab)))
#           break
# 
#       # This is a correct dimension, add it
#       nd <- nd + 1
# 
#       sc <- cbind(sc, mu)
#   }
# 
#   if(nd <= 0) {
#       mu <- parameters(model)[pickCoef(model, "MultHomog\\(")]
# 
#       if(length(mu) == nrow(tab)) {
#           nd <- 1
#           sc <- cbind(sc, mu)
#       }
#       else {
#           stop("No dimensions found. Are you sure this is a row-column association model with symmetric row and column scores?")
#       }
#   }
# 
#   ## Get skew association coefficients
#   scsk <- matrix(NA, nrow(tab), 0)
# 
#   ndsk <- 0
#   while(TRUE) {
#       musk <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(.*skew.*inst = %i\\)\\)(\\Q%s\\E)",
#                                                   nd + 1, paste(rownames(tab), collapse="\\E|\\Q")))]
# 
#       if(!(length(musk) == nrow(tab)))
#           break
# 
#       # This is a correct dimension, add it
#       ndsk <- ndsk + 1
# 
#       scsk <- cbind(scsk, musk)
#   }
# 
#   if(ndsk <= 0) {
#       musk <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(.*skew.*\\)(\\Q%s\\E)$",
#                                                   paste(rownames(tab), collapse="\\E|\\Q")))]
# 
#       if(length(musk) == nrow(tab)) {
#           ndsk <- 1
#           scsk <- cbind(scsk, musk)
#       }
#       else {
#           stop("No skew dimensions found. Are you sure this is a row-column association model with symmetric row and column scores plus skewness?")
#       }
#   }
# 
#   skew <- parameters(model)[pickCoef(model, "MultHomog\\(.*skew\\)$")]
#   if(length(skew) != ndsk)
#       stop("skew coefficients not found. Are you sure this is a row-column association model with symmetric row and column scores plus skewness?")
# 
#   if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
#       dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], ncol=nrow(tab))
#   else if(length(pickCoef(model, "Diag\\(")) > 0)
#       dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
#   else
#       dg <- numeric(0)
# 
# 
#   ## Normalize, cf. Clogg & Shihadeh (1994), eq. 5.3 et 5.4 (p. 83)
#   ## Symmetric part
#   # Center
#   sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")
# 
#   # Technique proposed in Goodman (1991), Appendix 4
#   # We use eigen() here rather than svd() because matrix is square and symmetric
#   lambda <- matrix(0, nrow(tab), ncol(tab))
#   for(i in 1:nd)
#       lambda <- lambda + sc[,i] %o% sc[,i]
#   lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
#   eigen <- eigen(lambda0, symmetric=TRUE)
#   sc[,1:nd] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
#   phi <- eigen$values[1:nd]
# 
#   ## Skew part
#   # Center
#   scsk <- sweep(scsk, 2, colSums(sweep(scsk, 1, p, "*")), "-")
# 
#   # Integrate skew coefficient, changing sign if needed
#   if(skew < 0) {
#       skew <- -skew
#       scsk <- -scsk
#   }
#   scsk <- sweep(scsk, 2, sqrt(skew), "*")
# 
#   # Technique proposed in Goodman (1991), Appendix 4
#   # We use eigen() here rather than svd() because matrix is square and symmetric
#   lambda <- matrix(0, nrow(tab), ncol(tab))
#   for(i in 1:ndsk)
#       lambda <- lambda + scsk[,i] %o% scsk[,i]
#   lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
#   eigen <- eigen(lambda0, symmetric=TRUE)
#   scsk[,1:ndsk] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:ndsk] # Eq. A.4.7
#   phisk <- eigen$values[1:ndsk]
# 
# 
#   ## Prepare objects
#   phi <- rbind(c(phi))
#   dim(sc)[3] <- dim(scsk)[3] <- 1
#   colnames(sc) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
#   colnames(scsk) <- paste("Dim", 1:ndsk, sep="")
#   rownames(sc) <- rownames(scsk) <- rownames(tab)
# 
#   if(length(dg) > 0) {
#       # Diag() sorts coefficients alphabetically!
#       dg[,order(rownames(tab))] <- dg
# 
#       colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
#                       else paste(rownames(tab), colnames(tab), sep=":")
# 
#       if(nrow(dg) > 1)
#           rownames(dg) <- dimnames(tab)[[3]]
#       else
#           rownames(dg) <- "All levels"
#   }
# 
#   row.weights <- apply(tab, 1, sum, na.rm=TRUE)
#   col.weights <- apply(tab, 2, sum, na.rm=TRUE)
#
#   obj <- list(phi = phi, row = sc, col = sc, phisk = phi, rowsk = scsk, colsk = scsk,
#               diagonal = dg, weighting = weighting, row.weights = row.weights, col.weights = col.weights,
#               vars = vars)
# 
#   class(obj) <- c("assoc.yrcskew.homog", "assoc.symm", "assoc")
#   obj
# }

assoc <- function(model, weighting, ...) UseMethod("assoc", model)
