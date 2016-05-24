## RC(M) model

rc <- function(tab, nd=1, symmetric=FALSE, diagonal=FALSE,
               weighting=c("marginal", "uniform", "none"), rowsup=NULL, colsup=NULL,
               se=c("none", "jackknife", "bootstrap"), nreplicates=100, ncpus=getOption("boot.ncpus"),
               family=poisson, weights=NULL, start=NULL, etastart=NULL, tolerance=1e-8, iterMax=5000,
               trace=FALSE, verbose=TRUE, ...) {
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(is.na(nd) || nd <= 0)
      stop("nd must be strictly positive")

  if(symmetric && nrow(tab) != ncol(tab))
      stop("tab must be a square table for symmetric model")

  if(symmetric && !all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for symmetric model")

  if(symmetric && nd/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!symmetric && nd > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  if(length(dim(tab)) > 2)
      tab <- margin.table(tab, 1:2)

  tab <- prepareTable(tab, TRUE, rowsup, colsup)
  vars <- names(dimnames(tab))


  if(diagonal)
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""


  nastart <- length(start) == 1 && is.na(start)

  if(nastart) {
      cat("Running base model to find starting values...\n")

      args <- list(formula=as.formula(sprintf("Freq ~ %s + %s %s", vars[1], vars[2], diagstr)),
                   data=tab, family=family, weights=weights,
                   tolerance=tolerance, iterMax=iterMax)

      base <- do.call("gnm", c(args, list(...)))

      res <- if(symmetric) residEVD(base, nd)
             else residSVD(base, nd)

      # Using NA for all linear parameters usually gives better results than linear parameter values
      start <- c(rep(NA, length(parameters(base))), res)

      if(is.null(etastart))
          etastart <- as.numeric(predict(base))

      cat("Running real model...\n")
  }

  if(symmetric)
      f <- sprintf("Freq ~ %s + %s %s+ instances(MultHomog(%s, %s), %i)",
                   vars[1], vars[2], diagstr, vars[1], vars[2], nd)
  else
      f <- sprintf("Freq ~ %s + %s %s+ instances(Mult(%s, %s), %i)",
                   vars[1], vars[2], diagstr, vars[1], vars[2], nd)

  args <- list(formula=as.formula(f), data=tab,
               family=family, weights=weights, start=start, etastart=etastart,
               tolerance=tolerance, iterMax=iterMax, verbose=verbose, trace=trace)

  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  newclasses <- if(symmetric) c("rc.symm", "rc", "assocmod") else c("rc", "assocmod")
  class(model) <- c(newclasses, class(model))

  model$call.gnm <- model$call
  model$call <- match.call()

  model$assoc <- assoc(model, weighting=weighting, rowsup=rowsup, colsup=colsup)


  if(se %in% c("jackknife", "bootstrap")) {
      jb <- jackboot(se, ncpus, nreplicates, tab, model,
                     assoc1=getS3method("assoc", class(model)), assoc2=NULL,
                     weighting, rowsup=rowsup, colsup=colsup,
                     family, weights, verbose, trace, start, etastart, ...)
      model$assoc$covmat <- jb$covmat
      model$assoc$adj.covmats <- jb$adj.covmats
      model$assoc$boot.results <- jb$boot.results
      model$assoc$jack.results <- jb$jack.results
  }
  else {
      model$assoc$covmat <- numeric(0)
      model$assoc$adj.covmats <- numeric(0)
      model$assoc$boot.results <- numeric(0)
      model$assoc$jack.results <- numeric(0)
  }

  model$assoc$covtype <- se

  model
}

assoc.rc <- function(model, weighting=c("marginal", "uniform", "none"),
                     rowsup=NULL, colsup=NULL, ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  tab <- prepareTable(model$data, TRUE, rowsup, colsup)
  vars <- names(dimnames(tab))

  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal") {
      rp <- prop.table(apply(tab, 1, sum, na.rm=TRUE))
      cp <- prop.table(apply(tab, 2, sum, na.rm=TRUE))
  }
  else if(weighting == "uniform") {
      rp <- rep(1/nrow(tab), nrow(tab))
      cp <- rep(1/ncol(tab), ncol(tab))
  }
  else {
      rp <- rep(1, nrow(tab))
      cp <- rep(1, ncol(tab))
  }


  # Prepare matrices before filling them
  row <- matrix(NA, nrow(tab), 0)
  col <- matrix(NA, ncol(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\., \\Q%s\\E, inst = %i\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[2], nd + 1, vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\Q%s\\E, \\., inst = %i\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], nd + 1, vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]

      if(!(length(mu) == nrow(tab) && length(nu) == ncol(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      row <- cbind(row, mu)
      col <- cbind(col, nu)
  }

  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\., \\Q%s\\E)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[2], vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\Q%s\\E, \\.)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]

      if(length(mu) == nrow(tab) && length(nu) == ncol(tab)) {
          nd <- 1

          row <- cbind(row, mu)
          col <- cbind(col, nu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model?")
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], ncol=nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  # Center
  row <- sweep(row, 2, colSums(sweep(row, 1, rp/sum(rp), "*")), "-")
  col <- sweep(col, 2, colSums(sweep(col, 1, cp/sum(cp), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd) {
      lambda <- lambda + row[,i] %o% col[,i]
  }
  lambda0 <- lambda * sqrt(rp %o% cp) # Eq. A.4.3
  sv <- svd(lambda0)
  row[] <- diag(1/sqrt(rp)) %*% sv$u[,1:nd] # Eq. A.4.7
  col[] <- diag(1/sqrt(cp)) %*% sv$v[,1:nd] # Eq. A.4.7
  phi <- sv$d[1:nd]

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(row[1,i] < 0) {
          row[,i] <- -row[,i]
          col[,i] <- -col[,i]
      }
  }


  ## Prepare objects
  phi <- rbind(c(phi))
  dim(row)[3] <- dim(col)[3] <- 1
  colnames(row) <- colnames(col) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)

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

  obj <- list(phi = phi, row = row, col = col, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)


  ## Supplementary rows/columns
  if(!is.null(rowsup) || !is.null(colsup)) {
      sup <- sup.scores.rc(model, tab, obj, rowsup, colsup)

      obj$row <- sup$row
      obj$col <- sup$col
      obj$row.weights <- sup$row.weights
      obj$col.weights <- sup$col.weights
  }

  class(obj) <- c("assoc.rc", "assoc")
  obj
}

## RC(M) model with symmetric row and column scores
assoc.rc.symm <- function(model, weighting=c("marginal", "uniform", "none"),
                          rowsup=NULL, colsup=NULL, ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  tab <- prepareTable(model$data, TRUE, rowsup, colsup)
  vars <- names(dimnames(tab))

  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(apply(tab, 1, sum, na.rm=TRUE) + apply(tab, 2, sum, na.rm=TRUE))
  else if(weighting == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))


  sc <- matrix(NA, nrow(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(\\Q%s\\E\\, \\Q%s\\E\\, inst = %i\\)(\\Q%s\\E)$",
                                                vars[1], vars[2], nd + 1,
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]

      if(!(length(mu) == nrow(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, mu)
  }

  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(\\Q%s\\E\\, \\Q%s\\E\\)(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]

      if(length(mu) == nrow(tab)) {
          nd <- 1
          sc <- cbind(sc, mu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with symmetric row and column scores?")
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], ncol=nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4, but with eigenvalues decomposition
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd)
      lambda <- lambda + sc[,i] %o% sc[,i]
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  eigen <- eigen(lambda0, symmetric=TRUE)
  sc[,1:nd] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
  phi <- eigen$values[1:nd]

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(sc[1,i] < 0)
          sc[,i] <- -sc[,i]
  }

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(sc)[3] <- 1
  colnames(sc) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
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

  obj <- list(phi = phi, row = sc, col= sc, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)

  ## Supplementary rows/columns
  if(!is.null(rowsup) || !is.null(colsup)) {
      sup <- sup.scores.rc(model, tab, obj, rowsup, colsup, symmetry="symmetric")

      obj$row <- sup$row
      obj$col <- sup$col
      obj$row.weights <- sup$row.weights
      obj$col.weights <- sup$col.weights
  }

  class(obj) <- c("assoc.rc", "assoc.symm", "assoc")
  obj
}

assoc <- function(model, weighting, rowsup, colsup, ...) UseMethod("assoc", model)
