## RC(M) with transitional layer effect

RCTrans <- function(row, col, layer, inst=NULL) {
  list(predictors = list(R1=substitute(row), C1=substitute(col), substitute(layer), R2=substitute(row), C2=substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("(%s + %s * (%s)^2) * (%s + %s * (%s)^2)",
                   predLabels[1], predLabels[4], predLabels[3],
                   predLabels[2], predLabels[5], predLabels[3])
       },
       call = as.expression(match.call()),
       match = c(1, 2, 3, 1, 2),
       start = function(theta) {
           theta[attr(theta, "assign") == 3] <- sqrt(seq(0, 1, length.out=sum(attr(theta, "assign") == 3)))
           theta
       })
   }
class(RCTrans) <- "nonlin"

RCTransSymm <- function(row, col, layer, inst=NULL) {
  list(predictors = list(R1=substitute(row), C1=substitute(col), substitute(layer), R2=substitute(row), C2=substitute(col)),
       term = function(predLabels, varLabels) {
           sprintf("(%s + %s * (%s)^2) * (%s + %s * (%s)^2)",
                   predLabels[1], predLabels[4], predLabels[3],
                   predLabels[2], predLabels[5], predLabels[3])
       },
       call = as.expression(match.call()),
       match = c(1, 2, 3, 1, 2),
       common = c(1, 1, 3, 2, 2),
       start = function(theta) {
           theta[attr(theta, "assign") == 3] <- sqrt(seq(0, 1, length.out=sum(attr(theta, "assign") == 3)))
           theta
       })
   }
class(RCTransSymm) <- "nonlin"


rcL.trans <- function(tab, nd=1, symmetric=FALSE, diagonal=c("none", "heterogeneous", "homogeneous"),
                      weighting=c("marginal", "uniform", "none"), se=c("none", "jackknife", "bootstrap"),
                      nreplicates=100, ncpus=getOption("boot.ncpus"),
                      family=poisson, weights=NULL, start=NULL, etastart=NULL, tolerance=1e-8, iterMax=5000,
                      trace=FALSE, verbose=TRUE, ...) {
  diagonal <- match.arg(diagonal)
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

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

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)

  tab <- prepareTable(tab, FALSE)
  vars <- names(dimnames(tab))


  if(diagonal == "homogeneous")
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else if(diagonal == "heterogeneous")
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])
  else
      diagstr <- ""

  f1 <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s",
                vars[1], vars[2], vars[3],
                vars[1], vars[3], vars[2], vars[3], diagstr)

  eliminate <- eval(parse(text=sprintf("quote(%s:%s)", vars[1], vars[3])))

  if(!is.null(start) && is.na(start)) {
      cat("Running base model to find starting values...\n")

      if(symmetric) {
          args <- list(formula=as.formula(sprintf("%s + instances(MultHomog(%s, %s), %i)", f1, vars[1], vars[2], nd)),
                       data=tab, family=family, weights=weights, eliminate=eliminate,
                       tolerance=1e-6, iterMax=iterMax)

          base <- do.call("gnm", c(args, list(...)))

          start <- parameters(base)[seq(1, length(parameters(base)) - nrow(tab) * nd)]
          for(i in 1:nd)
              start <- c(start, parameters(base)[pickCoef(base, paste("Mult.*inst =", i))],
                         sqrt(seq(0, 1, length.out=dim(tab)[3])), rep(NA, nrow(tab)))
      }
      else {
          args <- list(formula=as.formula(sprintf("%s + instances(Mult(%s, %s), %i)", f1, vars[1], vars[2], nd)),
                       data=tab, family=family, weights=weights, eliminate=eliminate,
                       tolerance=1e-6, iterMax=iterMax, verbose=verbose, trace=trace)

          base <- do.call("gnm", c(args, list(...)))

          start <- parameters(base)[seq(1, length(parameters(base)) - (nrow(tab) + ncol(tab)) * nd)]
          for(i in 1:nd)
              start <- c(start, parameters(base)[pickCoef(base, paste("Mult.*inst =", i))],
                         sqrt(seq(0, 1, length.out=dim(tab)[3])), rep(NA, nrow(tab) + ncol(tab)))
      }

      if(!is.null(etastart))
          etastart <- as.numeric(predict(base))

      cat("Running real model...\n")
    }

  f <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s+ instances(%s(%s, %s, %s), %i)",
               vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3], diagstr,
               if(symmetric) "RCTransSymm" else "RCTrans", vars[1], vars[2], vars[3], nd)

  args <- list(formula=as.formula(f), data=tab, family=family, weights=weights,
               # Both constraints are really needed: the computed scores are wrong without them
               # (rows are ordered along an oblique axis, and columns get weird values)
               constrain=sprintf("RCTrans.*\\).\\Q%s\\E(\\Q%s\\E|\\Q%s\\E)$",
                                 vars[3], head(dimnames(tab)[[3]], 1), tail(dimnames(tab)[[3]], 1), nd),
               constrainTo=rep(0:1, nd),
               start=start, etastart=etastart, eliminate=eliminate,
               tolerance=tolerance, iterMax=iterMax, verbose=verbose, trace=trace)

  model <- do.call("gnm", c(args, list(...)))


  if(is.null(model))
      return(NULL)

  newclasses <- if(symmetric) c("rcL.trans.symm", "rcL.trans", "rcL", "assocmod")
                else c("rcL.trans", "rcL", "assocmod")
  class(model) <- c(newclasses, class(model))

  model$call.gnm <- model$call
  model$call <- match.call()

  model$assoc <- assoc(model, weighting=weighting)


 if(se %in% c("jackknife", "bootstrap")) {
      jb <- jackboot(se, ncpus, nreplicates, tab, model, assoc, NULL,
                     weighting, NULL, NULL, family, weights,
                     verbose, trace, start, etastart, ...)
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

# Number of constraints applied:
#   - two layer coefficients (normally already present in the model)
#   - two for each dimension of the start scores
#   - two for each dimension of the end scores
assoc.rcL.trans <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  tab <- prepareTable(model$data, FALSE)
  vars <- names(dimnames(tab))

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal") {
      rp <- prop.table(apply(tab, 1, sum, na.rm=TRUE))
      cp <- prop.table(apply(tab, 2, sum, na.rm=TRUE))
  }
  else if(weighting == "uniform") {
      rp <- rep(1/nr, nr)
      cp <- rep(1/nc, nc)
  }
  else {
      rp <- rep(1, nr)
      cp <- rep(1, nc)
  }


  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, sprintf("RCTrans.*inst = %s\\)\\.[RC][12]%s", nd+1, vars[1]))) > 0)
      nd <- nd + 1

  # One dimension, or none
  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("RCTrans.*\\)\\.R1\\Q%s\\E(\\Q%s\\E)$", vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- parameters(model)[pickCoef(model, sprintf("RCTrans.*\\)\\.C1\\Q%s\\E(\\Q%s\\E)$", vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]

      mu1 <- parameters(model)[pickCoef(model, sprintf("RCTrans.*\\)\\.R2\\Q%s\\E(\\Q%s\\E)",vars[1],
                                                 paste(rownames(tab), collapse="\\E|\\Q")))]
      nu1 <- parameters(model)[pickCoef(model, sprintf("RCTrans.*\\)\\.C2\\Q%s\\E(\\Q%s\\E)", vars[2],
                                                 paste(colnames(tab), collapse="\\E|\\Q")))]

      phi <- parameters(model)[pickCoef(model, sprintf("RCTrans.*\\)\\.\\Q%s\\E(\\Q%s\\E)", vars[3],
                                                paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      if(length(mu) == nr && length(nu) == nc
         && length(mu1) == nr && length(nu1) == nc
         && length(phi) == nl) {
          nd <- 1

          row <- matrix(mu, nr, 1)
          col <- matrix(nu, nc, 1)
          row1 <- matrix(mu1, nr, 1)
          col1 <- matrix(nu1, nc, 1)
          layer <- matrix(phi, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with transitional layer effect?")
      }
  }
  else {
      # Several dimensions: prepare matrices before filling them
      row <- matrix(NA, nr, nd)
      col <- matrix(NA, nc, nd)
      row1 <- matrix(NA, nr, nd)
      col1 <- matrix(NA, nc, nd)
      layer <- matrix(NA, nl, nd)

      for(i in 1:nd) {
          mu <- parameters(model)[pickCoef(model, sprintf("RCTrans.*inst = %s\\)\\.R1\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]
          nu <- parameters(model)[pickCoef(model, sprintf("RCTrans.*inst = %s\\)\\.C1\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[2],
                                                    paste(colnames(tab), collapse="\\E|\\Q")))]

          mu1 <- parameters(model)[pickCoef(model, sprintf("RCTrans.*inst = %s\\)\\.R2\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]
          nu1 <- parameters(model)[pickCoef(model, sprintf("RCTrans.*inst = %s\\)\\.C2\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[2],
                                                    paste(colnames(tab), collapse="\\E|\\Q")))]

          phi <- parameters(model)[pickCoef(model, sprintf("RCTrans.*inst = %s\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[3],
                                                    paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

          if(length(mu) == nr && length(nu) == nc
               && length(mu1) == nr && length(nu1) == nc
               && length(phi) == nl) {
              row[,i] <- mu
              col[,i] <- nu
              row1[,i] <- mu1
              col1[,i] <- nu1
              layer[,i] <- phi
          }
          else {
              stop("Invalid dimensions found. Are you sure this is a row-column association model with transitional layer effect?")
          }
      }
  }

  if(length(pickCoef(model, "Diag\\(") > 0)) {
      dg <- matrix(NA, nl, nr)
      dg[] <- parameters(model)[pickCoef(model, "Diag\\(")]
  }
  else {
      dg <- numeric(0)
  }


  # Layer coefficients are squared internally by RCTrans
  layer <- layer^2

  # Center
  row <- sweep(row, 2, colSums(sweep(row, 1, rp/sum(rp), "*")), "-")
  col <- sweep(col, 2, colSums(sweep(col, 1, cp/sum(cp), "*")), "-")
  row1 <- sweep(row1, 2, colSums(sweep(row1, 1, rp/sum(rp), "*")), "-")
  col1 <- sweep(col1, 2, colSums(sweep(col1, 1, cp/sum(cp), "*")), "-")

  # Compute layer scores and scale
  row <- array(row, dim=c(nrow(row), nd, nl)) +
             sweep(array(row1, dim=c(nrow(row1), nd, nl)), 3:2, layer, "*")
  col <- array(col, dim=c(nrow(col), nd, nl)) +
             sweep(array(col1, dim=c(nrow(col1), nd, nl)), 3:2, layer, "*")

  phir <- sqrt(margin.table(sweep(row^2, 1, rp, "*"), 2:3))
  phic <- sqrt(margin.table(sweep(col^2, 1, cp, "*"), 2:3))
  row <- sweep(row, 2:3, phir, "/")
  col <- sweep(col, 2:3, phic, "/")
  layer <- t(phir * phic)

  # Order dimensions according to phi on first layer category
  ord <- order(abs(layer[1,]), decreasing=TRUE)
  layer <- layer[,ord, drop=FALSE]
  row <- row[,ord,, drop=FALSE]
  col <- col[,ord,, drop=FALSE]


  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(row)[3]) {
      for(j in 1:nd) {
          if(row[1,j,i] < 0) {
              row[,j,i] <- -row[,j,i]
              col[,j,i] <- -col[,j,i]
          }
      }
  }

  # Scores can switch sides from one layer to another
  # Reverse axes for each layer state so that the scores are positively correlated to the first layer scores for rows
  for(i in seq_len(dim(row)[3] - 1)) {
      chg <- ifelse(diag(as.matrix(cor(row[,,i+1], row[,,1]))) >= 0, 1, -1)
      row[,,i+1] <- sweep(row[,,i+1, drop=FALSE], 2, chg, "*")
      col[,,i+1] <- sweep(col[,,i+1, drop=FALSE], 2, chg, "*")
  }

  ## Prepare objects
  colnames(layer) <- colnames(row) <- colnames(col) <- paste("Dim", 1:nd, sep="")
  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)
  rownames(layer) <- dimnames(row)[[3]] <- dimnames(col)[[3]] <- dimnames(tab)[[3]]

  if(length(dg) > 0) {
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")
      rownames(dg) <- dimnames(tab)[[3]]
  }

  row.weights <- apply(tab, c(1, 3), sum, na.rm=TRUE)
  col.weights <- apply(tab, c(2, 3), sum, na.rm=TRUE)

  obj <- list(phi = layer, row = row, col = col, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)

  class(obj) <- c("assoc.rcL.trans", "assoc.rcL", "assoc")
  obj
}


assoc.rcL.trans.symm <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  tab <- prepareTable(model$data, FALSE)
  vars <- names(dimnames(tab))

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  stopifnot(nr == nc)

  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(apply(tab, 1, sum, na.rm=TRUE) + apply(tab, 2, sum, na.rm=TRUE))
  else if(weighting == "uniform")
      rp <- rep(1/nr, nr)
  else
      rp <- rep(1, nr)


  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, sprintf("RCTransSymm.*inst = %s\\)[RC][12]\\.\\Q%s\\E\\|\\Q%s\\E",
                                       nd+1, vars[1], vars[2]))) > 0)
      nd <- nd + 1

  # One dimension, or none
  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("RCTransSymm.*\\)[RC]1\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]

      mu1 <- parameters(model)[pickCoef(model, sprintf("RCTransSymm.*\\)[RC]2\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                 vars[1], vars[2],
                                                 paste(rownames(tab), collapse="\\E|\\Q")))]

      phi <- parameters(model)[pickCoef(model, sprintf("RCTransSymm.*\\)\\.\\Q%s\\E(\\Q%s\\E)$", vars[3],
                                                paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      if(length(mu) == nr && length(mu1) == nr && length(phi) == nl) {
          nd <- 1

          sc <- matrix(mu, nr, 1)
          sc1 <- matrix(mu1, nr, 1)
          layer <- matrix(phi, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with transitional layer effect?")
      }
  }
  else {
      # Several dimensions: prepare matrices before filling them
      sc <- matrix(NA, nr, nd)
      sc1 <- matrix(NA, nr, nd)
      layer <- matrix(NA, nl, nd)

      for(i in 1:nd) {
          mu <- parameters(model)[pickCoef(model, sprintf("RCTransSymm.*inst = %i\\)[RC]1\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1], vars[2],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]

          mu1 <- parameters(model)[pickCoef(model, sprintf("RCTransSymm.*inst = %i\\)[RC]2\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1], vars[2],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]

          phi <- parameters(model)[pickCoef(model, sprintf("RCTransSymm.*inst = %i\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[3],
                                                    paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

          if(length(mu) == nr && length(mu1) == nr && length(phi) == nl) {
              sc[,i] <- mu
              sc1[,i] <- mu1
              layer[,i] <- phi
          }
          else {
              stop("Invalid dimensions found. Are you sure this is a row-column association model with transitional layer effect?")
          }
      }
  }

  if(length(pickCoef(model, "Diag\\(") > 0)) {
      dg <- matrix(NA, nl, nr)
      dg[] <- parameters(model)[pickCoef(model, "Diag\\(")]
  }
  else {
      dg <- numeric(0)
  }


  # Layer coefficients are squared internally by RCTrans
  layer <- layer^2

  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")
  sc1 <- sweep(sc1, 2, colSums(sweep(sc1, 1, p/sum(p), "*")), "-")


   # Compute layer scores and scale
   sc <- array(sc, dim=c(nrow(sc), nd, nl)) +
             sweep(array(sc1, dim=c(nrow(sc1), nd, nl)), 3:2, layer, "*")

  phi <- sqrt(margin.table(sweep(sc^2, 1, p, "*"), 2:3))
  sc <- sweep(sc, 2:3, phi, "/")
  layer <- t(phi)

  # Order dimensions according to phi on first layer category
  ord <- order(abs(layer[1,]), decreasing=TRUE)
  layer <- layer[,ord, drop=FALSE]
  sc <- sc[,ord,, drop=FALSE]


  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(sc)[3]) {
      for(j in 1:nd) {
          if(sc[1,j,i] < 0)
              sc[,j,i] <- -sc[,j,i]
      }
  }

  # Scores can switch sides from one layer to another
  # Reverse axes for each layer state so that the scores are positively correlated to the first layer scores for rows
  for(i in seq_len(dim(sc)[3] - 1)) {
      chg <- ifelse(diag(as.matrix(cor(sc[,,i+1], sc[,,1]))) >= 0, 1, -1)
      sc[,,i+1] <- sweep(sc[,,i+1, drop=FALSE], 2, chg, "*")
  }

  ## Prepare objects
  colnames(layer) <- colnames(sc) <- paste("Dim", 1:nd, sep="")
  rownames(sc) <- rownames(tab)
  rownames(layer) <- dimnames(sc)[[3]] <- dimnames(tab)[[3]]

  if(length(dg) > 0) {
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")
      rownames(dg) <- dimnames(tab)[[3]]
  }

  row.weights <- apply(tab, c(1, 3), sum, na.rm=TRUE)
  col.weights <- apply(tab, c(2, 3), sum, na.rm=TRUE)

  obj <- list(phi = layer, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights)

  class(obj) <- c("assoc.rcL.trans", "assoc.rcL", "assoc.symm", "assoc")
  obj
}
