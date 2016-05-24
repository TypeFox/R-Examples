## RC(M)-L model (Wong, 2010)

rcL <- function(tab, nd=1, layer.effect=c("homogeneous.scores", "heterogeneous", "none"),
                symmetric=FALSE, diagonal=c("none", "heterogeneous", "homogeneous"),
                weighting=c("marginal", "uniform", "none"), se=c("none", "jackknife", "bootstrap"),
                nreplicates=100, ncpus=getOption("boot.ncpus"),
                family=poisson, weights=NULL, start=NULL, etastart=NULL, tolerance=1e-8, iterMax=5000,
                trace=FALSE, verbose=TRUE, ...) {
  layer.effect <- match.arg(layer.effect)
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


  if(diagonal == "heterogeneous")
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])
  else if(diagonal == "homogeneous")
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""

  f1 <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s",
                vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3])

  eliminate <- eval(parse(text=sprintf("quote(%s:%s)", vars[1], vars[3])))

  base <- NULL

  nastart <- length(start) == 1 && is.na(start)


  if(nastart) {
      cat("Running base model to find starting values...\n")

      args <- list(formula=as.formula(paste(f1, diagstr)), data=tab,
                   family=family, weights=weights, eliminate=eliminate,
                   tolerance=1e-6, iterMax=iterMax)

      base <- do.call("gnm", c(args, list(...)))

      res <- if(symmetric) residEVDL(base, nd, layer.effect)
             else residSVDL(base, nd, layer.effect)

      # Using NA for all linear parameters usually gives better results than linear parameter values
      if(layer.effect == "homogeneous.scores")
          start <- c(rep(NA, length(parameters(base))), rbind(matrix(NA, dim(tab)[3], ncol(res)), res))
      else
          start <- c(rep(NA, length(parameters(base))), res)
  }

  if(symmetric) {
      if(layer.effect == "homogeneous.scores") {
          f2 <- ""

          for(i in 1:nd)
              f2 <- paste(f2, sprintf("+ Mult(%s, MultHomog(%s, %s), inst = %i)",
                                      vars[3], vars[1], vars[2], i))
      }
      else if(layer.effect == "heterogeneous") {
          stop("Symmetric association with heterogeneous layer effect is currently not supported")

          f2 <- ""

          for(i in 1:nd)
              f2 <- paste(f2, sprintf("+ MultHomog(%s:%s, %s:%s, inst = %i)", 
                                      vars[3], vars[1], vars[3], vars[2], i))
      }
      else {
          f2 <- sprintf("+ instances(MultHomog(%s, %s), %i)", vars[1], vars[2], nd)
      }
  }
  else {
      if(layer.effect == "homogeneous.scores") {
          f2 <- sprintf("+ instances(Mult(%s, %s, %s), %i)",
                        vars[3], vars[1], vars[2], nd)
      }
      else if(layer.effect == "heterogeneous") {
          f2 <- sprintf("+ instances(Mult(%s:%s, %s:%s), %i)",
                        vars[3], vars[1], vars[3], vars[2], nd)
      }
      else {
          f2 <- sprintf("+ instances(Mult(%s, %s), %i)",
                        vars[1], vars[2], nd)
      }
  }

  if(!is.null(base) && is.null(etastart))
      etastart <- as.numeric(predict(base))

  if(!is.null(base))
      cat("Running real model...\n")

  args <- list(formula=as.formula(paste(f1, diagstr, f2)), data=tab,
               family=family, weights=weights, start=start, etastart=etastart,
               eliminate=eliminate,
               tolerance=tolerance, iterMax=iterMax, verbose=verbose, trace=trace)

  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  newclasses <- if(symmetric) c("rcL.symm", "rcL", "assocmod") else c("rcL", "assocmod")
  class(model) <- c(newclasses, class(model))

  model$call.gnm <- model$call
  model$call <- match.call()

  if(layer.effect == "none") {
      if(symmetric) {
          model$assoc <- assoc.rc.symm(model, weighting=weighting)
          assoc1 <- assoc.rc.symm
      }
      else {
          model$assoc <- assoc.rc(model, weighting=weighting)
          assoc1 <- assoc.rc
      }
  }
  else {
      model$assoc <- assoc(model, weighting=weighting)
      assoc1 <- getS3method("assoc", class(model))
  }

  class(model$assoc) <- if(symmetric) c("assoc.rcL", "assoc.symm", "assoc")
                        else c("assoc.rcL", "assoc")

  model$call.gnm <- model$call
  model$call <- match.call()


  if(se %in% c("jackknife", "bootstrap")) {
      jb <- jackboot(se, ncpus, nreplicates, tab, model, assoc1, NULL,
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


assoc.rcL <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
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
  while(length(pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*inst = %i\\)", vars[3], nd + 1))) > 0)
      nd <- nd + 1

  homogeneous <- TRUE

  # Only one dimension, or none
  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*\\).*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                vars[3], vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*\\).*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                vars[3], vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]
      phi <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                 vars[3],
                                                 paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      # Homogeneous scores for rows and/or columns
      if(length(phi) == nl &&
         (length(mu) == nr || length(mu) == nr * nl) &&
         (length(nu) == nc || length(nu) == nc * nl)) {
          nd <- 1

          if(length(mu) == nr) {
              row <- array(mu, dim=c(nr, 1, 1))
          }
          else {
              row <- array(mu, dim=c(nr, nl, 1))
              homogeneous <- FALSE
          }

          if(length(nu) == nc) {
              col <- array(nu, dim=c(nc, 1, 1))
          }
          else {
              col <- array(nu, dim=c(nc, nl, 1))
              homogeneous <- FALSE
          }

          layer <- matrix(phi, nl, 1)
      }
      # Fully heterogeneous scores
      else if(length(phi) == 0 &&
              length(mu) == nr * nl &&
              length(nu) == nc * nl) {
          homogeneous <- FALSE
          row <- array(mu, dim=c(nr, 1, nl))
          col <- array(nu, dim=c(nc, 1, nl))
          layer <- matrix(1, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with layer effect?")
      }
  }
  else {
      # Several dimensions: prepare arrays before filling them
      row <- array(NA, dim=c(nr, nd, nl))
      col <- array(NA, dim=c(nc, nd, nl))
      layer <- matrix(NA, nl, nd)

      for(i in 1:nd) {
          mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*inst = %i.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                    vars[3], i, vars[1],
                                                    paste(rownames(tab), collapse="\\E|\\Q")))]
          nu <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*inst = %i.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                    vars[3], i, vars[2],
                                                    paste(colnames(tab), collapse="\\E|\\Q")))]
          phi <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*inst = %i.*\\.\\Q%s\\E(\\Q%s\\E)$",
                                                     i, vars[3],
                                                     paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

          # Homogeneous scores for rows and/or columns
          if(length(phi) == nl &&
             (length(mu) == nr || length(mu) == nr * nl) &&
             (length(nu) == nc || length(nu) == nc * nl)) {
              if(length(mu) == nr)
                  row[,i,1] <- mu
              else {
                  row[,i,] <- t(matrix(mu, nl, nr))
                  homogeneous <- FALSE
              }

              if(length(nu) == nc)
                  col[,i,1] <- nu
              else {
                  col[,i,] <- t(matrix(nu, nl, nc))
                  homogeneous <- FALSE
              }

              layer[,i] <- phi
          }
          # Fully heterogeneous scores
          else if(length(phi) == 0 &&
                  length(mu) == nr * nl &&
                  length(nu) == nc * nl) {
              homogeneous <- FALSE
              row[,i,] <- t(matrix(mu, nl, nr))
              col[,i,] <- t(matrix(nu, nl, nc))
              layer[,i] <- 1
          }
          else {
              stop("Invalid dimensions found. Are you sure this is a row-column association model with layer effect?")
          }
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nr)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], nl, nr)
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nr)
  else
      dg <- numeric(0)


  # Center
  row <- sweep(row, 2:3, margin.table(sweep(row, 1, rp/sum(rp), "*"), 2:3), "-")
  col <- sweep(col, 2:3, margin.table(sweep(col, 1, cp/sum(cp), "*"), 2:3), "-")

  if(homogeneous) {
      # Scale
      phi.row <- sqrt(margin.table(sweep(row[,,1, drop=FALSE]^2, 1, rp, "*"), 2))
      phi.col <- sqrt(margin.table(sweep(col[,,1, drop=FALSE]^2, 1, cp, "*"), 2))
      row <- sweep(row[,,1, drop=FALSE], 2, phi.row, "/")
      col <- sweep(col[,,1, drop=FALSE], 2, phi.col, "/")
      layer <- sweep(layer, 2, phi.row * phi.col, "*")

      # Conventionally order dimensions according to phi on first layer category
      ord <- order(abs(colSums(sweep(layer, 1, apply(tab, 3, sum, na.rm=TRUE), "*"))), decreasing=TRUE)
      layer <- layer[,ord, drop=FALSE]
      row <- row[,ord,, drop=FALSE]
      col <- col[,ord,, drop=FALSE]
  }
  else {
      for(l in 1:nl) {
          # Technique proposed in Goodman (1991), Appendix 4
          lambda <- matrix(0, nr, nc)
          for(i in 1:nd) {
              lambda <- lambda + layer[l,i] * row[,i,l] %o% col[,i,l]
          }
          lambda0 <- lambda * sqrt(rp %o% cp) # Eq. A.4.3
          sv <- svd(lambda0)
          row[,,l] <- diag(1/sqrt(rp)) %*% sv$u[,1:nd] # Eq. A.4.7
          col[,,l] <- diag(1/sqrt(cp)) %*% sv$v[,1:nd] # Eq. A.4.7
          layer[l,] <- sv$d[1:nd]
      }
  }

  # By convention, keep the weighted average of layer coefficients positive
  if(homogeneous) {
      for(i in 1:nd) {
          if(sum(layer[,i] * apply(tab, 3, sum, na.rm=TRUE)) < 0) {
              layer[,i] <- -layer[,i]
              row[,i,] <- -row[,i,]
          }
      }
  }

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

  ## Prepare objects
  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)
  colnames(row) <- colnames(col) <- colnames(layer) <- paste("Dim", 1:nd, sep="")
  rownames(layer) <- dimnames(tab)[[3]]

  if(!homogeneous)
      dimnames(row)[[3]] <- dimnames(col)[[3]] <- dimnames(tab)[[3]]
  else
      dimnames(row)[[3]] <- dimnames(col)[[3]] <- "All levels"

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

  row.weights <- apply(tab, c(1, 3), sum, na.rm=TRUE)
  col.weights <- apply(tab, c(2, 3), sum, na.rm=TRUE)

  obj <- list(phi = layer, row = row, col = col, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)

  class(obj) <- c("assoc.rcL", "assoc")
  obj
}

assoc.rcL.symm <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  tab <- prepareTable(model$data, FALSE)
  vars <- names(dimnames(tab))

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  # Weight with marginal frequencies, cf. Clogg & Shihadeh (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(apply(tab, 1, sum, na.rm=TRUE) + apply(tab, 2, sum, na.rm=TRUE))
  else if(weighting == "uniform")
      p <- rep(1/nr, nr)
  else
      p <- rep(1, nr)


  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, paste("Mult\\(.*MultHomog\\(.*inst =", nd + 1))) > 0)
      nd <- nd + 1

  homogeneous <- TRUE

  # Only one dimension, or none
  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*MultHomog\\(.*\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]
      phi <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*MultHomog\\(.*\\.\\Q%s\\E(\\Q%s\\E)$", vars[3],
                                                 paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      # Homogeneous scores
      if(length(phi) == nl &&
         (length(mu) == nr || length(mu) == nr * nl)) {
          nd <- 1

          if(length(mu) == nr) {
              sc <- array(mu, dim=c(nr, 1, 1))
          }
          else {
              sc <- array(mu, dim=c(nr, nl, 1))
              homogeneous <- FALSE
          }

          layer <- matrix(phi, nl, 1)
      }
      # Fully heterogeneous scores
      else if(length(phi) == 0 &&
              length(mu) == nr * nl) {
          homogeneous <- FALSE
          sc <- array(mu, dim=c(nr, 1, nl))
          layer <- matrix(1, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a symmetric row-column association model with layer effect?")
      }
  }
  else {
      # Several dimensions: prepare arrays before filling them
      sc <- array(NA, dim=c(nr, nd, nl))
      layer <- matrix(NA, nl, nd)

      for(i in 1:nd) {
          mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*inst = %i.*MultHomog\\(.*\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                    i, vars[1], vars[2],
                                                    paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]
          phi <- parameters(model)[pickCoef(model, sprintf("Mult\\(.*MultHomog\\(.*inst = %i.*\\.\\Q%s\\E(\\Q%s\\E)$",
                                                     i, vars[3],
                                                     paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

          # Homogeneous scores
          if(length(phi) == nl &&
             (length(mu) == nr || length(mu) == nr * nl)) {
              if(length(mu) == nr)
                  sc[,i,1] <- mu
              else {
                  sc[,i,] <- t(matrix(mu, nl, nr))
                  homogeneous <- FALSE
              }

              layer[,i] <- phi
          }
          # Fully heterogeneous scores
          else if(length(phi) == 0 &&
                  length(mu) == nr * nl) {
              homogeneous <- FALSE
              sc[,i,] <- t(matrix(mu, nl, nr))
              layer[,i] <- 1
          }
          else {
              stop("Invalid dimensions found. Are you sure this is a symmetric row-column association model with layer effect?")
          }
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nr)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], nl, nr)
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nr)
  else
      dg <- numeric(0)

  # Center
  sc <- sweep(sc, 2:3, margin.table(sweep(sc, 1, p/sum(p), "*"), 2:3), "-")

  if(homogeneous) {
      # Scale
      phi <- sqrt(margin.table(sweep(sc[,,1, drop=FALSE]^2, 1, p/sum(p), "*"), 2))
      sc <- sweep(sc[,,1, drop=FALSE], 2, phi, "/")
      layer <- sweep(layer, 2, phi, "*")

      # Conventionally order dimensions according to weighted average of phi
      ord <- order(abs(colSums(sweep(layer, 1, apply(tab, 3, sum, na.rm=TRUE), "*"))), decreasing=TRUE)
      layer <- layer[,ord, drop=FALSE]
      sc <- sc[,ord,, drop=FALSE]
  }
  else {
      for(l in 1:nl) {
          # Technique proposed in Goodman (1991), Appendix 4, but with eigenvalues decomposition
          lambda <- matrix(0, nr, nc)
          for(i in 1:nd)
              lambda <- lambda + sc[,i,l] %o% sc[,i,l]
          lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
          eigen <- eigen(lambda0, symmetric=TRUE)
          sc[,,l] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
          layer[l,] <- eigen$values[1:nd]
      }
  }

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(sc)[3]) {
      for(j in 1:nd) {
          if(sc[1,j,i] < 0)
              sc[,j,i] <- -sc[,j,i]
      }
  }

  ## Prepare objects
  rownames(sc) <- rownames(tab)
  colnames(sc) <- colnames(layer) <- paste("Dim", 1:nd, sep="")
  rownames(layer) <- dimnames(tab)[[3]]

  if(!homogeneous)
      dimnames(sc)[[3]] <- dimnames(tab)[[3]]
  else
      dimnames(sc)[[3]] <- "All levels"


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

  row.weights <- apply(tab, c(1, 3), sum, na.rm=TRUE)
  col.weights <- apply(tab, c(2, 3), sum, na.rm=TRUE)

  obj <- list(phi = layer, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = row.weights, col.weights = col.weights,
              vars = vars)

  class(obj) <- c("assoc.rcL", "assoc.symm", "assoc")
  obj
}
