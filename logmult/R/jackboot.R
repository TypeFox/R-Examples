# Run jackknife or bootstrap replicates of the model
jackboot <- function(se, ncpus, nreplicates, tab, model, assoc1, assoc2,
                     weighting, rowsup=NULL, colsup=NULL, family, weights,
                     verbose, trace, start, etastart,
                     formula=NULL, design=NULL, Ntotal=NULL, exclude=c(NA, NaN), ...) {
  cat("Computing", se, "standard errors...\n")

  if(is.null(ncpus))
      ncpus <- if(requireNamespace("parallel")) min(parallel::detectCores(), 5)
               else 1

  if(ncpus > 1 && requireNamespace("parallel")) {
      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus), outfile="", methods=FALSE)
      on.exit(parallel::stopCluster(cl))

      libpaths <- .libPaths()
      parallel::clusterExport(cl, "libpaths", env=environment())
      parallel::clusterEvalQ(cl, library("logmult", lib.loc=libpaths,
                                         warn.conflicts=FALSE, quietly=TRUE))

      # Printing output from all nodes at the same time would be a mess, only print "."
      trace <- FALSE
  }
  else {
      cl <- NULL
      ncpus <- 1
  }

  ass1 <- assoc1(model, weighting=weighting, rowsup=rowsup, colsup=colsup)

  nl1 <- nrow(ass1$phi)
  nd1 <- ncol(ass1$phi)
  nr1 <- nrow(ass1$row)
  nlr1 <- dim(ass1$row)[3]
  nc1 <- nrow(ass1$col)
  nlc1 <- dim(ass1$col)[3]

  len1 <- nl1 * nd1 + nlr1 * nd1 * nr1 + nlc1 * nd1 * nc1
  end1 <- len1 + nl1 * nd1 * nr1 + nl1 * nd1 * nc1

  if(is.null(assoc2)) {
      int <- 1:len1
  }
  else {
      ass2 <- assoc2(model, weighting=weighting, rowsup=rowsup, colsup=colsup)

      nl2 <- nrow(ass2$phi)
      nd2 <- ncol(ass2$phi)
      nr2 <- nrow(ass2$row)
      nlr2 <- dim(ass2$row)[3]
      nc2 <- nrow(ass2$col)
      nlc2 <- dim(ass2$col)[3]
    
      len2 <- nl2 * nd2 + nlr2 * nd2 * nr2 + nlc2 * nd2 * nc2
      end2 <- end1 + len2 + nl2 * nd2 * nr2 + nl2 * nd2 * nc2
      int <- c(1:len1, end1 + 1:len2)
  }

  nacount <- 0

  # Add supplementary rows and columns to 'tab' so that they are included in the sample
  if(!is.null(rowsup) && !is.null(colsup)) {
      # Block matrix that is symmetric if rowsup and colsup are transposes of each other
      tab2 <- as.table(cbind(rbind(tab, rowsup),
                                  rbind(colsup, matrix(NA, nrow(rowsup), ncol(colsup)))))
      dimnames(tab2) <- list(c(rownames(tab), rownames(rowsup)),
                             c(colnames(tab), colnames(colsup)))
      names(dimnames(tab2)) <- names(dimnames(tab))
      tab <- tab2
  }
  else if(!is.null(rowsup)) {
      tab2 <- as.table(rbind(tab, rowsup))
      names(dimnames(tab2)) <- names(dimnames(tab))
      tab <- tab2
  }
  else if(!is.null(colsup)) {
      tab2 <- as.table(cbind(tab, colsup))
      names(dimnames(tab2)) <- names(dimnames(tab))
      tab <- tab2
  }


  if(se == "jackknife") {
      jack <- jackknife((1:length(tab))[!is.na(tab)], jackknife.assoc,
                        # Load-balance only when model is likely to take enough some time to fit
                        w=tab[!is.na(tab)], cl=cl, load.balancing=length(coef(model)) > 200,
                        model=model, tab=tab, assoc1=assoc1, assoc2=assoc2,
                        weighting=weighting, rowsup=rowsup, colsup=colsup,
                        family=family, weights=weights,
                        verbose=verbose, trace=trace,
                        start=start, etastart=etastart, ...)

      rowsna <- rowSums(is.na(jack$values))
      nacount <- sum(rowsna > 0)

      if(!any(rowsna == 0)) {
          warning("All model replicates failed. Cannot compute standard errors.")
          return(list())
      }

      jack.results <- list(bias=jack$bias, values=jack$values)

      # Build the matrix only for phi on all layers, and normalized scores for the layers they vary on
      tot <- sum(tab, na.rm=TRUE)
      get.covmat <- function(start, len) {
          int <- seq.int(start, length.out=len)
          (tot - 1)/tot * t(jack$dev[, int]) %*% sweep(jack$dev[, int], 1, tab[!is.na(tab)], "*")
      }

      boot.results <- numeric(0)
      svyrep.results <- numeric(0)
  }
  else if(se == "bootstrap") {
      boot.results <- boot::boot(1:sum(tab, na.rm=TRUE), boot.assoc, R=nreplicates,
                                 parallel="snow", cl=cl, ncpus=ncpus,
                                 args=list(model=model, tab=tab, assoc1=assoc1, assoc2=assoc2,
                                           weighting=weighting, rowsup=rowsup, colsup=colsup,
                                           family=family, weights=weights,
                                           verbose=verbose, trace=trace,
                                           start=start, etastart=etastart, ...))

      rowsna <- rowSums(is.na(boot.results$t))
      nacount <- sum(rowsna > 0)

      if(!any(rowsna == 0)) {
          warning("All model replicates failed. Cannot compute standard errors.")
          return(list())
      }


      get.covmat <- function(start, len) {
          int <- seq.int(start, length.out=len)
          cov(boot.results$t[, int], use="na.or.complete")
      }

      jack.results <- numeric(0)
      svyrep.results <- numeric(0)
  }
  else { # Survey replicate weights
      # Not very efficient, but allows checking for errors before starting the cluster
      thetahat <- as.numeric(svyrep.assoc(survey::svytable(formula, design, Ntotal=Ntotal, exclude=exclude),
                                          model=model, assoc1=assoc1, assoc2=assoc2,
                                          weighting=weighting, rowsup=rowsup, colsup=colsup,
                                          family=family, weights=weights,
                                          verbose=verbose, trace=trace,
                                          start=start, etastart=etastart, ...))

      svyrep.results <- svyrep(formula, design, svyrep.assoc, Ntotal=Ntotal, exclude=exclude,
                               # Load-balance only when model is likely to take enough some time to fit
                               cl=cl, load.balancing=length(coef(model)) > 200,
                               model=model, assoc1=assoc1, assoc2=assoc2,
                               weighting=weighting, rowsup=rowsup, colsup=colsup,
                               family=family, weights=weights,
                               verbose=verbose, trace=trace,
                               start=start, etastart=etastart, ...)

      rowsna <- rowSums(is.na(svyrep.results))
      nacount <- sum(rowsna > 0)

      if(!any(rowsna == 0)) {
          warning("All model replicates failed. Cannot compute standard errors.")
          return(list())
      }

      get.covmat <- function(start, len) {
          int <- seq.int(start, length.out=len)
          survey::svrVar(svyrep.results, design$scale, design$rscales, mse = design$mse, coef = thetahat)[int, int]
      }

      boot.results <- numeric(0)
      jack.results <- numeric(0)
  }

  # "." progress indicators need this
  if(verbose && !trace)
      cat("\n")

  if(nacount > 0)
      warning(ngettext(nacount,
                       "One model replicate failed and its results will be skipped. Standard errors will not be completely accurate. Consider raising the value of iterMax.",
                       sprintf("%i model replicates failed and their results will be skipped. Standard errors will not be completely accurate. Consider raising the value of iterMax.", nacount)))

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  submat.names <- function(ass) {
      nd <- ncol(ass$phi)

      outer(c(paste(vars[1], rownames(ass$row), sep=""),
              paste(vars[2], rownames(ass$col), sep="")),
            paste("Dim", 1:nd, "*", sep=""),
            paste, sep=":")
  }

  mat.names <- function(ass) {
      if(length(rownames(ass$phi)) > 0)
          lnames <- paste(":", vars[3], rownames(ass$phi), sep="")
      else
          lnames <- ""

      if(dim(ass$row)[3] > 1)
          lrcnames <- lnames
      else
          lrcnames <- ""

      nd <- ncol(ass$phi)

      scnames <- outer(c(paste(vars[1], rownames(ass$row), sep=""),
                         paste(vars[2], rownames(ass$col), sep="")),
                       outer(paste("Dim", 1:nd, sep=""),
                             lrcnames,
                             paste, sep=""),
                       paste, sep=":")

      c(outer(lnames, paste("Dim", 1:nd, sep=""), paste, sep=""), scnames)
  }

  values.names <- function(ass) {
      nd <- ncol(ass$phi)

      if(length(rownames(ass$phi)) > 0)
          lnames <- paste(":", vars[3], rownames(ass$phi), "*", sep="")
      else
          lnames <- "*"

      c(mat.names(ass),
        outer(c(paste(vars[1], rownames(ass$row), sep=""),
                paste(vars[2], rownames(ass$col), sep="")),
                outer(paste("Dim", 1:nd, sep=""),
                      lnames,
                      paste, sep=""),
                paste, sep=":"))
  }

  # Main matrix only for phi on all layers, and normalized scores for the layers they vary on
  tot <- sum(tab, na.rm=TRUE)
  covmat1 <- get.covmat(1, len1)
  rownames(covmat1) <- colnames(covmat1) <- mat.names(ass1)

  # Sub matrices with adjusted scores for all layers separately
  len <- nd1 * (nr1 + nc1)
  adj.covmats1 <- vapply(seq.int(len1 + 1, end1, by=len), get.covmat, len=len,
                         FUN.VALUE=matrix(0, len, len))

  rownames(adj.covmats1) <- colnames(adj.covmats1) <- submat.names(ass1)
  dimnames(adj.covmats1)[[3]] <- rownames(ass1$phi)

  if(!is.null(assoc2)) {
      # Main matrix only for phi on all layers, and normalized scores for the layers they vary on
      covmat2 <- get.covmat(end1 + 1, len2)
      rownames(covmat2) <- colnames(covmat2) <- mat.names(ass2)

      # Sub matrices with adjusted scores for all layers separately
      len <- nd2 * (nr2 + nc2)
      adj.covmats2 <- vapply(seq.int(end1 + len2 + 1, end2, by=len), get.covmat, len=len,
                             FUN.VALUE=matrix(0, len, len))
      rownames(adj.covmats2) <- colnames(adj.covmats2) <- submat.names(ass2)
      dimnames(adj.covmats2)[[3]] <- rownames(ass2$phi)
  }

  int <- seq.int(1, len1 + nl1 * nd1 * (nr1 + nc1))

  if(length(boot.results) > 0) {
      boot.results1 <- boot.results
      boot.results1$t0 <- boot.results$t0[int]
      boot.results1$t <- boot.results$t[,int]
      names(boot.results1$t0) <- colnames(boot.results1$t) <- values.names(ass1)
  }
  else {
      boot.results1 <- numeric(0)
  }

  if(length(jack.results) > 0) {
      jack.results1 <- jack.results
      jack.results1$bias <- jack.results$bias[int]
      jack.results1$values <- jack.results$values[,int]
      names(jack.results1$bias) <- colnames(jack.results1$values) <- values.names(ass1)
  }
  else {
      jack.results1 <- numeric(0)
  }

  if(length(svyrep.results) > 0) {
      svyrep.results1 <- svyrep.results[,int]
  }
  else {
      svyrep.results1 <- numeric(0)
  }

  if(!is.null(assoc2)) {
      if(length(boot.results) > 0) {
          boot.results2 <- boot.results
          boot.results2$t0 <- boot.results$t0[-int]
          boot.results2$t <- boot.results$t[,-int]
          names(boot.results2$t0) <- colnames(boot.results2$t) <- values.names(ass2)
      }
      else {
          boot.results2 <- numeric(0)
      }

      if(length(jack.results) > 0) {
          jack.results2 <- jack.results
          jack.results2$bias <- jack.results$bias[-int]
          jack.results2$values <- jack.results$values[,-int]
          names(jack.results2$bias) <- colnames(jack.results2$values) <- values.names(ass2)
      }
      else {
          jack.results2 <- numeric(0)
      }

      if(length(svyrep.results) > 0) {
          svyrep.results2 <- svyrep.results[,-int]
      }
      else {
          svyrep.results2 <- numeric(0)
      }
  }
  else {
      covmat2 <- numeric(0)
      boot.results2 <- numeric(0)
      jack.results2 <- numeric(0)
      svyrep.results2 <- numeric(0)
  }

  if(is.null(assoc2))
      list(covmat=covmat1, adj.covmats=adj.covmats1,
           boot.results=boot.results1, jack.results=jack.results1, svyrep.results=svyrep.results1)
  else
      list(covmat1=covmat1, adj.covmats1=adj.covmats1,
           boot.results1=boot.results1, jack.results1=jack.results1, svyrep.results1=svyrep.results1,
           covmat2=covmat2, adj.covmats2=adj.covmats2,
           boot.results2=boot.results2, jack.results2=jack.results2, svyrep.results2=svyrep.results2)
}

# Additional arguments are needed so that update() finds them even when using parLapply
jackknife.assoc <- function(x, tab, model, rowsup, colsup, ...) {
  if(sum(tab[-x], na.rm=TRUE) > 0) {
      mat <- tab
      mat[] <- -1
      mat[x] <- 0
      tab <- tab + mat
  }

  if(!is.null(rowsup) || !is.null(colsup))
      # Remove supplementary rows and columns and recreate them separately
      tab <- as.table(tab[seq(nrow(model$data)), seq(ncol(model$data))])
  else
      tab <- as.table(tab)

  if(!is.null(rowsup))
      rowsup <- tab[-seq(nrow(model$data)), seq(ncol(model$data))]

  if(!is.null(colsup))
      colsup <- tab[seq(nrow(model$data)), -seq(ncol(model$data))]

  replicate.assoc(model, tab, rowsup, colsup, ...)
}

boot.assoc <- function(data, indices, args) {
  tab <- args$tab

  # Create a table from the indices - one index identifies an observation in the original table,
  # following the cumulative sum, from 1 to sum(tab)
  tab[!is.na(tab) & tab > 0] <- tapply(tabulate(indices, nbins=sum(tab, na.rm=TRUE)),
                                      rep.int(1:length(tab), ifelse(is.na(tab), 0, tab)),
                                      sum)

  # Basic sanity check
  stopifnot(sum(tab, na.rm=TRUE) == sum(args$model$data, args$rowsup, args$colsup, na.rm=TRUE))

  if(!is.null(args$rowsup) || !is.null(args$colsup))
      # Remove supplementary rows and columns and recreate them separately
      args$tab <- as.table(tab[seq(nrow(args$model$data)), seq(ncol(args$model$data))])
  else
      args$tab <- as.table(tab)

  if(!is.null(args$rowsup))
      args$rowsup <- tab[-seq(nrow(args$model$data)), seq(ncol(args$model$data))]

  if(!is.null(args$colsup))
      args$colsup <- tab[seq(nrow(args$model$data)), -seq(ncol(args$model$data))]

  # We need to pass all arguments through "args" to prevent them
  # from being catched by boot(), especially "weights"
  do.call(replicate.assoc, args)
}

# Additional arguments are needed so that update() finds them even when using parLapply
svyrep.assoc <- function(tab, model, rowsup, colsup, ...) {
  if(!is.null(rowsup) || !is.null(colsup))
      # Remove supplementary rows and columns and recreate them separately
      tab <- as.table(tab[seq(nrow(model$data)), seq(ncol(model$data))])
  else
      tab <- as.table(tab)

  if(!is.null(rowsup))
      rowsup <- tab[-seq(nrow(model$data)), seq(ncol(model$data))]

  if(!is.null(colsup))
      colsup <- tab[seq(nrow(model$data)), -seq(ncol(model$data))]

  replicate.assoc(model, tab, rowsup, colsup, ...)
}

# Replicate model with new data, and combine assoc components into a vector
replicate.assoc <- function(model.orig, tab, assoc1, assoc2,
                            weighting, rowsup, colsup,
                            verbose, trace, start, etastart, ...) {

  if(!all(dim(tab) == dim(model.orig$data))) {
      cat("Dimensions of the table for replicate do not match that of the original table: ignoring the results of this replicate. Data was:\n")

      print(tab)

      # This NA value is skipped when computing variance-covariance matrix
      return(NA)
  }
  else if(any(sapply(seq_along(dim(tab)), function(i)
                     any((apply(tab, i, sum, na.rm=TRUE) == 0) != (apply(model.orig$data, i, sum, na.rm=TRUE) == 0))))) {
      cat("Some rows, columns or layers of the table for replicate have zero counts that do not appear in the original table: ignoring the results of this replicate. Data was:\n")

      print(tab)

      # This NA value is skipped when computing variance-covariance matrix
      return(NA)
  }

  # Models can generate an error if they fail repeatedly
  # Remove warnings because we handle them below
  model <- tryCatch(suppressWarnings(update(model.orig, tab=tab,
                                            rowsup=rowsup, colsup=colsup,
                                            start=parameters(model.orig),
                                            etastart=as.numeric(predict(model.orig)),
                                            verbose=trace, trace=trace, se="none")),
                    error=function(e) NULL)

  if(is.null(model) || !model$converged) {
      if(is.null(model)) {
          cat("Model replicate failed.\nData was:\n")
      }
      else {
          cat("Model replicate did not converge.\nData was:\n")
      }

      print(tab)
  }

  if(!is.null(start) && (is.null(model) || !model$converged)) {
      cat("Trying again with starting values from base model...\n")

      model <- tryCatch(suppressWarnings(update(model.orig, tab=tab,
                                                rowsup=rowsup, colsup=colsup,
                                                start=start, etastart=etastart,
                                                verbose=verbose, trace=verbose, se="none")),
                        error=function(e) NULL)
      if(is.null(model) || !model$converged)
          cat("Model failed again. ")
  }

  if(is.null(model) || !model$converged) {
      cat("Trying one last time with default starting values...\n")
      model <- tryCatch(suppressWarnings(update(model.orig, tab=tab,
                                                rowsup=rowsup, colsup=colsup,
                                               start=NA, etastart=NULL,
                                                verbose=verbose, trace=verbose, se="none")),
                        error=function(e) NULL)
  }

  if(is.null(model) || !model$converged) {
      cat("Model failed to converge three times: ignoring the results of this replicate.\n")

      # This NA value is skipped when computing variance-covariance matrix
      return(NA)
  }

  ass1 <- assoc1(model, weighting=weighting, rowsup=rowsup, colsup=colsup)
  ass1.orig <- assoc1(model.orig, weighting=weighting, rowsup=rowsup, colsup=colsup)

  ret <- if(inherits(ass1, c("assoc.hmskew", "assoc.hmskewL"))) find.stable.scores.hmskew(ass1, ass1.orig)
         else find.stable.scores(ass1, ass1.orig)

  # For double association models like some hmskew and yrcskew variants
  if(!is.null(assoc2)) {
      ass2 <- assoc2(model, weighting=weighting, rowsup=rowsup, colsup=colsup)
      ass2.orig <- assoc2(model.orig, weighting=weighting, rowsup=rowsup, colsup=colsup)

      ret <- c(ret, if(inherits(ass2, c("assoc.hmskew", "assoc.hmskewL"))) find.stable.scores.hmskew(ass2, ass2.orig)
                    else find.stable.scores(ass2, ass2.orig))
  }

  if(verbose && !trace)
      cat(".")

  ret
}

# Originally based on procrustes() from package vegan 2.0-4, by
# Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre,
# Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Helene Wagner.
# License GPL-2.
# In this very simplified version, we perform no centering nor scaling (handled manually with weighting)
procrustes <- function (X, Y) {
  XY <- crossprod(X, Y)
  sol <- svd(XY)

  A <- sol$v %*% t(sol$u)
  Yrot <- Y %*% A

  list(Yrot = Yrot, rotation = A, svd=sol)
}

# Compute distance between adjusted scores for this replicate to that of the original model
# We choose the permutation and sign of dimensions that minimizes the sum of squares,
# weighted by the inverse of row frequencies
find.stable.scores <- function(ass, ass.orig, detailed=FALSE) {
  weights <- 1

  nd <- ncol(ass$phi)
  nr <- nrow(ass$row)
  nc <- nrow(ass$col)
  nl <- nrow(ass$phi)
  nlr <- dim(ass$row)[3]
  nlc <- dim(ass$col)[3]
  probs <- get.probs(ass)

  sc <- adj.orig <- array(NA, dim=c(nr + nc, nd, nl))
  phi <- ass$phi

  sc[1:nr,,] <- ass$row
  sc[-(1:nr),,] <- ass$col

  adj.orig[1:nr,,] <- ass.orig$row
  adj.orig[-(1:nr),,] <- ass.orig$col

  adj <- sc
  adj <- sweep(adj, 3:2, sqrt(abs(ass$phi)), "*")
  adj.orig <- sweep(adj.orig, 3:2, sqrt(abs(ass.orig$phi)), "*")

  # If phi is negative, change sign of columns for adjusted scores so that the interpretation
  # is consistent with positive phi for plotting (plot.assoc() does the same)
  # For symmetric models, we change the sign of adjusted column scores since this makes computations easier
  adj[-(1:nr),,] <- sweep(adj[-(1:nr),, , drop=FALSE], 3:2, sign(ass$phi), "*")
  adj.orig[-(1:nr),,] <- sweep(adj.orig[-(1:nr),, , drop=FALSE], 3:2, sign(ass$phi), "*")

  # Transform arrays to matrices with one column per dimension
  adj <- aperm(adj, c(1, 3, 2))
  adj.orig <- aperm(adj.orig, c(1, 3, 2))
  dim(adj) <- dim(adj.orig) <- c((nr + nc) * nl, nd)
  procr <- procrustes(adj.orig, adj)

  adj <- procr$Yrot
  dim(adj) <- c(nr + nc, nl, nd)
  adj <- aperm(adj, c(1, 3, 2))

  for(l in 1:nl) {
      # phi for rows and column are identical since rotation is the same for both,
      # so take an average in case there are rounding errors
      phi[l,] <- margin.table(sweep(adj[,,l , drop=FALSE]^2, 1,
                                    c(probs$rp, probs$cp), "*"), 2)/2 * sign(ass$phi[l,])

      sc[1:nr,,l] <- sweep(adj[1:nr,,l, drop=FALSE], 2, sqrt(abs(phi[l,])), "/")
      sc[-(1:nr),,l] <- sweep(adj[-(1:nr),,l, drop=FALSE], 2, sqrt(abs(phi[l,])) * sign(ass$phi[l,]), "/")
  }

  # Sanity check: rebuild the association matrix and compare with the original
  lambda <- lambda.adj <- lambda.sav <- matrix(0, nr, nc)

  for(l in 1:nl) {
      lambda[] <- lambda.adj[] <- lambda.sav[] <- 0

      for(i in 1:nd) {
          lambda <- lambda + phi[l, i] * sc[1:nr, i, l] %o% sc[-(1:nr), i, l]
          lambda.adj <- lambda.adj + adj[1:nr, i, l] %o% adj[-(1:nr), i, l]

          if(dim(ass$row)[3] > 1) # Heterogeneous scores
              lambda.sav <- lambda.sav + ass$phi[l, i] * ass$row[, i, l] %o% ass$col[, i, l]
          else # Homogeneous scores
              lambda.sav <- lambda.sav + ass$phi[l, i] * ass$row[, i, 1] %o% ass$col[, i, 1]
      }

      stopifnot(isTRUE(all.equal(lambda, lambda.sav, check.attributes=FALSE, tolerance=1e-8)))
      stopifnot(isTRUE(all.equal(lambda.adj, lambda.sav, check.attributes=FALSE, tolerance=1e-8)))
  }

  if(detailed)
      list(phi=phi, row=sc[1:nr,,, drop=FALSE], col=sc[-(1:nr),,, drop=FALSE], rotation=procr$rotation)
  else
      c(phi, sc[,,1:nlr], adj)
}

# Simplified version of permutations() from the gtools 2.7.0 package,
# by Gregory R. Warnes.
# Original version by Bill Venables and cited by Matthew
# Wiener (mcw@ln.nimh.nih.gov) in an email to R-help dated
# Tue, 14 Dec 1999 09:11:32 -0500 (EST) in response to
# Alex Ahgarin <datamanagement@email.com>
perms <- function(n) {
  sub <- function(n, v) {
      if(n == 1) return(matrix(v, n, 1))

      X <- NULL

      for(i in 1:n)
        X <- rbind(X, cbind(v[i], Recall(n - 1, v[-i])))

      X
  }

  sub(n, 1:n)
}

# This class of models is special because dimensions are paired,
# and a rotation is needed rather than changing signs
find.stable.scores.hmskew <- function(ass, ass.orig) {
  weights <- 1

  nd <- ncol(ass$phi)
  nr <- nrow(ass$row)
  nc <- nrow(ass$col)
  nl <- nrow(ass$phi)
  nlr <- dim(ass$row)[3]
  nlc <- dim(ass$col)[3]

  phi <- ass$phi

  # Repeat scores once for each layer
  sc <- array(ass$row, dim=c(nr, nd, nl))
  adj <- array(ass$row, dim=c(nr, nd, nl))
  adj.orig <- array(ass.orig$row, dim=c(nr, nd, nl))

  # Compute adjusted scores
  adj <- sweep(adj, 3:2, sqrt(abs(ass$phi)), "*")
  # Where phi is negative, change signs on second dimension of each pair to get the same effect
  adj[,seq(1, nd, by=2),] <- sweep(adj[,seq(1, nd, by=2),, drop=FALSE], 3:2,
                                   sign(ass$phi)[,seq(1, nd, by=2)], "*")
  adj.orig[,seq(1, nd, by=2),] <- sweep(adj.orig[,seq(1, nd, by=2),, drop=FALSE], 3:2,
                                   sign(ass.orig$phi)[,seq(1, nd, by=2)], "*")

  # Rotate scores and return the sum of squares to the old scores
  rot <- function(angle, adj.tmp, dim) {
      rotmat <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
      adj.rot <- adj.tmp[, c(dim, dim + 1), , drop=FALSE]
      adj.rot[] <- apply(adj.rot, 3, "%*%", rotmat)
      sum(sweep((adj.rot - adj.orig[, c(dim,  dim + 1), , drop=FALSE])^2, 1, weights, "*"))
  }

  perms <- perms(nd/2)
  nperms <- nrow(perms)
  angles <- matrix(NA, nd/2, nperms)
  sq <- numeric(nperms)

  for(j in 1:nperms) {
      order <- rep(2 * perms[j,] - 1, each=2) + c(0, 1)
      adj.tmp <- adj[, order, , drop=FALSE]

      # Find optimal rotation for each pair of dimensions and apply it
      for(i in seq.int(1, nd, by=2)) {
          angle <- optim(0, rot, gr=NULL, adj.tmp, i, method="BFGS")$par
          angles[i, j] <- angle
          rotmat <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
          adj.tmp[, c(i, i + 1),] <- apply(adj.tmp[, c(i, i + 1), , drop=FALSE], 3, "%*%", rotmat)       
      }

      sq[i] <- sum(sweep((adj.orig - adj.tmp)^2, 1, weights, "*"))
  }

  # Find out which permutation allows for the smaller sum of squares,
  # choosing the best rotation for each pair of dimensions
  best.perm <- which.min(sq)
  order <- rep(2 * perms[best.perm,] - 1, 2) + c(0, 1)
  adj <- adj[, order, , drop=FALSE]
  # Phi does not need to be recomputed from adjusted scores since it does not change after a rotation
  phi <- phi[, order, drop=FALSE]

  for(i in seq.int(1, nd, by=2)) {
      best.angle <- angles[i, best.perm]
      rotmat <- matrix(c(cos(best.angle), sin(best.angle), -sin(best.angle), cos(best.angle)), 2, 2)

      sc[, c(i, i + 1),] <- apply(sc[, c(i, i + 1), , drop=FALSE], 3, "%*%", rotmat)
      adj[, c(i, i + 1),] <- apply(adj[, c(i, i + 1), , drop=FALSE], 3, "%*%",  rotmat)
  }

  # Sanity check: rebuild the association matrix and compare with the original
  lambda <- lambda.adj <- lambda.sav <- matrix(0, nr, nc)

  for(l in 1:nl) {
     lambda[] <- lambda.adj[] <- lambda.sav[] <- 0

      for(i in seq.int(1, nd, by=2)) {
          lambda <- lambda + phi[l, i] * (sc[, i + 1, l] %o% sc[, i, l] -
                                          sc[, i, l] %o% sc[, i + 1, l])
          lambda.adj <- lambda.adj + adj[, i + 1, l] %o% adj[, i, l] -
                                     adj[, i, l] %o% adj[, i + 1, l]

          if(dim(ass$row)[3] > 1) { # Heterogeneous scores
              lambda.sav <- lambda.sav + ass$phi[l, i] * (ass$row[, i + 1, l] %o% ass$row[, i, l] -
                                                          ass$row[, i, l] %o% ass$row[, i + 1, l])
          }
          else { # Homogeneous scores
              lambda.sav <- lambda.sav + ass$phi[l, i] * (ass$row[, i + 1, 1] %o% ass$row[, i, 1] -
                                                          ass$row[, i, 1] %o% ass$row[, i + 1, 1])
          }
      }

      stopifnot(isTRUE(all.equal(lambda, lambda.sav, check.attributes=FALSE, tolerance=1e-8)))
      stopifnot(isTRUE(all.equal(lambda.adj, lambda.sav, check.attributes=FALSE, tolerance=1e-8)))
  }

  # Repeat row scores to match a general association structure
  c(phi, sc[rep.int(1:nr, 2),,1:nlr], adj[rep.int(1:nr, 2),,])
}

