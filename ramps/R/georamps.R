georamps <- function(fixed, random, correlation, data = sys.frame(sys.parent()),
   subset, weights, variance = list(fixed = ~ 1, random = ~ 1, spatial = ~ 1),
   aggregate = list(grid = NULL, blockid = ""), kmat = NULL,
	control = ramps.control(...), contrasts = NULL, ...)
{
   ## Create data frame containing all relevant variables
   ## Random effects are added later
   call <- match.call()
   val <- c(all.vars(fixed), all.vars(variance$fixed),
            all.vars(variance$spatial), all.vars(variance$random))
   if (!is.null(aggregate$grid)) val <- c(val, aggregate$blockid)
   spvars <- all.vars(getCovariateFormula(correlation))
   val <- reformulate(c(val, spvars))
   mfargs <- list(formula = val, data = data, weights = call[["weights"]],
                  subset = call[["subset"]], na.action = na.pass)
   mfdata <- do.call("model.frame", mfargs)

   ## Zero-out the coordinates for aggregate measurements
   if (nchar(aggregate$blockid) > 0) {
      val <- c(aggregate$blockid, spvars)
      if (!all(val %in% colnames(aggregate$grid)))
         stop("Coordinates and 'blockid' must be given in 'grid'")
      aggregate$grid <- na.omit(aggregate$grid[,val])
      idx <- is.element(mfdata[,aggregate$blockid],
                        aggregate$grid[,aggregate$blockid])
      mfdata[idx, spvars] <- 0
   }

   ## Remove incomplete records from the data frame
   mfdata <- na.omit(mfdata)

   ## Extract weights
   weights <- model.weights(mfdata)

   ## Drop the model frame attributes that are no longer needed
   attr(mfdata, "terms") <- NULL

   # Extract response vector and main effects design matrix
   mf <- model.frame(fixed, data = mfdata, drop.unused.levels = TRUE)
   mt <- attr(mf, "terms")
   y <- model.response(mf, "numeric")
   val <- model.matrix(mt, mf, contrasts)
   xmat <- as(val, "sparseMatrix")
   attr(xmat, "contrasts") <- attr(val, "contrasts")

   ## Indices to map measurement error variances
   variance$fixed <- if (is.null(variance$fixed)) factor(rep(1, nrow(mfdata)))
                     else factor(getCovariate(mfdata, variance$fixed))

   ## Structures for latent spatial parameters
   if (missing(correlation)) {
      stop("Unspecified correlation structure")
   } else {
      ## Create matrix of unique coordinates for latent parameters
      spt <- terms(getCovariateFormula(correlation))
      attr(spt, "intercept") <- 0

      if (is.null(aggregate$grid)) {
         idx1 <- rep(TRUE, nrow(mfdata))
         val <- model.matrix(spt, mfdata)
         idx2 <- NULL
      } else {
         idx1 <- !is.element(mfdata[,aggregate$blockid],
                             aggregate$grid[,aggregate$blockid])
         val <- model.matrix(spt, mfdata[idx1,,drop=FALSE])
         idx2 <- is.element(aggregate$grid[,aggregate$blockid],
                            mfdata[,aggregate$blockid])
         val <- rbind(val, model.matrix(spt, aggregate$grid[idx2,,drop=FALSE]))
      }
      sites <- unique.sites(val)

      ## Logical vector indicating z values to be monitored
      if (is.logical(control$z$monitor)) {
         control$z$monitor <- rep(control$z$monitor,
                                  length.out = nrow(sites$coords))
      } else {
         idx <- colnames(sites$coords)
         if (!all(idx %in% colnames(control$z$monitor)))
            stop("Coordinate names not found in 'z' monitor")
         val <- unique.sites(control$z$monitor[,idx,drop=FALSE])
         val <- merge(cbind(sites$coords, 1:nrow(sites$coords)),
                      cbind(val$coords, 1), by=idx, all.x = TRUE)
         n <- length(idx)
         control$z$monitor <- !is.na(val[order(val[, n+1]), n+2])
      }

      ## Order latent parameters as (z$monitor == T, z$monitor == F)
      idx <- order(control$z$monitor, decreasing = TRUE)
      control$z$monitor <- control$z$monitor[idx]
      sites$coords <- sites$coords[idx,,drop=FALSE]
      sites$map <- sites$map[,idx,drop=FALSE]

      ## Initialize correlation structure
      correlation <- Initialize(correlation, data = as.data.frame(sites$coords))

      ## Matrix to map latent parameters to observed data
      if (is.null(kmat)) {
         k <- sites$map
         kmat <- Matrix(0, nrow(mfdata), nrow(sites$coords))
         kmat[idx1,] <- k[seq(length.out = sum(idx1)),]
         if (length(idx2) > 0) {
            idx <- aggregate$grid[idx2, aggregate$blockid]
            val <- sort(unique(idx))
            kmap <- Matrix(0, length(val), length(idx))
            kmap[nrow(kmap) * (seq(idx) - 1) + match(idx, val)] <- 1
            kmat[match(val, mfdata[, aggregate$blockid]),] <-
               (kmap / tabulate(idx)) %*%
               k[seq(sum(idx1) + 1, length.out = sum(idx2)),]
         }
      } else {
         n <- c(nrow(mfdata), nrow(sites$coords))
         if (!(is(kmat, "matrix") || is(kmat, "Matrix")) || any(dim(kmat) != n))
            stop("Supplied 'kmat' must be a matrix object of dimension ",
                 n[1], " x ", n[2])
         kmat <- as(kmat, "sparseMatrix")
         k <- abs(kmat)
         if (any(rowSums(k) == 0))
            stop("Supplied 'kmat' should not contain rows of zeros")
         if (any(colSums(k) == 0))
            stop("Supplied 'kmat' should not contain columns of zeros")
      }

      ## Indices to map spatial variances
      val <- if (is.null(variance$spatial)) factor(rep(1, nrow(mfdata)))
             else factor(getCovariate(mfdata, variance$spatial))
      idx <- unlist(apply(as.numeric(val) * (kmat != 0), 2, unique))
      idx <- idx[idx > 0]
      if (length(idx) != ncol(kmat))
         stop("Unsupported latent spatial structure.  Different spatial",
              " variances assigned to measurements from the same site.")
      else variance$spatial <- as.factor(levels(val)[idx])
   }

   ## Structures for random effects parameters
   if (missing(random)) {
      wmat <- Matrix(numeric(0), nrow(mfdata), 0)
   } else {
      ## Matrix to map random effects to observed data
      w <- factor(getGroups(data, random)[as.numeric(rownames(mfdata))])
      wmat <- Matrix(0, length(w), nlevels(w))
      wmat[na.omit(seq(w) + nrow(wmat) * (as.numeric(w) - 1))] <- 1

      ## Indices to map random effects variances
      val <- if (is.null(variance$random)) factor(rep(1, nrow(mfdata)))
             else factor(getCovariate(mfdata, variance$random))
      idx <- unlist(apply(as.numeric(val) * (wmat > 0), 2, unique))
      idx <- idx[idx > 0]
      if (length(idx) != ncol(wmat))
         stop("Unsupported random effects structure.  Different random effects",
              " variances assigned to measurements within the same group.")
      else variance$random <- as.factor(levels(val)[idx])
   }

   ## Default values for weights if not supplied
   if (is.null(weights)) weights <- rowSums(as(kmat, "lsparseMatrix"))

   ## Check parameter specifications against supplied data
   if (length(control$beta) != (n <- ncol(xmat)))
      stop("'beta' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$sigma2.e) != (n <- nlevels(variance$fixed)))
      stop("'sigma2.e' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$phi) != (n <- length(correlation)))
      stop("'phi' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$sigma2.z) != (n <- nlevels(variance$spatial)))
      stop("'sigma2.z' parameter specification in 'ramps.control' must be of",
           " length ", n)
   if (length(control$sigma2.re) != (n <- nlevels(variance$random)))
      stop("'sigma2.re' parameter specification in 'ramps.control' must be of",
           " length ", n)

   ## Set a single tuning parameter for the sigma2 parameters
   val <- min(sigma2tuning(control))
   if (length(control$sigma2.e)) control$sigma2.e$tuning[] <- val
   if (length(control$sigma2.z)) control$sigma2.z$tuning[] <- val
   if (length(control$sigma2.re)) control$sigma2.re$tuning[] <- val

   ## Obtain MCMC samples from ramps engine
   val <- ramps.engine(y, xmat, kmat, wmat, correlation, variance$fixed,
                       variance$spatial, variance$random, weights, control)

   structure(
      list(params = as.mcmc(val$params), z = as.mcmc(val$z),
           loglik = val$loglik, evals = val$evals, call = call, y = y,
           xmat = xmat, terms = attr(mf, "terms"),
           xlevels = .getXlevels(mt, mf), etype = variance$fixed,
           weights = weights, kmat = kmat, correlation = correlation,
           coords = sites$coords, ztype = variance$spatial, wmat = wmat,
           retype = variance$random, control = control),
      class = "ramps")
}


print.ramps <- function(x, ...)
{
   cat("\nCall: ", paste(deparse(x$call), collapse = "\n"), "\n")

   params <- colnames(x$params)
   cat("\nCoefficients:\n")
   if (length(tmp <- params2beta(params, x$control)) > 0)
      print.default(tmp, print.gap = 2, quote = FALSE)

   sigma2 <- params2kappa(params, x$control)

   n <- sum(rowSums(as(x$kmat, "lsparseMatrix")) > 1)
   cat("\nMeasurements\n",
       " N = ", length(x$y), "\n",
       " Point Source = ", length(x$y) - n, "\n",
       " Areal = ", n, "\n",
       " Error Variance: ",
       paste(kappa2kappa.e(sigma2, x$control), collapse = " "), "\n", sep = "")

   cat("\nLatent Spatial Process\n",
       " Sites = ", ncol(x$kmat), "\n",
       " Correlation: ", class(x$correlation)[1], "(",
          paste(params2phi(params, x$control), collapse = ", "), ")\n",
       " Variance: ",
       paste(kappa2kappa.z(sigma2, x$control), collapse = " "), "\n", sep = "")

   if (ncol(x$wmat) > 0) {
      cat("\nRandom Effects\n",
          " N = ", ncol(x$wmat), "\n",
          " Variance: ",
          paste(kappa2kappa.re(sigma2, x$control), collapse = " "), "\n",
          sep = "")
   }

   n <- nrow(x$params)
   rn <- rownames(x$params)[1:min(n, 3)]
   if (n > 4) rn <- c(rn, "...")
   if (n > 3) rn <- c(rn, rownames(x$params)[n])
   cat("\nMCMC Output\n",
       " Saved Samples = ", n, " (", paste(rn, collapse = ", "), ")\n",
       " Slice Evaluations = ", sum(x$evals), "\n", sep = "")

   invisible(x)
}


summary.ramps <- function(object, ...)
{
   summary(object$params, ...)
}
