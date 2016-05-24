predict.ramps <- function(object, newdata,
   type = c("response", "spatial", "error", "random"), ...)
{
   ## Type of prediction to produce
   type <- match.arg(type)

   ## Create data frame containing all relevant variables
   mt <- delete.response(object$terms)
   spf <- getCovariateFormula(object$correlation)
   variance <- eval(object$call[["variance"]])
   val <- reformulate(c(labels(mt), all.vars(spf), all.vars(variance$fixed),
                        all.vars(variance$spatial), all.vars(variance$random)))
   mfdata <- model.frame(val, newdata, xlev = object$xlevels)

   ## Extract spatial coordinates
   spt <- terms(spf)
   attr(spt, "intercept") <- 0
   newcoords <- model.matrix(spt, mfdata)

   ## Construct matrix of measured and unmeasured sites
   ## Initialize the correlation structure using all sites
   val <- merge(cbind(newcoords, 1:nrow(newcoords)),
                cbind(object$coords, 1:nrow(object$coords)),
                by=1:ncol(newcoords), all.x=TRUE, all.y=FALSE)
   val <- val[order(val[, ncol(newcoords)+1]),]
   monitor <- is.na(val[, ncol(newcoords)+2])

   n <- sum(monitor)

   if (n == 0)
      stop("No new sites supplied in 'newdata'.  Spatial prediction is not",
           " supported\n\tat point-source measurement sites or at areal grid",
           " sites used in 'georamps'.")

   mfdata <- mfdata[monitor,,drop=FALSE]
   newcoords <- newcoords[monitor,,drop=FALSE]

   sites <- unique.sites(newcoords)
   correlation <- Initialize(
              eval(object$call[match("correlation", names(object$call))][[1]]),
              data = as.data.frame(rbind(object$coords, sites$coords)))

   ## Extract main effects design matrix and spatial variances indices
   xmat <- Matrix(model.matrix(mt, mfdata, attr(object$xmat, "contrasts")))
   kmat <- sites$map
   etype <- if (is.null(variance$fixed)) rep(1, n)
            else factor(getCovariate(mfdata, variance$fixed), levels(object$etype))
   retype <- if (is.null(variance$random)) rep(1, n)
             else factor(getCovariate(mfdata, variance$random), levels(object$retype))
   ztype <- if (is.null(variance$spatial)) rep(1, n)
            else factor(getCovariate(mfdata, variance$spatial), levels(object$ztype))

   ## Extract sampled correlation and variance parameters 
   phi <- params2phi(object$params, object$control)
   kappa <- params2kappa(object$params, object$control)
   sigma.e <- sqrt(kappa2kappa.e(kappa, object$control))
   sigma.re <- sqrt(kappa2kappa.re(kappa, object$control))
   if (ncol(sigma.re) == 0) sigma.re <- cbind(sigma.re, 0)
   sigma.z <- sqrt(kappa2kappa.z(kappa, object$control))

   ## K matrix components needed for mpdensity
   val <- unique.sites(object$kmat)
   xk1mat <- cBind(object$xmat, val$map)
   k2mat <- val$coords

   ## Constants for the sampling algorithm
   allz <- all(object$control$z$monitor)
   nz <- nrow(k2mat)
   n1 <- nrow(object$coords)
   n2 <- nrow(sites$coords)
   p <- ncol(object$xmat)

   y <- structure(
      matrix(0, nrow(object$params), n,
         dimnames = list(object$control$iter,
                         name.ext("fit", rownames(newcoords)))),
      coords = newcoords,
      class = c("predict.ramps", "matrix")
   )

   cat("MCMC Sampler Progress (N = ", max(object$control$iter), "):\n", sep="")
   for(i in 1:nrow(y)) {
      if (allz) {
         beta <- params2beta(object$params[i,], object$control)
         z <- k2mat %*% object$z[i,]
      } else {
         mpd <- mpdbetaz(object$params[i,], object$y, xk1mat, k2mat,
                         object$wmat, object$correlation, object$etype,
                         object$ztype, object$retype, object$weights,
                         object$control)
         BETA <- mpd$betahat + solve(mpd$uXtSiginvX, rnorm(p + nz))
         beta <- BETA[seq(length.out = p)]
         z <- BETA[seq(p + 1, length.out = nz)]
      }

      coef(correlation) <- phi[i,]
      R <- corMatrix(correlation)
      KMAT <- k2mat %*% Diagonal(x = sigma.z[i, object$ztype])
      R11 <- symmpart(KMAT %*% tcrossprod(R[1:n1, 1:n1], KMAT))
      R11ui <- solve(chol(R11))
      val <- tcrossprod(R[(n1+1):(n1+n2), 1:n1], KMAT) %*% R11ui
      zp <- rmvnorm2(1, tcrossprod(val, R11ui) %*% z,
                        R[(n1+1):(n1+n2), (n1+1):(n1+n2)] - tcrossprod(val))[1,]

      y[i,] <- as.vector(xmat %*% beta + sigma.z[i, ztype] * kmat %*% zp) +
         switch(type,
            error    = rnorm(n, 0, sigma.e[i, etype]),
            random   = rnorm(n, 0, sigma.re[i, retype]),
            response = rnorm(n, 0, sigma.e[i, etype] + sigma.re[i, retype]),
            0
         )

      print.iter(object$control$iter[i])
   }
   cat("\n")

   y
}


print.predict.ramps <- function(x, ...)
{
   attr(x, "coords") <- NULL
   print.default(x[], ...)
}
