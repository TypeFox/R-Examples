################################################################################
## ramps.engine - MCMC engine
##
## Arguments:
##    y       - numerical response vector (n x 1)
##    xmat    - covariate design Matrix (n x p)
##    kmat    - design Matrix of spatial random effect (n by nz)
##    wmat    - design Matrix of non-spatial random effect (n by q)
##    spcor   - initialized (nlme) spatial correlation structure
##    etype   - factor indexing the measurement variances (n x 1)
##    ztype   - factor indexing the spatial variances (nz x 1)
##    retype  - factor indexing the random effect variances (q x 1)
##    weights - numerical vector by which to weight the measurement error
##              variance
##    control - ramps.control object
################################################################################

ramps.engine <- function(y, xmat, kmat, wmat, spcor, etype, ztype, retype,
                         weights, control)
{
   ## Global constants
   n <- length(y)
   p <- length(control$beta)
   nzp <- sum(control$z$monitor)
   iter <- control$iter

   ## Return objects
   val <- c(name.ext("phi", 1:length(control$phi)),
            name.ext("sigma2.e", levels(etype)),
            name.ext("sigma2.z", levels(ztype)),
            name.ext("sigma2.re", levels(retype)), colnames(xmat))
   params <- matrix(NA, length(iter), length(val), dimnames = list(iter, val))

   z <- matrix(NA, length(iter), nzp, dimnames = list(iter, NULL))
   if (nzp > 0) colnames(z) <- paste("z", 1:nzp, sep="") 

   loglik <- structure(rep(NA, length(iter)), names = iter)
   evals <- structure(rep(NA, length(iter)), names = iter)

   ## Initialize external output files
   if (!control$expand) {
      write.header(colnames(params), control$file$params)
      write.header(colnames(z), control$file$z)
   }

   ## Indicies of spatial parameters to monitor
   sites <- unique.sites(kmat)
   zidx <- seq(length.out = nzp)
   idx <- match(zidx, as.vector((sites$coords == 1) %*% 1:ncol(kmat)))

   ## Determine whether additional prediction is needed for spatial parameters
   if (any(is.na(idx))) {
      pred <- "mpdpred"

      ## Indices for sampled spatial paramters from mpdpred
      zidx <- p + zidx

      ## Extract matrix components for mpdpred
      val <- as.vector(kmat %*% as.numeric(!control$z$monitor)) == 0
      idx <- order(val, decreasing = TRUE)
      nr1 <- sum(val)
      r2 <- seq(nr1 + 1, length.out = nrow(kmat) - nr1)
      c1 <- seq(length.out = nzp)
      c2 <- seq(nzp + 1, length.out = ncol(kmat) - nzp)

      ## Reorder existing data structures by idx
      y <- y[idx]
      xmat <- xmat[idx, , drop = FALSE]
      sites$map <- sites$map[idx, , drop = FALSE]
      wmat <- wmat[idx, , drop = FALSE]
      weights <- weights[idx]
      etype <- etype[idx]

      ## Construct additional structures for mpdpred
      Y <- c(y, rep(0, nzp))
      X <- bdiag(xmat, Diagonal(x = -1, nzp))
      X[1:n, zidx] <- kmat[idx, c1]
      k22mat <- kmat[idx[r2], c2]

   } else if ((nzp > 0) && (control$mpdfun != "mpdbetaz")) {
      pred <- "mpdbetaz"

      ## Indices for sampled spatial paramters from mpdensity
      zidx <- p + idx

      ## Construct matrices for mpdensity
      xk1mat <- cBind(xmat, sites$map)
      k2mat <- sites$coords

   } else {
      pred <- ""

      ## Indices for slice sampled spatial paramters
      zidx <- p + idx
   }

   ## List of arguments for mpdensity call
   args <- list(theta = c(control$phi$init, prop.table(sigma2init(control))),
                y = y, wmat = wmat, spcor = spcor,
                etype = etype, ztype = ztype, retype = retype,
                weights = weights, control = control)
   switch(control$mpdfun,
      mpdbeta = {
         args$xmat <- xmat
         args$kmat <- kmat
      },
      mpdbetaz = {
         args$xk1mat <- cBind(xmat, sites$map)
         args$k2mat <- sites$coords
      }
   )

   ## First mpdensity evaluation
   curreval <- do.call(control$mpdfun, args)

   ## Arguments for subsequent slice sampler calls
   args$mpdfun <- get(control$mpdfun)
   args$log <- TRUE
   args$f.theta <- curreval$value

   ## Main MCMC sampler loop
   idx <- 1
   cat("MCMC Sampler Progress (N = ", control$expand + max(iter), "):\n",
       control$expand, sep = "")
   for (i in 1:max(iter)) {
      ## Draw phi and kappa
      curreval <- do.call("sliceSimplex", args)

      args$theta <- curreval$theta
      args$f.theta <- curreval$value

      if (i == iter[idx]) {
         ## Draw variance parameters sigma2
         kappa <- params2kappa(args$theta, control)
         as2 <- sum(sigma2shape(control)) + (n - p) / 2.0
         bs2 <- sum(sigma2scale(control) / kappa) + sum(curreval$quadform) / 2.0
         sigma2.tot <- 1.0 / rgamma(1, as2, bs2)

         ## Draw beta and z parameters
         switch(pred,
            mpdbetaz = {
              val <- mpdbetaz(args$theta, y, xk1mat, k2mat, wmat, spcor,
                              etype, ztype, retype, weights, control)
              BETA <- val$betahat +
                solve(val$uXtSiginvX, rnorm(length(val$betahat), sd=sqrt(sigma2.tot)))
              BETA0 <- head(BETA, p)
            },
            mpdpred = {
              val <- mpdpred(args$theta, Y, X, k22mat, wmat, spcor, etype,
                             ztype, retype, weights, control)
              BETA <- val$betahat +
                solve(val$uXtSiginvX, rnorm(length(val$betahat), sd=sqrt(sigma2.tot)))
              BETA0 <- head(BETA, p)
              if(length(curreval$betahat) > p) {
                BETA0 <- if (nzp != ncol(kmat)) NULL
                         else rBind(BETA0, kmat %*% tail(BETA, -p))
              }
            },
            {
              BETA <- curreval$betahat + solve(curreval$uXtSiginvX,
                        rnorm(length(curreval$betahat), sd=sqrt(sigma2.tot)))
              BETA0 <- BETA
            }
         )
         
         ## Calculate the likelihood
         if (!is.null(BETA0)) {
           loglik[idx] <- -0.5 * n * log(sigma2.tot) - curreval$logsqrtdet -
             (crossprod(curreval$uXtSiginvX %*% (curreval$betahat - BETA0))[1]
              + curreval$quadform[1]) / (2.0 * sigma2.tot)
         }
         
         ## Save model parameters
         val <- c(params2phi(args$theta, control),
                  sigma2.tot * kappa2kappa.e(kappa, control),
                  sigma2.tot * kappa2kappa.z(kappa, control),
                  sigma2.tot * kappa2kappa.re(kappa, control),
                  BETA[seq(length.out = p)])
         params[idx, ] <- val
         write.params(control$expand + i, val, control$file$params)
         ## Save latent spatial parameters
         val <- BETA[zidx]
         z[idx, ] <- val
         write.params(control$expand + i, val, control$file$z)

         ## Save number of new slice evaluations
         evals[idx] <- curreval$newevals

         ## Advance index for the output structures
         idx <- idx + 1
      }

      print.iter(control$expand + i)
   }
   cat("\n")

   list(params = params, z = z, loglik = loglik - 0.5 * n * log(2.0 * pi),
        evals = evals)
}
