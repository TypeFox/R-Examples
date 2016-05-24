"qua2ci.cov" <-
function(x,f, type=NULL, nsim=1000, interval=c("confidence", "none"), level=0.90,
              asnorm=FALSE, altlmoms=NULL, flip=NULL, dimless=TRUE,
              usefastlcov=FALSE, nmom=5, getsimlmom=FALSE, verbose=FALSE, ...) {
   interval <- match.arg(interval)
   if(is.null(type)) {
      message("must specify an lmomco distribution abbreviation, returning NA")
      return(NA)
   }
   if(! any(dist.list() == type)) {
      warning("invalid lmomco distribution abbreviation given, returning NA")
      return(NA)
   }
   if(interval == "none" & length(f) > 1) {
      warning("function is not vectorized in f for interval 'none', using ",
              "only the first value")
      f <- f[1]
   }
   if(! check.fs(f)) return(NA) # invalid probability
   if(! check.fs(level) | level == 1) { # invalid confidence level
      warning("level is not in [0,1), returning NA")
      return(NA)
   }

   # Notice that if either of the special extensions of alternative L-moments
   # and (or) the data have been flipped that dimensionless simulation will
   # take place.
   if(! is.null(altlmoms) & ! dimless) {
      warning("Altervative L-moments are given, so setting 'dimless=TRUE'")
      dimless <- TRUE
   }

   npar <- dist.list(type=type)
   if(nmom < npar) {
      warning("number of parameters distribution 'type' is smaller ",
              "than 'nmom' resetting 'nmom to' ", npar)
      nmom <- npar
   }

   # Note right off the bat that the data are actually flipped.
   if(! is.null(flip)) x <- flip - x  # flip the data if indicated

   xlmoms <- lmoms(x, nmom=nmom) # L-moments of the 'original' or
   # flipped data

   if(is.null(altlmoms)) altlmoms <- xlmoms # set altlmoms if not already
   if(verbose) {
      message(" Lambdas (alternate or of x) "); print(altlmoms$lambdas)
   }

   if(! are.lmom.valid(altlmoms)) {
      warning(" L-moments (either computed on 'x') or from 'altlmoms' ",
              "are invalid, returning NA")
      return(NA)
   }
   parent <- lmom2par(altlmoms, type=type)
   if(verbose) {
      message(" Parent parameters "); print(parent$para)
   }

   # Multi-variate Normal simulation of the L-moments
   if(dimless) { # Removal of dimensions and scale standardization are requested
      newx  <- (x - xlmoms$lambdas[1]) / xlmoms$lambdas[2] # standardize
      if(usefastlcov) {
         lmrcv <- lmoms.cov(newx, nmom=nmom)
      } else { # Lmoments package has a fast algorithm
         lmrcv <- Lmoments::Lmomcov(newx, rmax=nmom)
      }
      sLM <- NULL
      try(sLM   <- MASS::mvrnorm(n=nsim, c(0, 1, altlmoms$ratios[3:nmom]), lmrcv))
      if(is.null(sLM)) { print(lmrcv); stop("here") }
         MU <- altlmoms$lambdas[1] # desired mean, which could be the original mean
         L2 <- altlmoms$lambdas[2] # desired L2, which again could be the original
      # based on afore logic involving: if(is.null(altlmoms)) altlmoms <- xlmoms
      # notice that the LAMBDAS are converted, ratios will be RECOMPUTED later by
      # logic involving the standard operation using first two lambdas and the
      # trailing ratios in the vector-to-L-moments function vec2lmom().
      # vec2lmom(c(slams[1:2], slams[3:nmom]/slams[2]))
      for(i in 2:nmom) sLM[,i] <- sLM[,i]*L2      # rescale
                       sLM[,1] <- sLM[,1]*L2 + MU # rescale and shift
   } else {
      if(usefastlcov) {
         lmrcv <- lmoms.cov(x, nmom=nmom)
      } else { # Lmoments package has a fast algorithm
         lmrcv <- Lmoments::Lmomcov(x, rmax=nmom)
      }
      sLM <- NULL
      try(sLM <- MASS::mvrnorm(n=nsim, altlmoms$lambdas, lmrcv))
      if(is.null(sLM)) { print(x); print(lmrcv); stop("here") }
   }
   if(verbose) {
      message(" Var-covar matrix");              print(lmrcv)
      message(" Summary of simulated lambdas "); print(summary(sLM))
   }
   if(getsimlmom) return(sLM)
   # Confidence limit computation
   ci <- c((1-level)/2, 1-(1-level)/2) # division by 2 because two tailed
   n <- length(f)
   zz <- matrix(nrow=n, ncol=8)
   for(i in 1:n) {
      ff <- f[i]
      quasf <- sapply(1:nsim, function(i) {
               slams <- sLM[i,] # extract L-moments, then cast to lmomco style
               lmr  <- vec2lmom(c(slams[1:2], slams[3:nmom]/slams[2]))
               if(! are.lmom.valid(lmr)) return(NA) # silent warning
               par  <- lmom2par(lmr, type=type, ...)
               if(is.null(par)) return(NA)
               quaf <- NULL
               try(quaf <- ifelse(is.null(flip), par2qua(  ff, par),
                                          flip - par2qua(1-ff, par)))
               ifelse(is.null(quaf), return(NA), return(quaf))
            })
      if(interval == "none") return(quasf) # NAs will be returned by having this
      # conditional ahead of the stripping and subsequent statistics thereof

      if(any(is.na(quasf))) {
         num.na <- length(quasf[is.na(quasf)])
         warning("at least one of 'quasf' is NA, stripping ", num.na,
                 " NA values (entries) and continuing")
         quasf <- quasf[! is.na(quasf)] # This should be quite unusual
         # except for distributions having thin slices of valid L-moment
         # domain or simulations near the boundary of the distributions,
         # say above the Tau4 of the Generalized Logistic distribution.
      }
      lmr <- lmoms(quasf, nmom=3) # the L-moments of the simulated quantiles
      mu  <- lmr$lambdas[1]; lscale <- lmr$lambdas[2]
      ifelse(asnorm, z <- qnorm(ci, mean=mu, sd=sd(quasf)), # assume Normality
                     z <- dat2bernqua(ci,  quasf))  # nonparametric
      fit <- ifelse(is.null(flip), par2qua(  ff, parent),
                            flip - par2qua(1-ff, parent))
      zz[i,] <- c(ff, z[1], fit, z[2],
                  median(quasf), mu, var(quasf), lscale)
   }
   # A quick review if nsim is large enough is that the qua_med/qua_mean
   # columns should have values similar to the fit.  The question of symmetry
   # of the quantile estimates can be deduced by qua_med and qua_mean
   # being nearly identical numbers.
   colnames(zz) <- c("nonexceed", "lwr", "fit", "upr",
                     "qua_med", "qua_mean", "qua_var", "qua_lam2")
   zz <- as.data.frame(zz)
   return(zz)
}
