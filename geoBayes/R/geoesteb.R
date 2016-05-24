##' Empirical Bayes estimation for SGLMM
##'
##' Currently the following spatial correlation functions are
##' implemented. Below, \eqn{h} denotes the distance between
##' locations, \eqn{d} is the dimensionality of the locations,
##' \eqn{\phi}{phi} is the spatial range parameter and \eqn{\kappa}{kappa} is an
##' additional parameter. The correlation \eqn{r(u)} beween locations
##' with distance \eqn{u} apart is
##' \describe{
##' \item{Matern}{\deqn{r(h) =
##'    \frac{1}{2^{\kappa-1}\Gamma(\kappa)}(\frac{h}{\phi})^\kappa
##'    K_{\kappa}(\frac{h}{\phi})}{r(h) =
##'    (1/(2^(kappa-1) * Gamma(kappa))) * ((h/phi)^kappa) *
##' K_{kappa}(h/phi)}}
##'  \item{spherical}{
##'  \deqn{r(h) = \left\{ \begin{array}{ll}
##'    1 - 1.5\frac{h}{\phi} + 0.5(\frac{h}{\phi})^3
##'    \mbox{ , if $h$ < $\phi$} \cr
##'    0    \mbox{ ,  otherwise}
##'    \end{array} \right.}{r(h) = 
##'    1 - 1.5 * (h/phi) + 0.5(h/phi)^3   if h < phi , 
##'    0   otherwise}
##' Note that this is a valid correlation only for \eqn{d \leq 3}{d <=
##' 3}.}
##'  \item{powerexponential}{
##'  \deqn{r(h) = \exp\{-(\frac{h}{\phi})^\kappa\}
##'    }{r(h) = exp{-(h/phi)^kappa}}
##' Note that this is a valid correlation only for \eqn{0 < \kappa
##' \leq 2}{0 < kappa <= 2}.}
##' }
##'
##' The GEV (Generalised Extreme Value) link is defined by \deqn{\mu =
##' 1 - \exp\{-\max(0, 1 + \nu x)^{\frac{1}{\nu}}\}}{mu = 1 -
##' \exp[-max(0, 1 + nu x)^(1/nu)]} for any real \eqn{\nu}{nu}. At
##' \eqn{\nu = 0}{nu = 0} it reduces to the complementary log-log
##' link.
##' @title Empirical Bayes estimation for SGLMM
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations (see the example in
##' \code{\link{mcsglmm}} for how to do this using
##' \code{\link{stackdata}}). At the observed locations the response
##' is assumed to be a total of replicated measurements. The number of
##' replications is inputted using the argument \code{weights}.
##' @param family The distribution of the data. The
##' \code{"GEVbinomial"} family is the binomial family with link the
##' GEV link (see Details).
##' @param data An optional data frame containing the variables in the
##' model.
##' @param weights An optional vector of weights. Number of replicated
##' samples for Gaussian and gamma, number of trials for binomial,
##' time length for Poisson.
##' @param subset An optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##' with the coordinates of the sampled locations.
##' @param parskel A data frame with the components "linkp", "phi",
##' "omg", and "kappa", corresponding to the link function, the
##' spatial range, the relative nugget, and the spatial smoothness
##' parameters. The latter can be omitted if not used in the
##' correlation function. Let k denote the number of rows. Then, k
##' different MCMC samples will be taken from the models with
##' parameters fixed at those values. For a square grid the output
##' from the function \code{\link[base]{expand.grid}} can be used
##' here.
##' @param paroptim A named list with the components "linkp", "phi",
##' "omg", "kappa". Each component must be numeric with length 1, 2,
##' or 3 with elements in increasing order but for the binomial family
##' linkp is also allowed to be the character "logit" and "probit". If
##' its length is 1, then the corresponding parameter is considered to
##' be fixed at that value. If 2, then the two numbers denote the
##' lower and upper bounds for the optimisation of that parameter
##' (infinities are allowed). If 3, these correspond to lower bound,
##' starting value, upper bound for the estimation of that parameter.
##' @param corrfcn Spatial correlation function. See Details.
##' @param Nout A scalar or vector of size k. Number of MCMC samples
##' to take for each run of the MCMC algorithm for the estimation of
##' the Bayes factors. See argument \code{parskel}.
##' @param Nthin A scalar or vector of size k. The thinning of the
##' MCMC algorithm for the estimation of the Bayes factors.
##' @param Nbi A scalar or vector of size k. The burn-in of the MCMC
##' algorithm for the estimation of the Bayes factors.
##' @param Npro A scalar. The number of Gibbs samples to take for
##' estimation of the conjugate parameters and for prediction at the
##' unsampled locations while the other parameters are fixed at their
##' empirical Bayes estimates.
##' @param Nprt The thinning of the Gibbs algorithm for the estimation
##' of the conjugate parameters and for prediction.
##' @param Nprb The burn-in of the Gibbs algorithm for the estimation
##' of the conjugate parameters and for prediction.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##' matrix. Can be a scalar, vector or matrix. The first two imply a
##' diagonal with those elements. Set this to 0 to indicate a flat
##' improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##' partial sill parameter.
##' @param zstart Optional starting value for the MCMC for the GRF.
##' This can be either a scalar, a vector of size n where n is the
##' number of sampled locations, or a matrix with dimensions n by k
##' where k is the number of the skeleton points in \code{parskel}.
##' @param dispersion The fixed dispersion parameter.
##' @param bfsize1 A scalar or vector of length k with all integer
##' values or all values in (0, 1]. How many samples (or what
##' proportion of the sample) to use for estimating the Bayes factors
##' at the first stage. The remaining sample will be used for
##' estimating the Bayes factors in the second stage. Setting it to 1
##' will perform only the first stage.
##' @param reference An integer between 1 and k. Which model to be
##' used as a reference, i.e. the one that goes in the denominator of
##' the Bayes factors.
##' @param bfmethod Which method to use to calculate the Bayes
##' factors: Reverse logistic or Meng-Wong.
##' @param transf Whether to use the transformed sample mu for the
##' computations. Otherwise it uses z.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param control A list of control parameters for the optimisation.
##' See \code{\link[stats]{optim}}.
##' @param verbose Whether to print messages when completing each
##' stage on screen.
##' @return A list with components
##' \itemize{
##' \item \code{parest} The parameter estimates
##' \item \code{skeleton} The skeleton points used with the corresponding
##' logarithm of the Bayes factors at those points. 
##' \item \code{optim} The output from the \code{\link[stats]{optim}}
##' function. 
##' \item \code{mcmcsample} The MCMC samples for the remaining
##' parameters and the random field. These samples correspond to the
##' Gibbs and Metropolis-Hasting samples after fixing the parameters
##' estimated by empirical Bayes at their empirical Bayes estimates.
##' \item \code{sys_time} The time taken to complete the MCMC
##' sampling, calculation of the importance weights, the
##' optimization and the final MCMC sampling. 
##' }
##' @examples
##' \dontrun{
##' data(rhizoctonia)
##'  
##' ### Define the model
##' corrf <- "spherical"
##' kappa <- 0
##' ssqdf <- 1
##' ssqsc <- 1
##' betm0 <- 0
##' betQ0 <- .01
##'  
##' ### Skeleton points
##' philist <- c(100,140,180)
##' linkp <- "logit"
##' omglist <- c(0,.5,1)
##' parlist <- expand.grid(phi = philist, linkp = linkp, omg = omglist,
##'                        kappa = kappa)
##' paroptim <- list(linkp = linkp, phi = c(100, 200), omg = c(0, 2),
##'                  kappa = kappa)
##'  
##' ### MCMC sizes
##' Nout <- Npro <- 100
##' Nthin <- Nprt <- 1
##' Nbi <- Nprb <- 0
##'  
##' est <- ebsglmm(Infected ~ 1, 'binomial', rhizoctonia, weights = Total,
##'                atsample = ~ Xcoord + Ycoord, parskel = parlist,
##'                paroptim = paroptim, corrfcn = corrf, 
##'                Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                Npro = Npro, Nprt = Nprt, Nprb = Nprb, 
##'                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
##'                dispersion = 1, useCV=TRUE)
##'}
##' @references Roy, V., Evangelou, E., and Zhu, Z. (2015). Efficient estimation
##' and prediction for the Bayesian spatial generalized linear mixed
##' model with flexible link functions. \emph{Biometrics}.
##' \url{http://dx.doi.org/10.1111/biom.12371}
##' @importFrom sp spDists
##' @importFrom stats model.matrix model.response model.weights as.formula
##' @export 
ebsglmm <- function (formula,
                     family = c("gaussian", "binomial", "poisson", "Gamma",
                       "GEV.binomial", "GEVD.binomial",
                       "Wallace.binomial"),
                     data, weights, subset, atsample, parskel, paroptim,
                     corrfcn = c("matern", "spherical", "powerexponential"), 
                     Nout, Nthin = 1, Nbi = 0, Npro, Nprt = 1, Nprb = 0, 
                     betm0, betQ0, ssqdf, ssqsc,
                     zstart, dispersion = 1,
                     bfsize1 = 0.8, reference = 1, bfmethod = c("RL", "MW"), 
                     transf = FALSE, useCV = TRUE, longlat = FALSE, 
                     control = list(), verbose = TRUE) {

  ## Family
  family <- match.arg(family)
  ifam <- match(family, eval(formals()$family))

  ## Logical input
  useCV <- isTRUE(useCV)
  longlat <- isTRUE(longlat)
  verbose <- isTRUE(verbose)

  ## Correlation function
  corrfcn <- match.arg(corrfcn)
  icf <- match(corrfcn, eval(formals()$corrfcn))
  needkappa <- corrfcn %in% c("matern", "powerexponential")

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  if (!all(is.finite(FF))) stop ("Non-finite values in the design matrix")
  p <- NCOL(FF)
  yy <- unclass(model.response(mf))
  if (!is.vector(yy)) {
    stop ("The response must be a vector")
  }
  yy <- as.double(yy)
  ll <- model.weights(mf)

  ## All locations
  locvars <- all.vars(atsample)
  formula1 <- as.formula(paste('~', paste(c(locvars, all.vars(formula)),
                                          collapse = ' + ')))
  mfc1 <- mfc
  mfc1$formula <- formula1
  mf1 <- eval(mfc1, parent.frame())
  m <- match(locvars, names(mf1))
  loc <- as.matrix(mf1[, m])
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Check corrfcn with loc
  if (corrfcn == "spherical" & NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  n <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, n) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (family %in% c("binomial", "GEV.binomial", "GEVD.binomial",
                    "Wallace.binomial")) {
    l <- l - y # Number of failures
  }
  F <- FF[ii, , drop = FALSE]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  n0 <- sum(!ii)
  if (n0 > 0) {
    F0 <- FF[!ii, , drop = FALSE]
    dmdm0 <- sp::spDists(loc[ii, , drop = FALSE], loc[!ii, , drop = FALSE],
                         longlat = longlat)
  } else {
    F0 <- dmdm0 <- numeric(0)
    dim(F0) <- c(0, p)
    dim(dmdm0) <- c(n, 0)
  }

  ## Priors
  if (all(is.finite(betQ0[upper.tri(betQ0)]))) {
    if (length(betQ0) == 1 && betQ0[1] == 0) {
      ## Uniform prior
      betQ0 <- matrix(0, p, p)
      betm0 <- rep(0, p)
    } else if (length(betQ0) == 1 || length(betQ0) == p) {
      if (any(betQ0 <= 0)) stop ('betQ0 not > 0')
      betQ0 <- diag(betQ0, p, p)
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      }
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Read and check parskel
  parskel <- .check_pargrid(parskel, family, corrfcn)
  nruns <- nrow(parskel)
  linkp <- parskel$linkp
  phi <- parskel$phi
  omg <- parskel$omg
  kappa <- parskel$kappa
  nu <- parskel$nu
  parnmall <- c("linkp", "phi", "omg", "kappa")
  parnm <- c("linkp", "phi", "omg", if(needkappa) "kappa")

  ## Method for computing the Bayes factors
  bfmethod <- match.arg(bfmethod)
  imeth <- match(bfmethod, eval(formals()$bfmethod))

  ## Read and check paroptim argument
  if (!is.list(paroptim)) {
    stop ("Argument paroptim must be a list")
  }
  if (!all(parnm %in% names(paroptim))) {
    stop (paste("Named components in the argument paroptim",
                "must be", paste(parnm, collapse = " ")))
  } else {
    paroptim <- paroptim[parnm]
  }
  if (!needkappa) paroptim$kappa <- 0
  lower <- upper <- pstart <- rep.int(NA, 4)
  estim <- rep.int(FALSE, 4)
  linkpe <- paroptim[["linkp"]]
  if (is.character(linkpe) | is.factor(linkpe)) {
    if (family == "binomial") {
      if (linkpe == "logit") {
        if (all(linkp == "logit")) {
          pstart[1] <- -1
          estim[1] <- FALSE
        } else {
          stop ("Link in paroptim inconsistent with link in parskel")
        }
      } else if (linkpe == "probit") {
        if (all(linkp == "probit")) {
          pstart[1] <- 0
          estim[1] <- FALSE
        } else {
          stop ("Link in paroptim inconsistent with link in parskel")
        }
      } else {
        stop ("Character link for the binomial family can be either
\"logit\" or \"probit\" in argument paroptim")
      }
    } else stop ("Character link in argument paroptim is only used for
the binomial family")
  } else if (is.numeric(linkpe)) {
    if (length(linkpe) < 1 | length(linkpe) > 3) {
      stop ("The length for a numeric component in paroptim must be 1, 2, or 3")
    } else if (length(linkpe) == 1) {
      pstart[1] <- linkpe
      estim[1] <- FALSE
    } else if (length(linkpe) == 2) {
      pstart[1] <- NA
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[2]
      if (lower[1] >= upper[1]) {
        stop ("The lower bound must be less than the upper bound for linkp
in paroptim")
      }
    } else {
      pstart[1] <- linkpe[2]
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[3]
      if (lower[1] > estim[1] | estim[1] > upper[1] | lower[1] == upper[1]) {
        stop ("The elements in the component linkp in paroptim must be ordered")
      }
    }
  } else {
    stop ("The element linkp in paroptim must be either numeric or character")
  }
  for (i in 2:4) {
    ppp <- paroptim[[parnmall[i]]]
    if (!is.numeric(ppp)) {
      stop(paste("The element", parnmall[i], "in paroptim must be numeric"))
    }
    lppp <- length(ppp)
    if (lppp < 1 | lppp > 3) {
      stop (paste("The element", parnmall[i], "in paroptim must have 1, 2, or 3
components"))
    }
    if (lppp == 1) {
      pstart[i] <- ppp
      estim[i] <- FALSE
    } else if (lppp == 2) {
      pstart[i] <- NA
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[2]
      if (lower[i] >= upper[i]) {
        stop (paste("The lower bound must be less than the upper bound for",
                     parnmall[i], "in paroptim"))
      }
    } else {
      pstart[i] <- ppp[2]
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[3]
      if (lower[i] > pstart[i] | pstart[i] > upper[i] | lower[i] == upper[i]) {
        stop (paste("The elements in the component", parnmall[i],
                     "in paroptim must be ordered"))
      }
    }
  }

  ## MCMC samples
  Nout <- rep(as.integer(Nout), length.out = nruns)
  Nbi <- rep(as.integer(Nbi), length.out = nruns)
  Nthin <- rep(as.integer(Nthin), length.out = nruns)
  Ntot <- sum(Nout)
  z <- matrix(0, n, Ntot)
  lglk <- numeric(Ntot)
  tsqsc <- tsqdf <- 0

  ## Starting value for z
  if (missing(zstart)) {
    zstart <- switch(family,
                     binomial =, Wallace.binomial =,
                     GEV.binomial = (y+.5)/(y+l+1),
                     GEVD.binomial = (l+.5)/(y+l+1), 
                     poisson = (y+.5)/(l+1), gaussian = y/l)
    zstart <- sapply(1:nruns,
                     function (i) linkfcn(zstart, linkp[i], family))
  }
  z[, 1 + c(0, cumsum(Nout[-nruns]))] <- zstart

  bfsize1 <- as.double(bfsize1)
  if (length(bfsize1) > nruns) {
    warning ("The number of elements in bfsize1 exceeds the number of runs;
the excess elements will be discarded")
    bfsize1 <- bfsize1[1:nruns]
  }
  if (any(bfsize1 <= 0)) {
    stop ("Argument bfsize1 must be positive")
  } else if (all(bfsize1 <= 1)) {
    Nout1 <- as.integer(bfsize1*Nout)
    if (any(Nout1 == 0)) stop ("Calculated 0 sizes; give a larger bfsize1")
  } else if (all(bfsize1 >= 1)) {
    Nout1 <- as.integer(bfsize1)
    if (any(Nout1 > Nout)) stop ("The bfsize1 exceeds the number of samples")
  } else {
    stop ("Argument bfsize1 is a mix of proportions and sizes")
  }
  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)

  tsq <- if (ifam > 0) dispersion else tsqsc

  ## RUN MCMC
  tm1 <- system.time({
    RUN1 <- .Fortran("samplemulti",
                     lglk = lglk, z = z, mu = z, 
                     as.double(phi), as.double(omg), as.double(y),
                     as.double(l),
                     as.double(F), as.double(betm0), as.double(betQ0),
                     as.double(ssqdf), as.double(ssqsc), as.double(kappa),
                     as.integer(icf), as.double(nu), as.double(tsqdf),
                     as.double(tsq),
                     as.double(dm), as.integer(Ntot), as.integer(Nout),
                     as.integer(Nbi), as.integer(Nthin), as.integer(n),
                     as.integer(p), as.integer(nruns), as.integer(ifam))
  })

  if (verbose) {
    message ("Completed MCMC sampling: ", round(tm1[1]), " sec")
  }

  ## Prepare data for first stage
  z <- RUN1$z
  mu <- RUN1$mu
  lz <- unlist(lapply(1:nruns, function(i)
                      c(rep.int(TRUE, Nout1[i]), rep.int(FALSE, Nout2[i]))))
  if (transf) {
    bfroutine <- "bfspmu"
    s1 <- mu[, lz, drop = FALSE]
    s2 <- mu[, !lz, drop = FALSE]
  } else {
    bfroutine <- "bfspz"
    s1 <- z[, lz, drop = FALSE]
    s2 <- z[, !lz, drop = FALSE]
  }
  logbf <- numeric(nruns)
  lglk1 <- matrix(0., Ntot1, nruns)
  lglk2 <- matrix(0., Ntot2, nruns)
  isweights <- numeric(Ntot2)
  zcv <- matrix(0., Ntot2, nruns)

  ## RUN estimation of Bayes factors
  tm2 <- system.time({
    RUN2 <- .Fortran(bfroutine,
                     isweights = isweights, zcv = zcv, logbf = logbf,
                     lglk1 = lglk1, lglk2 = lglk2,
                     as.double(phi), as.double(omg),
                     as.double(nu), as.double(s1), as.integer(Nout1),
                     as.integer(Ntot1), as.double(s2), as.integer(Nout2),
                     as.integer(Ntot2),
                     as.double(y), as.double(l), as.double(F),
                     as.double(dm), as.double(betm0), as.double(betQ0),
                     as.double(ssqdf), as.double(ssqsc), as.double(tsqdf),
                     as.double(tsq), as.double(kappa), as.integer(icf),
                     as.integer(n), as.integer(p), as.integer(nruns),
                     as.integer(ifam), as.integer(imeth))
  })

  if (verbose) {
    message ("Calculated Bayes factors and importance weights at skeleton
points: ", round(tm2[1]), " sec")
  }

  refbf <- RUN2$logbf[reference]
  logbf <- RUN2$logbf - refbf
  isweights <- RUN2$isweights
  zcv <- RUN2$zcv

  ## Check for non-finite values in control variates
  if (useCV && any(!is.finite(zcv))) {
    warning ("Computed non-finite control variates, probably due to numerical
overflow. Control variates corrections will not be used.")
    useCV <- FALSE
  }

  ## Optimisation
  imaxlogbf <- which.max(logbf)
  wbf <- exp(logbf - logbf[imaxlogbf] - log(sum(exp(logbf - logbf[imaxlogbf]))))
  pstart.d <- colSums(cbind(nu, phi, omg, kappa)*wbf)
  i <- is.na(pstart) & estim
  pstart[i] <- pmax(pmin(upper[i], pstart.d[i]), lower[i])
  ## Function to optimise
  if (useCV) {
    if (transf) {
      froutine <- "calcbmu_cv"
    } else {
      froutine <- "calcbz_cv"
    }
    fn <- function (par) {
      parin <- as.double(pstart)
      parin[estim] <- par
      RUN <- .Fortran(froutine, 0.0, parin[2], parin[1], parin[3], parin[4],
                      as.integer(icf), 1L, 1L, as.integer(Ntot2),
                      as.double(s2), as.double(isweights), as.double(zcv),
                      as.integer(n), as.integer(p), as.integer(nruns),
                      as.double(betm0),
                      as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                      as.double(tsqdf), as.double(tsq), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.integer(ifam))
      -RUN[[1]][1]
    }
  } else {
    if (transf) {
      froutine <- "calcbmu_st"
    } else {
      froutine <- "calcbz_st"
    }
    fn <- function (par) {
      parin <- as.double(pstart)
      parin[estim] <- par
      RUN <- .Fortran(froutine, 0.0, parin[2], parin[1], parin[3], parin[4],
                      as.integer(icf), 1L, 1L, as.integer(Ntot2), as.double(s2),
                      as.double(isweights), as.integer(n), as.integer(p),
                      as.double(betm0),
                      as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                      as.double(tsqdf), as.double(tsq), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.integer(ifam))
      -RUN[[1]][1]
    }
  }
  method <- if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  tm3 <- system.time({
    op <- stats::optim(pstart[estim], fn, method = method,
                       lower = lower[estim], upper = upper[estim],
                       control = control)
  })

  if (verbose) {
    message ("Finished optimization: ", round(tm3[1]), " sec")
  }

  parest <- as.double(pstart)
  parest[estim] <- op$par
  names(parest) <- parnmall

  ## Perform prediction
  Npro <- as.integer(Npro)
  if (Npro > 0) {
    Nprt <- as.integer(Nprt)
    Nprb <- as.integer(Nprb)
    lglks <- numeric(Npro)
    zs <- gmus <- matrix(0, n, Npro)
    zs[, 1] <- rowMeans(z[, seq(1 + c(0, Nout[-nruns])[imaxlogbf],
                                Nout[imaxlogbf])])
    z0s <- gmu0s <- matrix(0, n0, Npro)
    beta <- matrix(0, p, Npro)
    ssq <- numeric(Npro)
    acc <- 0L
    phis <- rep.int(parest[2], Npro)
    omgs <- rep.int(parest[3], Npro)
    tmppars <- rep.int(0, 4)
    tm4 <- system.time({
      RUN4 <- .Fortran("mcspsample", ll = lglks, z = zs, z0 = z0s,
                       mu = gmus, mu0 = gmu0s, 
                       beta = beta, ssq = ssq, as.double(phis), as.double(omgs),
                       acc = acc,
                       as.double(y), as.double(l), as.double(F), as.double(F0),
                       as.double(betm0), as.double(betQ0), as.double(ssqdf),
                       as.double(ssqsc),
                       as.double(tmppars), 0, as.double(tmppars), 0, parest[4],
                       as.integer(icf), parest[1], 
                       as.double(tsq), as.double(dm), as.double(dmdm0),
                       as.integer(Npro), as.integer(Nprb), as.integer(Nprt),
                       as.integer(n), as.integer(n0), as.integer(p),
                       as.integer(ifam))
    })
    if (verbose) {
      message ("Performed Gibbs sampling: ", round(tm4[1]), " sec")
    }
    ll <- RUN4$ll
    zz0 <- mm0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN4$z
    zz0[!ii, ] <- RUN4$z0
    mm0[ii, ] <- RUN4$mu
    mm0[!ii, ] <- RUN4$mu0
    beta <- RUN4$beta
    ssq <- RUN4$ssq
    acc_ratio <- RUN4$acc/Npro
    sample <- list(z = zz0, mu = mm0, beta = beta, ssq = ssq,
                   acc_ratio = acc_ratio, whichobs = ii)
  } else {
    sample <- NULL
    tm4 <- NULL
  }
  
  times <- rbind(sampling = tm1, importance = tm2, optimization = tm3,
                 MCMC = tm4)
  out <- list(parest = parest, skeleton = cbind(parskel[parnm], logbf = logbf),
              optim = op, mcmcsample = sample, sys_time = times)
  out
}


##' Empirical Bayes estimation for the spatial transformed Gaussian model.
##'
##' Runs the MCMC sampling, computes the importance weights, and
##' estimates the parameters.
##' @title Empirical Bayes estimation for the TGRFM
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations (see the example in
##' \code{\link{mcsglmm}} for how to do this using
##' \code{\link{stackdata}}). At the observed locations the response
##' is assumed to be a total of replicated measurements. The number of
##' replications is inputted using the argument \code{weights}. See
##' the Note for cases where overflow may occur.
##' @param data An optional data frame containing the variables in the
##' model.
##' @param weights An optional vector of weights. Number of replicated
##' samples for Gaussian and gamma, number of trials for binomial,
##' time length for Poisson.
##' @param subset An optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##' with the coordinates of the sampled locations.
##' @param parskel A data frame with the components "linkp", "phi",
##' "omg", and "kappa", corresponding to the link function, the
##' spatial range, the relative nugget, and the spatial smoothness
##' parameters. The latter can be omitted if not used in the
##' correlation function. Let k denote the number of rows. Then, k
##' different MCMC samples will be taken from the models with
##' parameters fixed at those values. For a square grid the output
##' from the function \code{\link[base]{expand.grid}} can be used
##' here.
##' @param paroptim A named list with the components "linkp", "phi",
##' "omg", "kappa". The latter can be omitted if not used in the
##' correlation function. Each component must be numeric with length
##' 1, 2, or 3 with elements in increasing order but for the binomial
##' family linkp is also allowed to be the character "logit" and
##' "probit". If its length is 1, then the corresponding parameter is
##' considered to be fixed at that value. If 2, then the two numbers
##' denote the lower and upper bounds for the optimisation of that
##' parameter (infinities are allowed). If 3, these correspond to
##' lower bound, starting value, upper bound for the estimation of
##' that parameter.
##' @param corrfcn Spatial correlation function. See
##' \code{\link{ebsglmm}} for details.
##' @param Nout A scalar or vector of size k. Number of MCMC samples
##' to take for each run of the MCMC algorithm for the estimation of
##' the Bayes factors. See argument \code{parskel}.
##' @param Nthin A scalar or vector of size k. The thinning of the
##' MCMC algorithm for the estimation of the Bayes factors.
##' @param Nbi A scalar or vector of size k. The burn-in of the MCMC
##' algorithm for the estimation of the Bayes factors.
##' @param Npro A scalar. The number of Gibbs samples to take for
##' estimation of the conjugate parameters and for prediction at the
##' unsampled locations while the other parameters are fixed at their
##' empirical Bayes estimates.
##' @param Nprt The thinning of the Gibbs algorithm for the estimation
##' of the conjugate parameters and for prediction.
##' @param Nprb The burn-in of the Gibbs algorithm for the estimation
##' of the conjugate parameters and for prediction.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##' matrix. Can be a scalar, vector or matrix. The first two imply a
##' diagonal with those elements. Set this to 0 to indicate a flat
##' improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##' partial sill parameter.
##' @param tsqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the measurement error parameter.
##' @param tsqsc Scale for the scaled inverse chi-square prior for the
##' measurement error parameter.
##' @param zstart Optional starting value for the MCMC for the GRF.
##' This can be either a scalar, a vector of size n where n is the
##' number of sampled locations, or a matrix with dimensions n by k
##' where k is the number of the skeleton points in \code{parskel}.
##' @param dispersion The fixed dispersion parameter.
##' @param bfsize1 A scalar or vector of the same length as \code{...}
##' with all integer values or all values in (0, 1]. How many samples
##' (or what proportion of the sample) to use for estimating the Bayes
##' factors at the first stage. The remaining sample will be used for
##' estimating the Bayes factors in the second stage. Setting it to 1
##' will perform only the first stage.
##' @param reference Which model goes in the denominator of the Bayes
##' factors.
##' @param bfmethod Which method to use to calculate the Bayes
##' factors: Reverse logistic or Meng-Wong.
##' @param transf Whether to use the transformed sample mu for the
##' computations. Otherwise it uses z.
##' @param useCV Whether to use control variates for finer
##' corrections.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param control A list of control parameters for the optimisation.
##' See \code{\link[stats]{optim}}.
##' @param verbose Whether to print messages when completing each
##' stage on screen.
##' @return A list with components
##' \itemize{
##' \item \code{parest} The parameter paroptims
##' \item \code{skeleton} The skeleton points used with the corresponding
##' logarithm of the Bayes factors at those points. 
##' \item \code{optim} The output from the \code{\link[stats]{optim}}
##' function. 
##' \item \code{mcmcsample} The MCMC samples for the remaining
##' parameters and the random field. These samples correspond to the
##' Gibbs and Metropolis-Hasting samples after fixing the parameters
##' estimated by empirical Bayes at their empirical Bayes estimates.
##' \item \code{sys_time} The time taken to complete the MCMC
##' sampling, calculation of the importance weights, the
##' optimization and the final MCMC sampling.
##' }
##' @note To avoid numerical overflow it is recommended to scale the
##' response variable by its geometric mean. This is strongly
##' recommended when the exponent of the transformation is 0 or less.
##' Numerical overflow may occur when samples which are simulated from
##' a smaller exponent are being used to compute the likelihood
##' corresponding to a larger exponent.
##' @examples \dontrun{
##' ### Load the data
##' data(rhizoctonia)
##' rhiz <- na.omit(rhizoctonia)
##' rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
##'                               # rhizoctonia disease
##' ### Define the model
##' corrf <- "spherical"
##' ssqdf <- 1
##' ssqsc <- 1
##' tsqdf <- 1
##' tsqsc <- 1
##' betm0 <- 0
##' betQ0 <- diag(.01, 2, 2)
##' ### Skeleton points
##' philist <- seq(120, 280, 40)
##' linkp <- 1
##' omglist <- seq(0, 2.5, .5)
##' parlist <- expand.grid(phi = philist, linkp = linkp, omg = omglist)
##' paroptim <- list(linkp = linkp, phi = c(100, 300), omg = c(0, 2))
##' ### MCMC sizes
##' Nout <- Npro <- 100
##' Nthin <- Nprt <- 1
##' Nbi <- Nprb <- 0
##' ### Estimate
##' est <- ebstrga(Yield ~ IR, rhiz, 
##'                atsample = ~ Xcoord + Ycoord, parskel = parlist,
##'                paroptim = paroptim, corrfcn = corrf, 
##'                Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                Npro = Npro, Nprt = Nprt, Nprb = Nprb, 
##'                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
##'                tsqdf = tsqdf, tsqsc = tsqsc,
##'                useCV=TRUE)
##' est$parest
##' }
##' @references Roy, V., Evangelou, E. and Zhu, Z. (2014). Empirical
##' Bayes methods for the transformed Gaussian random fields model
##' with additive measurement errors. In Upadhyay, S. K., Singh, U.,
##' Dey, D. K., and Loganathan, A., editors, \emph{Current Trends in
##' Bayesian Methodology with Applications}, Boca Raton, FL, USA, CRC
##' Press.
##' @importFrom sp spDists
##' @importFrom stats model.matrix model.response model.weights as.formula
##' @export 
ebstrga <- function (formula,
                     data, weights, subset, atsample, parskel, paroptim,
                     corrfcn = c("matern", "spherical", "powerexponential"), 
                     Nout, Nthin = 1, Nbi = 0, Npro, Nprt = 1, Nprb = 0, 
                     betm0, betQ0, ssqdf, ssqsc,
                     tsqdf, tsqsc, zstart, dispersion = 1,
                     bfsize1 = 0.8, reference = 1, bfmethod = c("RL", "MW"), 
                     transf = FALSE, useCV = TRUE, longlat = FALSE, 
                     control = list(), verbose = TRUE) {

  ## Family
  family <- "transformed-gaussian"
  ifam <- 0L

  ## Logical input
  useCV <- isTRUE(useCV)
  longlat <- isTRUE(longlat)
  verbose <- isTRUE(verbose)

  ## Correlation function
  corrfcn <- match.arg(corrfcn)
  icf <- match(corrfcn, eval(formals()$corrfcn))
  needkappa <- corrfcn %in% c("matern", "powerexponential")

  ## Design matrix and data
  if (missing(data)) data <- environment(formula)
  mfc <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"),
             names(mfc), 0L)
  mfc <- mfc[c(1L, m)]
  mfc$drop.unused.levels <- TRUE
  mfc$na.action <- "na.pass"
  mfc[[1L]] <- quote(stats::model.frame)
  mf <- eval(mfc, parent.frame())
  mt <- attr(mf, "terms")
  FF <- model.matrix(mt,mf)
  if (!all(is.finite(FF))) stop ("Non-finite values in the design matrix")
  p <- NCOL(FF)
  yy <- unclass(model.response(mf))
  if (!is.vector(yy)) {
    stop ("The response must be a vector")
  }
  yy <- as.double(yy)
  ll <- model.weights(mf)

  ## All locations
  locvars <- all.vars(atsample)
  formula1 <- as.formula(paste('~', paste(c(locvars, all.vars(formula)),
                                          collapse = ' + ')))
  mfc1 <- mfc
  mfc1$formula <- formula1
  mf1 <- eval(mfc1, parent.frame())
  m <- match(locvars, names(mf1))
  loc <- as.matrix(mf1[, m])
  if (!all(is.finite(loc))) stop ("Non-finite values in the locations")

  ## Check corrfcn with loc
  if (corrfcn == "spherical" & NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  n <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, n) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  F <- FF[ii, , drop = FALSE]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  n0 <- sum(!ii)
  if (n0 > 0) {
    F0 <- FF[!ii, , drop = FALSE]
    dmdm0 <- sp::spDists(loc[ii, , drop = FALSE], loc[!ii, , drop = FALSE],
                         longlat = longlat)
  } else {
    F0 <- dmdm0 <- numeric(0)
    dim(F0) <- c(0, p)
    dim(dmdm0) <- c(n, 0)
  }

  ## Priors
  if (all(is.finite(betQ0[upper.tri(betQ0)]))) {
    if (length(betQ0) == 1 && betQ0[1] == 0) {
      ## Uniform prior
      betQ0 <- matrix(0, p, p)
      betm0 <- rep(0, p)
    } else if (length(betQ0) == 1 || length(betQ0) == p) {
      if (any(betQ0 <= 0)) stop ('betQ0 not > 0')
      betQ0 <- diag(betQ0, p, p)
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      }
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(n + ssqdf)
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## Read and check parskel
  parskel <- .check_pargrid(parskel, family, corrfcn)
  nruns <- nrow(parskel)
  linkp <- parskel$linkp
  phi <- parskel$phi
  omg <- parskel$omg
  kappa <- parskel$kappa
  nu <- parskel$nu
  parnmall <- c("linkp", "phi", "omg", "kappa")
  parnm <- c("linkp", "phi", "omg", if(needkappa) "kappa")

  ## Method for computing the Bayes factors
  bfmethod <- match.arg(bfmethod)
  imeth <- match(bfmethod, eval(formals()$bfmethod))

  ## Read and check paroptim argument
  if (!is.list(paroptim)) {
    stop ("Argument paroptim must be a list")
  }
  if (!all(parnm %in% names(paroptim))) {
    stop (paste("Named components in the argument paroptim",
                "must be", paste(parnm, collapse = " ")))
  } else {
    paroptim <- paroptim[parnm]
  }
  if (!needkappa) paroptim$kappa <- 0
  lower <- upper <- pstart <- rep.int(NA, 4)
  estim <- rep.int(FALSE, 4)
  linkpe <- paroptim[["linkp"]]
  if (is.character(linkpe) | is.factor(linkpe)) {
    if (family == "binomial") {
      if (linkpe == "logit") {
        if (linkp == "logit") {
          pstart[1] <- -1
          estim[1] <- FALSE
        } else {
          stop ("Link in paroptim inconsistent with link in parskel")
        }
      } else if (linkpe == "probit") {
        if (linkp == "probit") {
          pstart[1] <- 0
          estim[1] <- FALSE
        } else {
          stop ("Link in paroptim inconsistent with link in parskel")
        }
      } else {
        stop ("Character link for the binomial family can be either
\"logit\" or \"probit\" in argument paroptim")
      }
    } else stop ("Character link in argument paroptim is only used for
the binomial family")
  } else if (is.numeric(linkpe)) {
    if (length(linkpe) < 1 | length(linkpe) > 3) {
      stop ("The length for a numeric component in paroptim must be 1, 2, or 3")
    } else if (length(linkpe) == 1) {
      pstart[1] <- linkpe
      estim[1] <- FALSE
    } else if (length(linkpe) == 2) {
      pstart[1] <- NA
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[2]
      if (lower[1] >= upper[1]) {
        stop ("The lower bound must be less than the upper bound for linkp
in paroptim")
      }
    } else {
      pstart[1] <- linkpe[2]
      estim[1] <- TRUE
      lower[1] <- linkpe[1]
      upper[1] <- linkpe[3]
      if (lower[1] > estim[1] | estim[1] > upper[1] | lower[1] == upper[1]) {
        stop ("The elements in the component linkp in paroptim must be ordered")
      }
    }
  } else {
    stop ("The element linkp in paroptim must be either numeric or character")
  }
  for (i in 2:4) {
    ppp <- paroptim[[parnmall[i]]]
    if (!is.numeric(ppp)) {
      stop(paste("The element", parnmall[i], "in paroptim must be numeric"))
    }
    lppp <- length(ppp)
    if (lppp < 1 | lppp > 3) {
      stop (paste("The element", parnmall[i], "in paroptim must have 1, 2, or 3
components"))
    }
    if (lppp == 1) {
      pstart[i] <- ppp
      estim[i] <- FALSE
    } else if (lppp == 2) {
      pstart[i] <- NA
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[2]
      if (lower[i] >= upper[i]) {
        stop (paste("The lower bound must be less than the upper bound for",
                     parnmall[i], "in paroptim"))
      }
    } else {
      pstart[i] <- ppp[2]
      estim[i] <- TRUE
      lower[i] <- ppp[1]
      upper[i] <- ppp[3]
      if (lower[i] > estim[i] | estim[i] > upper[i] | lower[i] == upper[i]) {
        stop (paste("The elements in the component", parnmall[i],
                     "in paroptim must be ordered"))
      }
    }
  }

  ## MCMC samples
  Nout <- rep(as.integer(Nout), length.out = nruns)
  Nbi <- rep(as.integer(Nbi), length.out = nruns)
  Nthin <- rep(as.integer(Nthin), length.out = nruns)
  Ntot <- sum(Nout)
  z <- matrix(0, n, Ntot)
  lglk <- numeric(Ntot)
  dispersion <- NA

  ## Starting value for z
  if (missing(zstart)) {
    zstart <- y/l
    zstart <- sapply(1:nruns,
                     function (i) linkfcn(zstart, linkp[i], "gaussian"))
    ## zstart <- pmax(zstart, -1e8) + rnorm(n*nruns, 0, sqrt(ssqsc))
  }
  z[, 1 + c(0, cumsum(Nout[-nruns]))] <- zstart


  bfsize1 <- as.double(bfsize1)
  if (length(bfsize1) > nruns) {
    warning ("The number of elements in bfsize1 exceeds the number of runs;
the extra elements will be discarded")
    bfsize1 <- bfsize1[1:nruns]
  }
  if (any(bfsize1 <= 0)) {
    stop ("Argument bfsize1 must be positive")
  } else if (all(bfsize1 <= 1)) {
    Nout1 <- as.integer(bfsize1*Nout)
    if (any(Nout1 == 0)) stop ("Calculated 0 sizes; give a larger bfsize1")
  } else if (all(bfsize1 >= 1)) {
    Nout1 <- as.integer(bfsize1)
    if (any(Nout1 > Nout)) stop ("The bfsize1 exceeds the number of samples")
  } else {
    stop ("Argument bfsize1 is a mix of proportions and sizes")
  }
  Ntot1 <- sum(Nout1)
  Nout2 <- Nout - Nout1
  Ntot2 <- sum(Nout2)

  tsq <- if (ifam > 0) dispersion else tsqsc

  ## RUN MCMC
  tm1 <- system.time({
    RUN1 <- .Fortran("samplemulti", lglk = lglk, z = z, mu = z,
                     as.double(phi), as.double(omg),
                     as.double(y), as.double(l),
                     as.double(F), as.double(betm0), as.double(betQ0),
                     as.double(ssqdf), as.double(ssqsc), as.double(kappa),
                     as.integer(icf), as.double(nu), as.double(tsqdf),
                     as.double(tsq),
                     as.double(dm), as.integer(Ntot), as.integer(Nout),
                     as.integer(Nbi), as.integer(Nthin), as.integer(n),
                     as.integer(p), as.integer(nruns), as.integer(ifam))
  })

  if (verbose) {
    message ("Completed MCMC sampling: ", round(tm1[1]), " sec")
  }

  ## Prepare data for first stage
  z <- RUN1$z
  mu <- RUN1$mu
  lz <- unlist(lapply(1:nruns, function(i)
                      c(rep.int(TRUE, Nout1[i]), rep.int(FALSE, Nout2[i]))))
  if (transf) {
    bfroutine <- "bfspmu"
    s1 <- mu[, lz, drop = FALSE]
    s2 <- mu[, !lz, drop = FALSE]
  } else {
    bfroutine <- "bfspz"
    s1 <- z[, lz, drop = FALSE]
    s2 <- z[, !lz, drop = FALSE]
  }
  logbf <- numeric(nruns)
  lglk1 <- matrix(0., Ntot1, nruns)
  lglk2 <- matrix(0., Ntot2, nruns)
  isweights <- numeric(Ntot2)
  zcv <- matrix(0., Ntot2, nruns)

  ## RUN estimation of Bayes factors
  tm2 <- system.time({
    RUN2 <- .Fortran(bfroutine,
                     isweights = isweights, zcv = zcv, logbf = logbf,
                     lglk1 = lglk1, lglk2 = lglk2,
                     as.double(phi), as.double(omg),
                     as.double(nu), as.double(s1), as.integer(Nout1),
                     as.integer(Ntot1), as.double(s2), as.integer(Nout2),
                     as.integer(Ntot2),
                     as.double(y), as.double(l), as.double(F), as.double(dm),
                     as.double(betm0), as.double(betQ0),
                     as.double(ssqdf), as.double(ssqsc), as.double(tsqdf),
                     as.double(tsq), as.double(kappa), as.integer(icf),
                     as.integer(n), as.integer(p), as.integer(nruns),
                     as.integer(ifam), as.integer(imeth))
  })

  if (verbose) {
    message ("Calculated Bayes factors and importance weights at skeleton
points: ", round(tm2[1]), " sec")
  }

  refbf <- RUN2$logbf[reference]
  logbf <- RUN2$logbf - refbf
  isweights <- RUN2$isweights
  zcv <- RUN2$zcv

  ## Check for non-finite values in control variates
  if (useCV && any(!is.finite(zcv))) {
    warning ("Computed non-finite control variates, probably due to numerical
overflow. Control variates corrections will not be used.")
    useCV <- FALSE
  }

  ## Optimisation
  imaxlogbf <- which.max(logbf)
  pmaxlogbf <- c(nu[imaxlogbf], phi[imaxlogbf], omg[imaxlogbf],
                 kappa[imaxlogbf])
  i <- is.na(pstart) & estim
  pstart[i] <- pmax(pmin(upper[i], pmaxlogbf[i]), lower[i])
  ## Function to optimise
  if (useCV) {
    if (transf) {
      froutine <- "calcbmu_cv"
    } else {
      froutine <- "calcbz_cv"
    }
    fn <- function (par) {
      parin <- as.double(pstart)
      parin[estim] <- par
      RUN <- .Fortran(froutine, 0.0, parin[2], parin[1], parin[3], parin[4],
                      as.integer(icf), 1L, 1L, as.integer(Ntot2),
                      as.double(s2), as.double(isweights), as.double(zcv),
                      as.integer(n), as.integer(p), as.integer(nruns),
                      as.double(betm0),
                      as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                      as.double(tsqdf), as.double(tsq), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.integer(ifam))
      -RUN[[1]][1]
    }
  } else {
    if (transf) {
      froutine <- "calcbmu_st"
    } else {
      froutine <- "calcbz_st"
    }
    fn <- function (par) {
      parin <- as.double(pstart)
      parin[estim] <- par
      RUN <- .Fortran(froutine, 0.0, parin[2], parin[1], parin[3], parin[4],
                      as.integer(icf), 1L, 1L, as.integer(Ntot2), as.double(s2),
                      as.double(isweights), as.integer(n), as.integer(p),
                      as.double(betm0),
                      as.double(betQ0), as.double(ssqdf), as.double(ssqsc),
                      as.double(tsqdf), as.double(tsq), as.double(y),
                      as.double(l), as.double(F), as.double(dm),
                      as.integer(ifam))
      -RUN[[1]][1]
    }
  }
  method <- if (sum(estim) == 1) "Brent" else "L-BFGS-B"
  tm3 <- system.time({
    op <- stats::optim(pstart[estim], fn, method = method,
                       lower = lower[estim], upper = upper[estim],
                       control = control)
  })

  if (verbose) {
    message ("Finished optimization: ", round(tm3[1]), " sec")
  }

  ## Output
  parest <- as.double(pstart)
  parest[estim] <- op$par
  names(parest) <- parnmall

  ## Perform prediction
  Npro <- as.integer(Npro)
  if (Npro > 0) {
    Nprt <- as.integer(Nprt)
    Nprb <- as.integer(Nprb)
    lglks <- numeric(Npro)
    zs <- gmus <- matrix(0, n, Npro)
    zs[, 1] <- rowMeans(z[, seq(1 + c(0, Nout[-nruns])[imaxlogbf],
                                Nout[imaxlogbf])])
    z0s <- gmu0s <- matrix(0, n0, Npro)
    beta <- matrix(0, p, Npro)
    ssq <- numeric(Npro)
    tsqs <- numeric(Npro)
    acc <- 0L
    phis <- rep.int(parest[2], Npro)
    omgs <- rep.int(parest[3], Npro)
    tmppars <- rep.int(0, 4)
    tm4 <- system.time({
      RUN4 <- .Fortran("trgasample", ll = lglks, z = zs, z0 = z0s,
                       mu = gmus, mu0 = gmu0s, 
                       beta = beta, ssq = ssq, tsq = tsqs,
                       as.double(phis), as.double(omgs), acc = acc,
                       as.double(y), as.double(l), as.double(F), as.double(F0),
                       as.double(betm0), as.double(betQ0), as.double(ssqdf),
                       as.double(ssqsc), as.double(tsqdf), as.double(tsqsc), 
                       as.double(tmppars), 0, as.double(tmppars), 0, parest[4],
                       as.integer(icf), parest[1],
                       as.double(dm), as.double(dmdm0), as.integer(Npro),
                       as.integer(Nprb), as.integer(Nprt), as.integer(n),
                       as.integer(n0), as.integer(p))
    })
    if (verbose) {
      message ("Performed Gibbs sampling: ", round(tm4[1]), " sec")
    }
    ll <- RUN4$ll
    zz0 <- mm0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN4$z
    zz0[!ii, ] <- RUN4$z0
    mm0[ii, ] <- RUN4$mu
    mm0[!ii, ] <- RUN4$mu0
    beta <- RUN4$beta
    ssq <- RUN4$ssq
    tsq <- RUN4$tsq
    acc_ratio <- RUN4$acc/Npro
    sample <- list(z = zz0, mu = mm0, beta = beta, ssq = ssq, tsq = tsq,
                   acc_ratio = acc_ratio,
                   whichobs = ii)
  } else {
    sample <- NULL
    tm4 <- NULL
  }
  
  times <- rbind(sampling = tm1, importance = tm2, optimization = tm3,
                 MCMC = tm4)
  out <- list(parest = parest, skeleton = cbind(parskel[parnm], logbf = logbf),
              optim = op, mcmcsample = sample, sys_time = times)
  out
}


## Check and return the parameter grid
.check_pargrid <- function(pargrid, family, corrfcn)
{
  pargrid <- data.frame(pargrid, stringsAsFactors = FALSE)
  ##if(!is.list(pargrid)) stop ("Argument pargrid must be a list")
  parnm <- c("linkp", "phi", "omg", "kappa")
  needkappa <- corrfcn %in% c("matern", "powerexponential")
  if (!needkappa) pargrid$kappa <- 0
  if (!all(parnm %in% names(pargrid))) {
    stop (paste("Argument pargrid must have the names",
                paste(c("linkp", "phi", "omg", if (needkappa) "kappa"),
                      collapse=" ")))
  } else {
    pargrid <- pargrid[parnm]
  }
  nruns <- nrow(pargrid)
  phi <- pargrid[["phi"]]
  omg <- pargrid[["omg"]]
  kappa <- pargrid[["kappa"]]
  linkp <- pargrid[["linkp"]]
  if (any (phi < 0)) stop ("Element phi given must be non-negative")
  if (any (omg < 0)) stop ("Element omg given must be non-negative")

  ## Check kappa and corrfcn
  kappa <- as.double(kappa)
  if (any(kappa < 0) & corrfcn %in% c("matern", "powerexponential")) {
    stop ("Argument kappa cannot be negative")
  }
  if (any(kappa > 2) & corrfcn == "powerexponential") {
    stop ("Argument kappa cannot be more than 2")
  }

  ## Check if linkp conforms with family
  if (is.character(linkp) | is.factor(linkp)) {
    if (family == "binomial") {
      if (all(linkp == "logit")) {
        nu <- rep.int(-1, nruns)
      } else if (all(linkp == "probit")) {
        nu <- rep.int(0, nruns)
      } else stop ("Cannot recognise character link for binomial")
    } else stop ("Character link is only allowed for binomial")
  } else if (!is.numeric(linkp)) {
    stop ("Element linkp in list must be numeric, or in the case of
the binomial can also be the character \"logit\" or \"probit\"")
  } else {
    nu <- as.double(linkp)
    if ((family %in% c("binomial", "Wallace.binomial")) && any(nu <= 0)) {
      stop ("The link parameter must be positive")
    }
  }

  ## Output
  data.frame(linkp = linkp, phi = phi, omg = omg, kappa = kappa, nu = nu)
}
