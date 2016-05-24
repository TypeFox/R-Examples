##' Draw MCMC samples from the Spatial GLMM with known link function
##'
##' The four-parameter prior for \code{phi} is defined by
##' \deqn{\propto (\phi - \theta_4)^{\theta_2 -1} \exp\{-(\frac{\phi -
##' \theta_4}{\theta_1})^{\theta_3}\}}{propto (phi -
##' phiprior[4])^(phiprior[2]-1) *
##' exp(-((phi-phiprior[4])/phiprior[1])^phiprior[3])} for \eqn{\phi >
##' \theta_4}{phi > phiprior[4]}. The prior for \code{omg} is similar.
##' The prior parameters correspond to scale, shape, exponent, and
##' location. See \code{arXiv:1005.3274} for details of this
##' distribution.
##'
##' The GEV (Generalised Extreme Value) link is defined by \deqn{\mu =
##' 1 - \exp\{-\max(0, 1 + \nu x)^{\frac{1}{\nu}}\}}{mu = 1 -
##' \exp[-max(0, 1 + nu x)^(1/nu)]} for any real \eqn{\nu}{nu}. At
##' \eqn{\nu = 0}{nu = 0} it reduces to the complementary log-log
##' link.
##' @title MCMC samples from the Spatial GLMM
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations (see the examples on how to do this
##' using the function \code{\link{stackdata}}). At the observed
##' locations the response is assumed to be a total of replicated
##' measurements. The number of replications is inputted using the
##' argument \code{weights}.
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
##' @param Nout Number of MCMC samples to return.
##' @param Nthin The thinning of the MCMC algorithm.
##' @param Nbi The burn-in of the MCMC algorithm.
##' @param betm0 Prior mean for beta (a vector or scalar).
##' @param betQ0 Prior standardised precision (inverse variance)
##' matrix. Can be a scalar, vector or matrix. The first two imply a
##' diagonal with those elements. Set this to 0 to indicate a flat
##' improper prior.
##' @param ssqdf Degrees of freedom for the scaled inverse chi-square
##' prior for the partial sill parameter.
##' @param ssqsc Scale for the scaled inverse chi-square prior for the
##' partial sill parameter.
##' @param phipars Parameters for the generalized inverse gamma prior
##' for the spatial range parameter \code{phi}. A four dimensional
##' vector with parameters scale, shape, exponent, location in that
##' order. See Details.
##' @param omgpars Parameters for the generalized inverse gamma prior
##' for the relative nugget parameter \code{omg}. A four dimensional
##' vector with parameters scale, shape, exponent, location in that
##' order. See Details.
##' @param corrfcn Spatial correlation function. See
##' \code{\link{ebsglmm}} for details.
##' @param kappa Spatial correlation parameter. Smoothness parameter
##' for Matern, exponent for the power family.
##' @param linkp Parameter of the link function. For binomial, a
##' positive number for the degrees of freedom of the robit family or
##' "logit" or "probit". For the other families any number for the
##' exponent of the Box-Cox transformation.
##' @param phisc Random walk parameter for \code{phi}. Smaller values
##' increase the acceptance ratio. Set this to 0 for fixed \code{phi}.
##' In this case the fixed value is given in the argument
##' \code{phistart}.
##' @param omgsc Random walk parameter for \code{omg}. Smaller values
##' increase the acceptance ratio. Set this to 0 for fixed \code{omg}.
##' In this case the fixed value is given in the argument
##' \code{omgstart}.
##' @param zstart Optional starting value for the MCMC for the GRF.
##' This can be either a scalar, a vector of size n where n is the
##' number of sampled locations.
##' @param phistart Optional starting value for the MCMC for the
##' spatial range parameter \code{phi}. Defaults to the mean of its
##' prior. If \code{phisc} is 0, then this argument is required and it
##' corresponds to the fixed value of \code{phi}.
##' @param omgstart Optional starting value for the MCMC for the relative
##' nugget parameter \code{omg}. Defaults to the mean of its prior. If
##' \code{omgsc} is 0, then this argument is required and
##' itcorresponds to the fixed value of \code{omg}.
##' @param dispersion The fixed dispersion parameter.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param test Whether this is a trial run to monitor the acceptance
##' ratio of the random walk for \code{phi} and \code{omg}. If set to
##' \code{TRUE}, the acceptance ratio will be printed on the screen
##' every 100 iterations of the MCMC. Tune the \code{phisc} and
##' \code{omgsc} parameters in order to achive 20 to 30\% acceptance.
##' Set this to a positive number to change the default 100. No
##' thinning or burn-in are done when testing.
##' @return A list containing the MCMC samples and other variables as
##' follows:
##' \itemize{
##'  \item \code{z} A matrix containing the MCMC samples for the
##' spatial random field. Each column is one sample. 
##'  \item \code{mu} A matrix containing the MCMC samples for the
##' mean response (a transformation of z). Each column is one sample. 
##'  \item \code{beta} A matrix containing the MCMC samples for the
##' regressor coefficients. Each column is one sample. 
##'  \item \code{ssq} A vector with the MCMC samples for the partial
##' sill parameter. 
##'  \item \code{phi} A vector with the MCMC samples for the spatial
##' range parameter. 
##'  \item \code{omg} A vector with the MCMC samples for the relative
##' nugget parameter. 
##'  \item \code{nu} The link function parameter translated to
##' numeric code used internally. 
##'  \item \code{logLik} A vector containing the value of the
##' log-likelihood evaluated at each sample. 
##'  \item \code{acc_ratio} The acceptance ratio for the joint update
##' of the parameters \code{phi} and \code{omg}. 
##'  \item \code{sys_time} The total computing time for the MCMC sampling.
##'  \item \code{Nout}, \code{Nbi},  \code{Nthin} As in input. Used
##' internally in other functions. 
##'  \item \code{response} The value of the response variable at the
##' observed locations. Used internally in other functions. 
##'  \item \code{weights} The response weights at the observed
##' locations. Used internally in other functions. 
##'  \item \code{modelmatrix} The model matrix at the observed
##' locations. Used internally in other functions. 
##'  \item \code{family} As in input. Used internally in other functions.
##'  \item \code{betm0}, \code{betQ0}, \code{ssqdf}, \code{ssqsc},
##' \code{corrfcn}, \code{kappa}, \code{dispersion} As in
##' input. Used internally in other functions. 
##'  \item \code{locations} Coordinates of the observed locations.
##' Used internally in other functions. 
##'  \item \code{whichobs} A logical vector indicated which rows in
##' the data and in the MCMC samples for the spatial random field
##' correspond to the observed locations.
##' }
##' @examples \dontrun{
##' data(rhizoctonia)
##'  
##' ### Create prediction grid
##' predgrid <- mkpredgrid2d(rhizoctonia[c("Xcoord", "Ycoord")],
##'                          par.x = 100, chull = TRUE, exf = 1.2)
##'  
##' ### Combine observed and prediction locations
##' rhizdata <- stackdata(rhizoctonia, predgrid$grid)
##' 
##' ### Define the model
##' corrf <- "spherical"
##' kappa <- 0
##' ssqdf <- 1
##' ssqsc <- 1
##' betm0 <- 0
##' betQ0 <- .01
##' phiprior <- c(100, 1, 1000, 100) # U(100, 200)
##' phisc <- 3
##' omgprior <- c(2, 1, 1, 0)        # Exp(mean = 2)
##' omgsc <- .1
##' linkp <- "probit"
##' 
##' ### MCMC sizes
##' Nout <- 100
##' Nthin <- 1
##' Nbi <- 0
##' 
##' ### Trial run
##' emt <- mcsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
##'                atsample = ~ Xcoord + Ycoord,
##'                Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
##'                phipars = phiprior, omgpars = omgprior, linkp = linkp, 
##'                corrfcn = corrf, kappa = kappa, phisc = phisc, omgsc = omgsc, 
##'                dispersion = 1, test = 10)
##' 
##' ### Full run
##' emc <- mcsglmm(Infected ~ 1, 'binomial', rhizdata, weights = Total,
##'                atsample = ~ Xcoord + Ycoord,
##'                Nout = Nout, Nthin = Nthin, Nbi = Nbi,
##'                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
##'                phipars = phiprior, omgpars = omgprior, linkp = linkp, 
##'                corrfcn = corrf, kappa = kappa, phisc = phisc, omgsc = omgsc, 
##'                dispersion = 1, test = FALSE)
##' 
##' 
##' plot.ts(cbind(phi = emc$phi, omg = emc$omg, beta = c(emc$beta),
##'               ssq = emc$ssq), nc = 2)
##' 
##' emcmc <- mcmcmake(emc)
##' summary(emcmc[, c("phi", "omg", "beta", "ssq")])
##' }
##' @importFrom sp spDists
##' @importFrom stats model.matrix model.response model.weights as.formula
##' @export 
mcsglmm <- function (formula,
                     family = c("gaussian", "binomial", "poisson", "Gamma",
                       "GEV.binomial", "GEVD.binomial", "Wallace.binomial"),
                     data, weights, subset, atsample,
                     Nout, Nthin = 1, Nbi = 0, betm0, betQ0, ssqdf, ssqsc,
                     phipars, omgpars,
                     corrfcn = c("matern", "spherical", "powerexponential"), 
                     kappa, linkp, phisc, omgsc,
                     zstart, phistart, omgstart,
                     dispersion = 1, longlat = FALSE, test = FALSE) {
  ## Family
  family <- match.arg(family)
  ifam <- match(family, eval(formals()$family))

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

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  k <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, k) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  if (family %in% c("binomial", "GEV.binomial", "GEVD.binomial",
                    "Wallace.binomial")) {
    l <- l - y # Number of failures
  }
  F <- FF[ii, , drop = FALSE]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  k0 <- sum(!ii)
  if (k0 > 0) {
    F0 <- FF[!ii, , drop = FALSE]
    dmdm0 <- sp::spDists(loc[ii, , drop = FALSE], loc[!ii, , drop = FALSE],
                         longlat = longlat)
  } else {
    F0 <- dmdm0 <- numeric(0)
    dim(F0) <- c(0, p)
    dim(dmdm0) <- c(k, 0)
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
      modeldf <- as.double(k + ssqdf)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      } 
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(k + ssqdf)
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")

  ## phi and omg
  phisc <- as.double(phisc)
  if (phisc < 0) stop ("Argument phisc must be non-negative")
  omgsc <- as.double(omgsc)
  if (omgsc < 0) stop ("Argument omgsc must be non-negative")
  if (phisc > 0) {
    if (missing(phipars)) stop ("Argument phipars not provided")
    phipars <- as.double(phipars)
    if (length(phipars) != 4) {
      stop ("Argument phipars must be a vector of length 4")
    }
    if (!all(is.finite(phipars))) stop ("Non-finite values in phipars")
    if (phipars[1] <= 0 || phipars[4] < 0 || phipars[2]*phipars[3] <= 0) {
      stop ("Invalid values in phipars")
    }
  } else {
    phipars <- rep.int(0, 4)
  }
  if (omgsc > 0) {
    if (missing(omgpars)) stop ("Argument omgpars not provided")
    omgpars <- as.double(omgpars)
    if (length(omgpars) != 4) {
      stop ("Argument omgpars must be a vector of length 4")
    }
    if (!all(is.finite(omgpars))) stop ("Non-finite values in omgpars")
    if (any(omgpars[1:3] <= 0) || omgpars[4] < 0) {
      stop ("Invalid values in omgpars")
    }
  } else {
    omgpars <- rep.int(0, 4)
  }

  ## Other fixed parameters
  kappa <- if (needkappa) as.double(kappa) else 0
  if (kappa < 0 && corrfcn %in% c("matern", "powerexponential")) {
    stop ("Argument kappa cannot be negative")
  }
  if (kappa > 2 && corrfcn == "powerexponential") {
    stop ("Argument kappa cannot be more than 2")
  }
  if (corrfcn == "spherical" && NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }
  dispersion <- as.double(dispersion)
  if (dispersion <= 0) stop ("Invalid argument dispersion")

  if (is.character(linkp) || is.factor(linkp)) {
    linkp <- as.character(linkp)
    if (family == "binomial") {
      if (all(linkp == "logit")) {
        nu <- -1
      } else if (all(linkp == "probit")) {
        nu <- 0
      } else {
        nu <- as.numeric(linkp)
        if (any(is.na(nu))) {
          stop ("Cannot recognise character link for binomial")
        } else if (any(nu <= 0)) {
          stop ("The robit link parameter must be positive")
        }
      }
    } else stop ("Character link is only allowed for binomial")
  } else if (!is.numeric(linkp)) {
    stop ("Element linkp in parameters must be numeric, or in the case of
the binomial can also be the character \"logit\" or \"probit\"")
  } else {
    nu <- as.double(linkp)
    if ((family %in% c("binomial", "Wallace.binomial")) && any(nu <= 0)) {
      stop ("The link parameter must be positive")
    }
  }

  ## MCMC samples
  Nout <- as.integer(Nout)
  Nbi <- as.integer(Nbi)
  Nthin <- as.integer(Nthin)
  lglk <- numeric(Nout)
  z <- matrix(0, k, Nout)
  z0 <- matrix(0, k0, Nout)
  beta <- matrix(0, p, Nout)
  ssq <- numeric(Nout)
  phi <- numeric(Nout)
  omg <- numeric(Nout)
  acc <- 0L

  ## Starting values
  if (missing(zstart)) {
    zstart <- switch(family,
                     binomial =, Wallace.binomial =,
                     GEV.binomial = (y+.5)/(y+l+1),
                     GEVD.binomial = (l+.5)/(y+l+1),
                     poisson = (y+.5)/(l+1), Gamma =, gaussian = y/l)
    zstart <- linkfcn(zstart, linkp, family)
    ## zstart <- pmax(zstart, -1e8) + rnorm(k, 0, sqrt(ssqsc))
  }
  z[, 1] <- zstart
  if (missing(phistart)) {
    if (phisc == 0) {
      stop ("Argument phistart needed for fixed phi")
    } else {
      if(phipars[2] == -1) {
        tmp <- .1/abs(phipars[3])
      } else {
        tmp <- abs((phipars[2]+1)/phipars[3])
      }
      phistart <- phipars[4] + phipars[1]*gamma(tmp)/
        gamma(phipars[2]/phipars[3])
    }
  } else {
    phistart <- as.double(phistart)
    if (phisc > 0 && phistart <= phipars[4]) {
      stop ("Starting value for phi not in the support of its prior")
    }
  }
  phi[1] <- phistart
  if (missing(omgstart)) {
    if (omgsc == 0) {
      stop ("Argument omgstart needed for fixed omg")
    } else {
      if(omgpars[2] == -1) {
        tmp <- .1/abs(omgpars[3])
      } else {
        tmp <- abs((omgpars[2]+1)/omgpars[3])
      }
      omgstart <- omgpars[4] + omgpars[1]*gamma(tmp)/
        gamma(omgpars[2]/omgpars[3])
    }
  } else {
    omgstart <- as.double(omgstart)
    if (omgsc > 0 && omgstart <= omgpars[4]) {
      stop ("Starting value for omg not in the support of its prior")
    }
  }
  omg[1] <- omgstart

  ## Run code
  if (test > 0) { # Running a test
    if (is.logical(test)) test <- 100
    test <- as.integer(test)
    tm <- system.time({
      RUN <- .Fortran("mcspsamtry", ll = lglk, z = z, phi = phi, omg = omg,
                      acc = acc,
                      as.double(y), as.double(l), as.double(F), 
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(phipars), as.double(phisc),
                      as.double(omgpars),
                      as.double(omgsc), as.double(kappa), as.integer(icf), 
                      as.double(nu), as.double(dispersion), as.double(dm),
                      as.integer(Nout), as.integer(test), as.integer(k),
                      as.integer(p), as.integer(ifam))
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    beta <- NULL
    ssq <- NULL
    phi <- RUN$phi
    attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
    attr(omg, 'fixed') <- omgsc == 0
    attr(nu, 'fixed') <- TRUE
    acc_ratio <- RUN$acc/Nout
    out <- list(z = zz0, beta = beta, ssq = ssq, phi = phi, omg = omg, nu = nu, 
                logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
                Nout = Nout, Nbi = Nbi, Nthin = Nthin,
                response = y, weights = l, modelmatrix = F, family = family,
                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
                corrfcn = corrfcn, kappa = kappa, 
                dispersion = dispersion, locations = loc[ii, , drop = FALSE],
                whichobs = ii)
  } else {
    tm <- system.time({
      RUN <- .Fortran("mcspsample", ll = lglk, z = z, z0 = z0,
                      mu = z, mu0 = z0, 
                      beta = beta, ssq = ssq,
                      phi = phi, omg = omg, acc = acc,
                      as.double(y), as.double(l), as.double(F), as.double(F0),
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(phipars), as.double(phisc),
                      as.double(omgpars),
                      as.double(omgsc), as.double(kappa), as.integer(icf), 
                      as.double(nu), as.double(dispersion), as.double(dm),
                      as.double(dmdm0), as.integer(Nout), as.integer(Nbi),
                      as.integer(Nthin), as.integer(k), as.integer(k0),
                      as.integer(p), as.integer(ifam))
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- mm0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    zz0[!ii, ] <- RUN$z0
    mm0[ii, ] <- RUN$mu
    mm0[!ii, ] <- RUN$mu0
    beta <- RUN$beta
    ssq <- RUN$ssq
    phi <- RUN$phi
    attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
    attr(omg, 'fixed') <- omgsc == 0
    attr(nu, 'fixed') <- TRUE
    acc_ratio <- RUN$acc/(Nout*Nthin + max(Nthin, Nbi))
    out <- list(z = zz0, mu = mm0,
                beta = beta, ssq = ssq, phi = phi, omg = omg, nu = nu, 
                logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
                Nout = Nout, Nbi = Nbi, Nthin = Nthin,
                response = y, weights = l, modelmatrix = F, family = family,
                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
                corrfcn = corrfcn, kappa = kappa, 
                dispersion = dispersion, locations = loc[ii, , drop = FALSE],
                whichobs = ii)
  }
  class(out) <- "geomcmc"
  out
}


##' Draw MCMC samples from the transformed Gaussian model with known
##' link function
##'
##' Simulates from the posterior distribution of this model.
##' @title MCMC samples from the transformed Gaussian model
##' @param formula A representation of the model in the form
##' \code{response ~ terms}. The response must be set to \code{NA}'s
##' at the prediction locations (see the example in
##' \code{\link{mcsglmm}} for how to do this using
##' \code{\link{stackdata}}). At the observed locations the response
##' is assumed to be a total of replicated measurements. The number of
##' replications is inputted using the argument \code{weights}.
##' @param data An optional data frame containing the variables in the
##' model.
##' @param weights An optional vector of weights. Number of replicated
##' samples.
##' @param subset An optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param atsample A formula in the form \code{~ x1 + x2 + ... + xd}
##' with the coordinates of the sampled locations.
##' @param Nout Number of MCMC samples to return.
##' @param Nthin The thinning of the MCMC algorithm.
##' @param Nbi The burn-in of the MCMC algorithm.
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
##' @param phipars Parameters for the generalized inverse gamma prior
##' for the range parameter \code{phi}. A four dimensional vector with
##' parameters scale, shape, exponent, location in that order. See
##' \code{\link{mcsglmm}}. 
##' @param omgpars Parameters for the generalized inverse gamma prior
##' for the relative nugget parameter \code{omg}. A four dimensional
##' vector with parameters scale, shape, exponent, location in that
##' order. See \code{\link{mcsglmm}}.
##' @param corrfcn Spatial correlation function. See
##' \code{\link{ebsglmm}} for details.
##' @param kappa Spatial correlation parameter. Smoothness parameter
##' for Matern, exponent for the power family.
##' @param linkp The exponent of the Box-Cox transformation.
##' @param phisc Random walk parameter for \code{phi}. Smaller values
##' increase the acceptance ratio. Set this to 0 for fixed \code{phi}.
##' In this case the fixed value is given in the argument
##' \code{phistart}.
##' @param omgsc Random walk parameter for \code{omg}. Smaller values
##' increase the acceptance ratio. Set this to 0 for fixed \code{omg}.
##' In this case the fixed value is given in the argument
##' \code{omgstart}.
##' @param zstart Optional starting value for the MCMC for the GRF.
##' This can be either a scalar, a vector of size n where n is the
##' number of sampled locations.
##' @param phistart Optional starting value for the MCMC for the
##' spatial range parameter \code{phi}. Defaults to the mean of its
##' prior. If \code{phisc} is 0, then this argument is required and it
##' corresponds to the fixed value of \code{phi}.
##' @param omgstart Optional starting value for the MCMC for the relative
##' nugget parameter \code{omg}. Defaults to the mean of its prior. If
##' \code{omgsc} is 0, then this argument is required and
##' itcorresponds to the fixed value of \code{omg}.
##' @param longlat How to compute the distance between locations. If
##' \code{FALSE}, Euclidean distance, if \code{TRUE} Great Circle
##' distance. See \code{\link[sp]{spDists}}.
##' @param test Whether this is a trial run to monitor the acceptance
##' ratio of the random walk for \code{phi} and \code{omg}. If set to
##' \code{TRUE}, the acceptance ratio will be printed on the screen
##' every 100 iterations of the MCMC. Tune the \code{phisc} and
##' \code{omgsc} parameters in order to achive 20 to 30\% acceptance.
##' Set this to a positive number to change the default 100. No
##' thinning or burn-in are done when testing.
##' @return A list containing the MCMC samples and other variables as
##' follows:
##' \itemize{
##'  \item \code{z} A matrix containing the MCMC samples for the
##' spatial random field. Each column is one sample. 
##'  \item \code{mu} A matrix containing the MCMC samples for the
##' mean response (a transformation of z). Each column is one sample. 
##'  \item \code{beta} A matrix containing the MCMC samples for the
##' regressor coefficients. Each column is one sample. 
##'  \item \code{ssq} A vector with the MCMC samples for the partial
## sill parameter. 
##'  \item \code{tsq} A vector with the MCMC samples for the
##' measurement error variance. 
##'  \item \code{phi} A vector with the MCMC samples for the spatial
##' range parameter. 
##'  \item \code{omg} A vector with the MCMC samples for the relative
##' nugget parameter. 
##'  \item \code{nu} The link function parameter translated to
##' numeric code used internally. 
##'  \item \code{logLik} A vector containing the value of the
##' log-likelihood evaluated at each sample. 
##'  \item \code{acc_ratio} The acceptance ratio for the joint update
##' of the parameters \code{phi} and \code{omg}. 
##'  \item \code{sys_time} The total computing time for the MCMC sampling.
##'  \item \code{Nout}, \code{Nbi},  \code{Nthin} As in input. Used
##' internally in other functions. 
##'  \item \code{response} The average of the response variable at the
##' observed locations, i.e. its value divided by the corresponding
##' weight. Used internally in other functions.  
##'  \item \code{weights} The response weights at the observed
##' locations. Used internally in other functions. 
##'  \item \code{modelmatrix} The model matrix at the observed
##' locations. Used internally in other functions. 
##'  \item \code{family} As in input. Used internally in other functions.
##'  \item \code{betm0}, \code{betQ0}, \code{ssqdf}, \code{ssqsc},
##' \code{corrfcn}, \code{kappa}, \code{tsqdf}, \code{tsqsc} As in
##' input. Used internally in other functions. 
##'  \item \code{locations} Coordinates of the observed locations.
##' Used internally in other functions. 
##'  \item \code{whichobs} A logical vector indicated which rows in
##' the data and in the MCMC samples for the spatial random field
##' correspond to the observed locations.
##' }
##' @examples \dontrun{
##' ### Load the data
##' data(rhizoctonia)
##' rhiz <- na.omit(rhizoctonia)
##' rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
##'                               # rhizoctonia disease
##' 
##' ### Define the model
##' corrf <- "spherical"
##' ssqdf <- 1
##' ssqsc <- 1
##' tsqdf <- 1
##' tsqsc <- 1
##' betm0 <- 0
##' betQ0 <- diag(.01, 2, 2)
##' phiprior <- c(200, 1, 1000, 100) # U(100, 300)
##' phisc <- 1
##' omgprior <- c(3, 1, 1000, 0) # U(0, 3)
##' omgsc <- 1.3
##' linkp <- 1
##' 
##' ## MCMC parameters
##' Nout <- 100
##' Nbi <- 0
##' Nthin <- 1
##' 
##' samplt <- mcstrga(Yield ~ IR, data = rhiz, 
##'                   atsample = ~ Xcoord + Ycoord, corrf = corrf, 
##'                   Nout = Nout, Nthin = Nthin,
##'                   Nbi = Nbi, betm0 = betm0, betQ0 = betQ0,
##'                   ssqdf = ssqdf, ssqsc = ssqsc,
##'                   tsqdf = tsqdf, tsqsc = tsqsc,
##'                   phipars = phiprior, omgpars = omgprior,
##'                   linkp = linkp,
##'                   phisc = phisc, omgsc = omgsc, test=10)
##' 
##' sample <- mcstrga(Yield ~ IR, data = rhiz, 
##'                   atsample = ~ Xcoord + Ycoord, corrf = corrf, 
##'                   Nout = Nout, Nthin = Nthin,
##'                   Nbi = Nbi, betm0 = betm0, betQ0 = betQ0,
##'                   ssqdf = ssqdf, ssqsc = ssqsc,
##'                   tsqdf = tsqdf, tsqsc = tsqsc,
##'                   phipars = phiprior, omgpars = omgprior,
##'                   linkp = linkp,
##'                   phisc = phisc, omgsc = omgsc, test=FALSE)
##' }
##' @importFrom sp spDists
##' @importFrom stats model.matrix model.response model.weights as.formula
##' @export 
mcstrga <- function (formula,
                     data, weights, subset, atsample,
                     Nout, Nthin = 1, Nbi = 0, betm0, betQ0, ssqdf, ssqsc,
                     tsqdf, tsqsc, phipars, omgpars,
                     corrfcn = c("matern", "spherical", "powerexponential"), 
                     kappa, linkp, phisc, omgsc,
                     zstart, phistart, omgstart, longlat = FALSE,
                     test = FALSE) {

  family <- "transformed-gaussian"

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

  ## Split sample, prediction
  ii <- is.finite(yy)
  y <- yy[ii]
  k <- sum(ii)
  l <- ll[ii]
  l <- if (is.null(l)) rep.int(1.0, k) else as.double(l)
  if (any(!is.finite(l))) stop ("Non-finite values in the weights")
  if (any(l <= 0)) stop ("Non-positive weights not allowed")
  ybar <- y/l
  F <- FF[ii, , drop = FALSE]
  dm <- sp::spDists(loc[ii, , drop = FALSE], longlat = longlat)
  k0 <- sum(!ii)
  if (k0 > 0) {
    F0 <- FF[!ii, , drop = FALSE]
    dmdm0 <- sp::spDists(loc[ii, , drop = FALSE], loc[!ii, , drop = FALSE],
                         longlat = longlat)
  } else {
    F0 <- dmdm0 <- numeric(0)
    dim(F0) <- c(0, p)
    dim(dmdm0) <- c(k, 0)
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
      modeldf <- as.double(k + ssqdf)
    } else if (length(betQ0) == p*p) {
      betQ0 <- matrix(as.double(betQ0), p, p)
      betQ0[lower.tri(betQ0)] <- 0
      betQ0eig <- eigen(t(betQ0), 1, 1)$values
      if (any (betQ0eig < sqrt(.Machine$double.eps))) {
        stop ('betQ0 not > 0 within tolerance')
      }
      betm0 <- rep(as.double(betm0), length.out = p)
      modeldf <- as.double(k + ssqdf)
    } else stop ('Bad betQ0')
  } else stop ('Non-finite betQ0')
  ssqdf <- as.double(ssqdf)
  if (ssqdf <= 0) stop ("Argument ssqdf must > 0")
  ssqsc <- as.double(ssqsc)
  if (ssqsc <= 0) stop ("Argument ssqsc must > 0")
  tsqdf <- as.double(tsqdf)
  if (tsqdf <= 0) stop ("Argument tsqdf must > 0")
  tsqsc <- as.double(tsqsc)
  if (tsqsc <= 0) stop ("Argument tsqsc must > 0")

  ## phi and omg
  phisc <- as.double(phisc)
  if (phisc < 0) stop ("Argument phisc must be non-negative")
  omgsc <- as.double(omgsc)
  if (omgsc < 0) stop ("Argument omgsc must be non-negative")
  if (phisc > 0) {
    if (missing(phipars)) stop ("Argument phipars not provided")
    phipars <- as.double(phipars)
    if (length(phipars) != 4) {
      stop ("Argument phipars must be a vector of length 4")
    }
    if (!all(is.finite(phipars))) stop ("Non-finite values in phipars")
    if (phipars[1] <= 0 || phipars[4] < 0 || phipars[2]*phipars[3] <= 0) {
      stop ("Invalid values in phipars")
    }
  } else {
    phipars <- rep.int(0, 4)
  }
  if (omgsc > 0) {
    if (missing(omgpars)) stop ("Argument omgpars not provided")
    omgpars <- as.double(omgpars)
    if (length(omgpars) != 4) {
      stop ("Argument omgpars must be a vector of length 4")
    }
    if (!all(is.finite(omgpars))) stop ("Non-finite values in omgpars")
    if (any(omgpars[1:3] <= 0) || omgpars[4] < 0) {
      stop ("Invalid values in omgpars")
    }
  } else {
    omgpars <- rep.int(0, 4)
  }

  ## Other fixed parameters
  kappa <- if (needkappa) as.double(kappa) else 0
  if (kappa < 0 && corrfcn %in% c("matern", "powerexponential")) {
    stop ("Argument kappa cannot be negative")
  }
  if (kappa > 2 && corrfcn == "powerexponential") {
    stop ("Argument kappa cannot be more than 2")
  }
  if (corrfcn == "spherical" && NCOL(loc) > 3) {
    stop ("Cannot use the spherical correlation for dimensions
grater than 3.")
  }

  nu <- as.double(linkp)
  
  ## MCMC samples
  Nout <- as.integer(Nout)
  Nbi <- as.integer(Nbi)
  Nthin <- as.integer(Nthin)
  z <- matrix(0, k, Nout)
  z0 <- matrix(0, k0, Nout)
  beta <- matrix(0, p, Nout)
  ssq <- numeric(Nout)
  tsq <- numeric(Nout)
  phi <- numeric(Nout)
  omg <- numeric(Nout)
  acc <- 0L
  lglk <- numeric(Nout)

  ## Starting values
  if (missing(zstart)) {
    zstart <- y/l
    zstart <- linkfcn(zstart, linkp, "gaussian")
    ## zstart <- pmax(zstart, -1e8) + rnorm(k, 0, sqrt(ssqsc))
  }
  z[, 1] <- zstart
  if (missing(phistart)) {
    if (phisc == 0) {
      stop ("Argument phistart needed for fixed phi")
    } else {
      if(phipars[2] == -1) {
        tmp <- .1/abs(phipars[3])
      } else {
        tmp <- abs((phipars[2]+1)/phipars[3])
      }
      phistart <- phipars[4] + phipars[1]*gamma(tmp)/
        gamma(phipars[2]/phipars[3])
    }
  } else {
    phistart <- as.double(phistart)
  }
  phi[1] <- phistart
  if (missing(omgstart)) {
    if (omgsc == 0) {
      stop ("Argument omgstart needed for fixed omg")
    } else {
      if(omgpars[2] == -1) {
        tmp <- .1/abs(omgpars[3])
      } else {
        tmp <- abs((omgpars[2]+1)/omgpars[3])
      }
      omgstart <- omgpars[4] + omgpars[1]*gamma(tmp)/
        gamma(omgpars[2]/omgpars[3])
    }
  } else {
    omgstart <- as.double(omgstart)
  }
  omg[1] <- omgstart

  ## Run code
  if (test > 0) { # Running a test
    if (is.logical(test)) test <- 100
    test <- as.integer(test)
    tm <- system.time({
      RUN <- .Fortran("trgasamtry", ll = lglk, z = z, phi = phi, omg = omg,
                      acc = acc,
                      as.double(ybar), as.double(l), as.double(F), 
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsqsc),
                      as.double(phipars),
                      as.double(phisc), as.double(omgpars), as.double(omgsc),
                      as.double(kappa), as.integer(icf), 
                      as.double(nu), as.double(dm), as.integer(Nout),
                      as.integer(test), as.integer(k), as.integer(p))
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    beta <- NULL
    ssq <- NULL
    phi <- RUN$phi
    attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
    attr(omg, 'fixed') <- omgsc == 0
    attr(nu, 'fixed') <- TRUE
    acc_ratio <- RUN$acc/Nout
    out <- list(z = zz0, beta = beta, ssq = ssq, phi = phi, omg = omg, nu = nu, 
                logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
                Nout = Nout, Nbi = Nbi, Nthin = Nthin,
                response = y, weights = l, modelmatrix = F, family = family,
                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
                corrfcn = corrfcn, kappa = kappa, 
                tsqdf = tsqdf, tsqsc = tsqsc,
                locations = loc[ii, , drop = FALSE],
                whichobs = ii)
  } else {
    tm <- system.time({
      RUN <- .Fortran("trgasample", ll = lglk, z = z, z0 = z0,
                      mu = z, mu0 = z0, beta = beta, ssq = ssq,
                      tsq = tsq, phi = phi, omg = omg, acc = acc,
                      as.double(ybar),
                      as.double(l), as.double(F), as.double(F0),
                      as.double(betm0), as.double(betQ0), as.double(ssqdf),
                      as.double(ssqsc), as.double(tsqdf), as.double(tsqsc),
                      as.double(phipars), as.double(phisc),
                      as.double(omgpars), as.double(omgsc), as.double(kappa),
                      as.integer(icf), 
                      as.double(linkp), as.double(dm), as.double(dmdm0),
                      as.integer(Nout), as.integer(Nbi), as.integer(Nthin),
                      as.integer(k), as.integer(k0), as.integer(p))
    })
    ## Store samples
    ll <- RUN$ll
    zz0 <- mm0 <- matrix(NA, NROW(yy), Nout)
    zz0[ii, ] <- RUN$z
    zz0[!ii, ] <- RUN$z0
    mm0[ii, ] <- RUN$mu
    mm0[!ii, ] <- RUN$mu0
    beta <- RUN$beta
    ssq <- RUN$ssq
    tsq <- RUN$tsq
    phi <- RUN$phi
    attr(phi, 'fixed') <- phisc == 0
    omg <- RUN$omg
    attr(omg, 'fixed') <- omgsc == 0
    attr(nu, 'fixed') <- TRUE
    acc_ratio <- RUN$acc/(Nout*Nthin + max(Nthin, Nbi))
    out <- list(z = zz0, mu = mm0, beta = beta, ssq = ssq, tsq = tsq,
                phi = phi, omg = omg, nu = nu,
                logLik = ll, acc_ratio = acc_ratio, sys_time = tm,
                Nout = Nout, Nbi = Nbi, Nthin = Nthin,
                response = ybar, weights = l, modelmatrix = F, family = family,
                betm0 = betm0, betQ0 = betQ0, ssqdf = ssqdf, ssqsc = ssqsc,
                corrfcn = corrfcn, kappa = kappa, 
                tsqdf = tsqdf, tsqsc = tsqsc,
                locations = loc[ii, , drop = FALSE],
                whichobs = ii)
  }
  class(out) <- "geomcmc"
  out
}


##' Convert to an \code{\link[coda]{mcmc}} object.
##'
##' This function takes as input the one or more output(s) from
##' function \code{\link{mcsglmm}} or \code{\link{mcstrga}} and
##' returns an \code{\link[coda]{mcmc}} object or an
##' \code{\link[coda]{mcmc.list}} object for coda. The function
##' requires the \code{coda} package to be installed.

##' The spatial random field components are assigned the names
##' \code{z_*} where \code{*} is a number beginning at 1. Similarly,
##' the regressor coefficients are assigned the names \code{beta_*} if
##' not unique, or simply \code{beta} if there is only one regressor.
##' The names \code{ssq}, \code{tsq}, \code{phi}, \code{omg}
##' correspond to the partial sill, measurement error variance,
##' spatial range, and relative nugget parameters respectively.
##' @title Convert to an \code{\link[coda]{mcmc}} object
##' @param ... Output(s) from the functions mentioned in the Details.
##' @return An mcmc object.
##' @examples \dontrun{
##' ### Load the data
##' data(rhizoctonia)
##' rhiz <- na.omit(rhizoctonia)
##' rhiz$IR <- rhiz$Infected/rhiz$Total # Incidence rate of the
##'                               # rhizoctonia disease
##' 
##' ### Define the model
##' corrf <- "spherical"
##' ssqdf <- 1
##' ssqsc <- 1
##' tsqdf <- 1
##' tsqsc <- 1
##' betm0 <- 0
##' betQ0 <- diag(.01, 2, 2)
##' phiprior <- c(200, 1, 1000, 100) # U(100, 300)
##' phisc <- 1
##' omgprior <- c(3, 1, 1000, 0) # U(0, 3)
##' omgsc <- 1.3
##' linkp <- 1
##' 
##' ## MCMC parameters
##' Nout <- 100
##' Nbi <- 0
##' Nthin <- 1
##'
##' ### Run MCMC
##' sample <- mcstrga(Yield ~ IR, data = rhiz, 
##'                   atsample = ~ Xcoord + Ycoord, corrf = corrf, 
##'                   Nout = Nout, Nthin = Nthin,
##'                   Nbi = Nbi, betm0 = betm0, betQ0 = betQ0,
##'                   ssqdf = ssqdf, ssqsc = ssqsc,
##'                   tsqdf = tsqdf, tsqsc = tsqsc,
##'                   phipars = phiprior, omgpars = omgprior,
##'                   linkp = linkp,
##'                   phisc = phisc, omgsc = omgsc, test=FALSE)
##' 
##' mcsample <- mcmcmake(sample)
##' plot(mcsample[, c("phi", "omg", "beta_1", "beta_2", "ssq", "tsq")],
##'      density = FALSE)
##' summary(mcsample[, c("phi", "omg", "beta_1", "beta_2", "ssq", "tsq")])
##' }
##' @seealso Functions such as \code{\link[coda]{plot.mcmc}} and
##' \code{\link[coda]{summary.mcmc}} in the \code{coda} package. The
##' function \code{\link[base]{do.call}} can be used to pass arguments
##' stored in a list.
##' @importFrom coda mcmc mcmc.list
##' @export 
mcmcmake <- function (...) {
  ### This function takes as input the output from function mcmcrun
  ### and returns an mcmc object or an mcmc.list object for coda. The
  ### function requires the coda package to be installed.
  varnm <- function(x,prefix) {
    # Function to determine variable name.
    if (is.null(x)) return ()
    dx <- dim(x)
    if (length(dx) >= 2) {
      dnmx2 <- dimnames(x)[[2]]
      if (is.null(dnmx2)) {
        out <- if (dx[2] > 1) paste(prefix,'_',seq(dx[2]),sep='') else prefix
      } else {
        out <- dnmx2
      }
    } else out <- prefix
    out
  }
  vnm <- c('z','beta','ssq','tsq','phi','omg')
  input <- list(...)
  nruns <- length(input)
  mcl <- list(); length(mcl) <- nruns
  for (j in seq_len(nruns)) {
    chain <- input[[j]][vnm]
    names(chain) <- vnm
    if (is.matrix(chain$z)) chain$z <- t(chain$z)
    if (is.matrix(chain$z0)) chain$z0 <- t(chain$z0)
    if (is.matrix(chain$beta)) chain$beta <- t(chain$beta)
    dnm <- unlist(sapply(seq_along(vnm),function(i) varnm(chain[[i]],vnm[i])))
    thin <- max(1, input[[j]]$Nthin)
    start <- max(thin, input[[j]]$Nbi)
    mcl[[j]] <- coda::mcmc(matrix(unlist(chain),ncol = length(dnm),
                                  dimnames = list(NULL,dnm)),
                           start = start, thin = thin)
  }
  if (nruns == 1) {
    return (mcl[[1]])
  } else if (nruns > 1) {
    return (do.call(coda::mcmc.list, mcl))
  } else return (NULL)
}
