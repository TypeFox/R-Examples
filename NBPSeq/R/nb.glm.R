## List of functions:
##
## irls.nb.1 = function(y, s, x, phi, beta0=rep(NA,p),
##   maxit=50, tol.mu=1e-3/length(y), print.level=0);
##
## irls.nb = function(y, s, x, phi, beta0, maxit=50, tol.mu=0.01, print.level=1);
##
##
## irls.nbp.1 = function(y, s, x, phi0, alpha1, beta0=rep(NA, p),
##   maxit=50, tol.mu=1e-3/length(y), print.level=1);
##
## irls.nbp = function(y, s, x, phi0, alpha1, maxit=50, tol.mu=0.01, print.level=0);
##
## pl.phi.1 = function(phi, y, s, x, beta0); 
##


##' Estimate the regression coefficients in an NB GLM model with known
##' dispersion parameters
##'
##' This function estimates the regression coefficients using iterative
##' reweighted least squares (IRLS) algorithm, which is equivalent to
##' Fisher scoring. The implementation is based on \code{glm.fit}.
##'
##' Users can choose to fix some regression coefficients by specifying
##' \code{beta0}. (This is useful when fitting a model under a null
##' hypothesis.)
##'
##' @title Estimate the regression coefficients in an NB GLM model
##' @export
##' @param y an n vector of counts
##' @param s a scalar or an n vector of effective library sizes
##' @param x an n by p design matrix
##' @param phi a scalar or an n-vector of dispersion parameters
##' @param mustart starting values for the vector of means
##' @param beta0 a vector specifying known and unknown components of
##' the regression coefficients: non-NA components are hypothesized
##' values of beta, NA components are free components
##' @param maxit maximum number of iterations 
##' @param tol.mu a number, convergence criteria
##' @param print.level a number, print level
##' @return a list of the following components:
##'  \item{beta}{a p-vector of estimated regression coefficients}
##'  \item{mu}{an n-vector of estimated mean values}
##'  \item{conv}{logical. Was the IRLS algorithm judged to have converged?}
##'  \item{zero}{logical. Was any of the fitted mean close to 0?}
irls.nb.1 = function(y, s, x, phi, 
  beta0=rep(NA,p),
  mustart=NULL,
  maxit=50, tol.mu=1e-3/length(y), print.level=0) {

  nobs = as.integer(dim(x)[1]);
  p = as.integer(dim(x)[2]);

  ## Indices to fixed and free components of beta
  id1 = (1:p)[is.na(beta0)];
  id0 = (1:p)[!is.na(beta0)];
  q = length(id0);
  nvars = p - q;

  ## Offset
  beta = beta0;
  offset = matrix(x[, id0], nobs, q) %*% beta[id0];

  ## Initial estiamte of beta
  ## eta = log((y+0.5)/s);
  ## beta[id1] = qr.solve(x[, id1], eta - offset);
  ## eta = drop(x %*% beta);

  if (is.null(mustart)) {
    mu = y + (y==0)/6;
  } else {
    mu = mustart;
  }
  eta = log(mu/s);

  ##------------- THE Iteratively Reweighting L.S. iteration -----------
  conv = FALSE;
  for (iter in 1L:maxit) {
    varmu = mu + phi * mu^2;
    if (any(is.na(varmu)))
      stop("NAs in V(mu)")
    if (any(varmu == 0))
      stop("0s in V(mu)")

    ## mu.eta.val = mu;
    ## if (any(is.na(mu.eta.val))) stop("NAs in d(mu)/d(eta)")
    ## drop observations for which w will be zero
    ## good = (mu.eta.val != 0)
    z = eta - offset + (y - mu)/mu;
    w = drop(mu/sqrt(varmu));

    ## call Fortran code to perform weighted least square
    epsilon = 1e-7;

    ## call Fortran code via C wrapper
    fit = .Call(Cdqrls, x[, id1, drop=FALSE] * w,  w * z, epsilon);

    ##    fit = .Fortran("dqrls",
    ##                    qr = x[,id1] * w, n = nobs,
    ##                    p = nvars, y = w * z, ny = 1L,
    ##                    tol = epsilon,
    ##                    coefficients = double(nvars),
    ##                    residuals = double(nobs),
    ##                    effects = double(nobs),
    ##                    rank = integer(1L),
    ##                    pivot = 1L:nvars,
    ##                    qraux = double(nvars),
    ##                    work = double(2 * nvars),
    ##                    PACKAGE = "base");

    if (any(!is.finite(fit$coefficients))) {
      warning(gettextf("non-finite coefficients at iteration %d", iter), domain = NA);
      break
    }

    ## stop if not enough parameters
    if (nobs < fit$rank)
      stop(gettextf("X matrix has rank %d, but only %d observations",
                    fit$rank, nobs), domain = NA)

    ## calculate updated values of eta and mu with the new coef:
    muold = mu;
    beta[id1[fit$pivot]] = fit$coefficients;
    eta = drop(x %*% beta);
    mu = s*exp(eta);

    ## check for convergence of mu
    ## ! when mu converges, beta may not
    if (max(abs(mu - muold)) < tol.mu) {
      conv = TRUE;
      break
    }
  }

  ## Check whether the any of the means is close to 0
  zero = any(mu < tol.mu * 2); 

  ##-------------- end IRLS iteration -------------------------------
  list(mu=mu, beta=beta, iter=iter, conv=conv, zero=zero);
}


##' Estimate the regression coefficients in an NBP GLM model for each gene
##'
##' @title (private) Estiamte the regression coefficients in an NB GLM model
##' @param y an m*n matrix of counts
##' @param s an n vector of effective library sizes
##' @param x an n*p design matrix
##' @param phi a scalar or an m*n matrix of NB2 dispersion coefficients
##' @param beta0 a K vector, non NA components are hypothesized values of beta, NA components are free components
##' @param mustart an m*n matrix of starting values of the means
##' @param ...  other parameters 
##' @param print.level a number, print level
##' @return beta a K vector, the MLE of the regression coefficients.
irls.nb = function(y, s, x, phi, beta0, mustart=NULL, ..., print.level=0) {
  m = dim(y)[1];
  n = dim(y)[2];
  p = dim(x)[2];

  phi = matrix(phi, m, n);

  if (print.level > 0)
    print("Estimating NB regression coefficients using IRLS.");
  
  res=list(mu=matrix(NA, m, n), beta=matrix(NA, m, p),
    conv=logical(m), iter=numeric(m));

  if (print.level > 1) {
    ## Set up progress bar
    pb=txtProgressBar(style=3);
  }

  for (i in 1:m) {
    if (print.level > 1) {
      setTxtProgressBar(pb, i/m);
    }

    if (is.null(mustart)) {
      res0 = irls.nb.1(y[i,], s, x, phi[i,], beta0,
        ...,
        print.level=print.level-1);
    } else {
      res0 = irls.nb.1(y[i,], s, x, phi[i,], beta0,
        mustart=mustart[i,],
        ...,
        print.level=print.level-1);
    }

    res$mu[i,] =  res0$mu;
    res$beta[i,] = res0$beta;
    res$conv[i] = res0$conv;
    res$iter[i] = res0$iter;
  }

  if (print.level>1) close(pb);

  res;
}



##' Find the MLE of the regression coefficients in an NB regression
##' model using the iteratively reweighted least squres (ILRS)
##' algorithm and compute the Fisher and/or observed information.
##'
##' Under the NB regression model, the components of y follow a NB
##' distribution with means mu = s exp(x' beta) and dispersion
##' parameters phi.
##'
##' The function will call \code{\link{irls.nb.1}} to find MLE of the
##' regression coefficients using the iteratively reweighted least
##' squres (ILRS) algorithm.
##'
##' @note
##'
##' The information matries, i and j, will be computed for all all components
##' of beta---including known components.
##'
##' @title (private) Fit a single negative binomial (NB) log-linear regression model
##' with known dispersion paramreters
##' @noRd
##'
##' @param nb.data output from prepare.nb.data
##' @param dispersion output from estimate.dispersion
##' @param x an n by p design matrix.
##' @param beta0 a p-vector specifying the known and unknown
##' components of beta, the regression coefficients. NA values
##' indicate unknown components and non-NA values specify the values
##' of the known components. The default is that all components of
##' beta are unknown.
##' @param ... furhter arguements to be passed to \code{\link{irls.nb.1}}.
##' @return a list
fit.nb.regression = function(nb.data, dispersion, x, beta0=rep(NA, dim(x)[2]), ...) {

  ## Find MLE of beta
  res = irls.nb(nb.data$counts, nb.data$eff.lib.sizes, dispersion$estimates, x, beta0, ...);

  res
}

