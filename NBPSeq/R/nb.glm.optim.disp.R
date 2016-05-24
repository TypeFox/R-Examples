##' Estimate the parameters in a dispersion model.
##'
##' The function will call the R funciton \code{optim} to mimimize the
##' negative log likelihood of the dipserison model.
##'
##' @title (private) Estimate the parameters in a dispersion model
##'
##' @param disp a list, output from \code{\link{disp.nbp}}, \code{\link{disp.nbq}}, \code{\link{disp.nbs}}, and so on.
##' @param counts a matrix, the nb counts
##' @param eff.lib.sizes effective library sizes
##' @param x a desing matrix
##' @param method the optimization method to be used by \code{optim}
##' @param ... additional parareters, will be passed to optim().
##' @return a list with components:
optim.disp.pl = function(disp, counts, eff.lib.sizes, x, method="L-BFGS-B", ...) {

  y = counts[disp$subset, ,drop=FALSE];

  ## mustart = NULL;
  mustart = irls.nb(y, eff.lib.sizes, x, phi=0.1, beta0 = rep(NA, dim(x)[2]))$mu;

  ## Negative log likelihood of the dispersion model
  nll = function(par) {
    ## Compute phi as a function of par
    phi = disp$fun(par)[disp$subset,,drop=FALSE];

    ## Fit NB regression models to the rows of the count matrix
    mu.hat = irls.nb(y, eff.lib.sizes, x, phi=phi, beta0 = rep(NA, dim(x)[2]), mustart=mustart)$mu;

    ## This line has no effect!
    ## mustart = mu.hat;

    ## Compute the negative log likelihood of the fitted model
    ## kappa = 1/phi;

    - ll.nb(1/phi, mu.hat, y); 
  }


  if (is.null(disp$lower) | is.null(disp$upper)) {
    res = optim(disp$par.init,  nll, method=method, ...);
  } else {
    res = optim(disp$par.init,  nll, lower=disp$lower, upper=disp$upper,
      method=method, ...);
  }

  res
}

##' Estimate the parameters in a dispersion model.
##'
##' The function will call the R funciton \code{optim} to mimimize the
##' negative log adjusted profile likelihood of the dipserison model.
##'
##' @title (private) Estimate the parameters in a dispersion model
##'
##' @param disp a list, output from \code{\link{disp.nbp}}, \code{\link{disp.nbq}}.
##' @param counts a matrix, the nb counts
##' @param eff.lib.sizes effective library sizes
##' @param x a desing matrix
##' @param method the optimization method to be used by \code{optim}
##' @param print.level print level
##' @param ... additional parareters, will be passed to optim().
##' @return a list with components:
optim.disp.apl = function(disp, counts, eff.lib.sizes, x, method="L-BFGS-B", print.level=1, ...) {

  y = counts[disp$subset, ,drop=FALSE];

  mustart = irls.nb(y, eff.lib.sizes, x, phi=0.1, beta0 = rep(NA, dim(x)[2]))$mu;

  ## Adjustment term in apl
  adj = function(y, mu.hat, phi) {
    kappa = 1/phi;
    v.hat = drop(mu.hat + phi * mu.hat^2);
    j.hat = t(x) %*% diag(mu.hat^2 * (y/mu.hat^2 - (y + kappa)/(mu.hat + kappa)^2) - (y - mu.hat) * mu.hat/v.hat) %*% x;
    d = det(j.hat);

    ## Sometimes d can be negative 
    if (d > 0) {
      return (-0.5 * log(d))
    } else {
      print(y);
      print(mu.hat);
      print(phi);
      return (0)
    }
  }
  
  ## Negative log likelihood of the dispersion model
  nll = function(par) {
    ## Compute phi as a function of par
    phi = disp$fun(par)[disp$subset,,drop=FALSE];

    ## Fit NB regression models to the rows of the count matrix
    mu.hat = irls.nb(y, eff.lib.sizes, x, phi=phi, beta0 = rep(NA, dim(x)[2]),
      mustart=mustart)$mu;

    ## Compute the negative log likelihood of the fitted model
    l.hat = ll.nb(1/phi, mu.hat, y);

    m = dim(mu.hat)[1];
    for (i in (1:m)) {
      l.hat = l.hat + adj(y[i,], mu.hat[i,], phi[i,]);
    }

    -l.hat
  }


  if (is.null(disp$lower) | is.null(disp$upper)) {
    res = optim(disp$par.init,  nll, method=method, ...);
  } else {
    res = optim(disp$par.init,  nll, lower=disp$lower, upper=disp$upper,
      method=method, ...);
  }

  res
}

