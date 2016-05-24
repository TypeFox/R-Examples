##' Specify a NBP dispersion model. The parameters of the specified
##' model are to be estimated from the data using the function
##' \code{optim.disp.apl} or \code{optim.disp.pl}.
##'
##' Under this NBP model, the log dispersion is modeled as a linear
##' function of the preliminary estimates of the log mean realtive
##' frequencies (\code{pi.pre}):
##'
##' log(phi) =  par[1] + par[2] * log(pi.pre/pi.offset),
##'
##' where \code{pi.offset} is 1e-4.
##'
##' Under this parameterization, par[1] is the dispersion value when
##' the estimated relative frequency is 1e-4 (or 100 RPM).
##'
##' @title (private) Specify a NBP dispersion model
##' @internal
##' @param counts an \eqn{m \times n} matrix of NB counts
##' @param eff.lib.sizes a \eqn{n}-vector of estimated effective library sizes
##' @param x a nxp matrix, design matrix (specifying the treatment structure)
##' @param phi.pre a number, a preliminary constant dispersion value,
##' will be used to get preliminary estimates of mean counts (mu.pre).
##' @param mu.lower  a number, rows with any component of mu.pre < mu.lower will not be used for estimating the dispersion model 
##' @param mu.upper  a number, rows with any component of mu.pre > mu.upper will not be used for estimating the dispersion model
##' @return a list
##'
##' \item{fun}{a function that takes a vector, par, as
##' input and outputs a matrix of dispersion values (same dimension as
##' counts)}
##' \item{par.init}{a numeric vector of length 2, initial values of par}
##' \item{lower}{a numeric vector of length 2, lower bounds of the parameter values}
##' \item{upper}{a numeric vector of length 2, upper bounds of the parameter values}
##' \item{subset}{a logical vector of length \eqn{m}, specifying the subset of rows to be used when estimating the dispersion model parameters.}
##' \item{pi.pre}{a m by n matrix, preliminary estimates of the relative frequencies.}
##' \item{pi.offset}{a scalar, fixed to be 1e-4, an offset used in the NBP model (see Details)}
disp.nbp = function(counts, eff.lib.sizes, x, phi.pre=0.1, mu.lower=1, mu.upper=Inf) {
  pi.offset = 1e-4;

  ## Fit NB regression models to the rows of the count matrix using a
  ## preliminary constant dispersion parameter, phi0.
  obj = irls.nb(counts, eff.lib.sizes, x, phi=phi.pre, beta0 = rep(NA, dim(x)[2]));

  ## Obtain preliminary estimates of the mean relative frequencies.
  mu.pre = obj$mu;
  pi.pre = t(t(mu.pre)/eff.lib.sizes);

  ## Replace extremely small pi values 
  ## 11-22-2013
  eps = 1/sum(eff.lib.sizes);
  pi.pre[pi.pre<eps] = eps;

  ## Specify the subset of rows to be used when estimating the dispersion model
  subset = rowSums(is.na(mu.pre) | mu.pre<mu.lower | mu.pre>mu.upper)==0;

  ## Contruct predicting vector of the NBP model: a linear function of
  ## the preliminary estimtes of the log mean realtive frequencies
  z = log(pi.pre / pi.offset);

  ## Specify the dispersion function as a function of par
  env = new.env(parent=baseenv());
  assign("z", z, envir=env)
  fun=function(par){
    exp(par[1] + par[2] * z)
  }
  environment(fun) = env;

  ## Specify the lower and upper bounds for the components of par
  lower = c(log(1e-20), -1.1);
  upper = c(0, 0.1);

  ## Sepcify the initial values of par
  par.init = c(log(0.1), 0);

  list(name="NBP", fun=fun, par.init=par.init, lower=lower, upper=upper, subset=subset, pi.pre=pi.pre, pi.offset = pi.offset);
}

##' Specify a NBQ dispersion model. The unknown parameters in the
##' specified model are to be estimated using the function
##' \code{optim.disp.pl} or \code{optim.disp.apl}.
##'
##' Under this NBQ model, the dispersion is modeled as a quadratic
##' function of the preliminary estimates of the log mean realtive
##' frequencies (pi.pre):
##'
##' log(phi) =  par[1] + par[2] * z + par[3] * z^2,
##'
##' where z = log(pi.pre/pi.offset). By default, pi.offset is the median of pi.pre[subset,].
##'
##' @title (private) Specify a NBQ dispersion model
##' @internal
##' @param counts a mxn matrix of NB counts
##' @param eff.lib.sizes a n-vector of estimated effective library sizes
##' @param x a nxp matrix, design matrix (specifying the treatment structure)
##' @param phi.pre a number, a preliminary constant dispersion value
##' @param mu.lower a number, rows with mu.pre < mu.lower will not be used for estimating the dispersion model 
##' @param mu.upper a number, rows with mu.pre > mu.upper will not be used for estimating the dispersion model
##' @param pi.offset  a scalar, an offset used in the NBQ model (see Details).
##' @return a list
##'
##' \item{fun}{a function that takes a vector, par, as
##' input and outputs a matrix of dispersion values (same dimension as
##' counts)}
##'
##' \item{par.init}{a vector of length 3, initial values of par}
##' \item{lower}{a vector of length 3, lower bounds of the parameter values}
##' \item{upper}{a vector of length 3, upper bounds of the parameter values}
##' \item{subset}{a logical vector of length \eqn{m}, specifying the subset of rows to be used when estimating the dispersion model parameters.}
##' \item{pi.pre}{a m by n matrix, preliminary estimates of the relative frequencies.}
##' \item{pi.pre}{a m by n matrix, preliminary estimates of the relative frequencies.}
##' \item{pi.offset}{a scalar used as an offset in the NBQ model (see Details)}
disp.nbq = function(counts, eff.lib.sizes, x, phi.pre=0.1, mu.lower=1, mu.upper=Inf,
  pi.offset=median(pi.pre[subset,])) {

  ## Fit NB regression models to the rows of the count matrix using a
  ## preliminary constant dispersion parameter, phi0.
  obj = irls.nb(counts, eff.lib.sizes, x, phi=phi.pre, beta0 = rep(NA, dim(x)[2]));

  ## Obtain preliminary estimates of the mean relative frequencies.
  mu.pre = obj$mu;
  pi.pre = t(t(mu.pre)/eff.lib.sizes);

  ## Replace extremely small pi values 
  ## 11-22-2013
  eps = 1/sum(eff.lib.sizes);
  pi.pre[pi.pre<eps] = eps;

  ## Specify the subset of rows to be used when estimating the dispersion model
  subset = rowSums(is.na(mu.pre) | mu.pre<mu.lower | mu.pre>mu.upper)==0;

  ## Contruct predicting vectors of the NBQ model: they are linear and
  ## quadratic functions of the preliminary estimtes of the log mean
  ## realtive frequencies
  z = log(pi.pre / pi.offset);

  ## Specify the dispersion function as a function of par
  env = new.env(parent=baseenv());
  assign("z", z, envir=env)
  fun = function(par) {
    exp(par[1] + par[2] * z + par[3] * z^2);
  }
  environment(fun) = env;

  ## Specify the lower and upper bounds for the components of par
  lower = c(log(1e-20), -1.0, -0.2);
  upper = c(0, 1.0, 0.2);

  ## Sepcify the initial values of par
  par.init = c(log(0.1), 0, 0);

  list(name="NBQ", fun=fun, par.init=par.init, lower=lower, upper=upper, subset=subset, pi.pre=pi.pre, pi.offset=pi.offset);
  ## list(fun=fun, lower=lower, upper=upper, par.init=par.init, subset=subset);
}


##' Specify a NBS dispersion model. The specified model are to be
##' estimated using the function \code{optim.disp.pl} or
##' \code{optim.disp.apl}.  Under this NBS model, the dispersion is
##' modeled as a smooth function (a natural cubic spline function) of
##' the preliminary estimates of the log mean realtive frequencies
##' (pi.pre).
##'
##' \code{disp.nbs} calls the function \code{ns} to generate a set of
##' spline bases, using log(pi.pre) (converted to a vector) as the
##' predictor variable. Linear combinations of these spline bases are
##' smooth functions of log(pi.pre). The return value includes a
##' function, \code{fun}, to be optimized by \code{optim.disp.pl} or
##' \code{optim.disp.apl}. The parameter of that function is a vector
##' of linear combination coefficients of the spline bases.
##'
##' \code{df+2} nodes are used when constructing the splint bases. The
##' Boundary.nodes are placed at the min and max values of
##' log(pi.pre). Two nodes are placed at the 0.05 and 0.95th quantiles
##' of log(pi.pre) and an additional \code{df-2} inner nodes are
##' equally spaced between the two nodes.
##'
##' It is a challenging issue to determine the optimal number and
##' placement of the nodes.
##'
##' @title (private) Specify a NBS dispersion model
##' @internal
##' @param counts a mxn matrix of NB counts
##' @param eff.lib.sizes a n-vector of estimated effective library sizes
##' @param x a nxp matrix, design matrix (specifying the treatment structure)
##' @param df the number of interior nodes
##' @param phi.pre a number, a preliminary constant dispersion value
##' @param mu.lower a number, rows with mu.pre < mu.lower will not be used for estimating the dispersion model 
##' @param mu.upper a number, rows with mu.pre > mu.upper will not be used for estimating the dispersion model
##' @return a list
##'
##' \item{fun}{a function that takes a vector, par, as input and
##' outputs a matrix (same dimension as \code{counts}) of dispersion
##' values. \code{par} will be used as linear-combination efficients
##' for the spline bases, the estimated dispersion values are a spline
##' function of log(pi.pre).  }
##'
##' \item{par.init}{initial values of par}
##'
##' \item{subset}{a logical vector of length \eqn{m}, specifying the
##' subset of genes to be used when estimating the model
##' parameters. Note that the estimated model will be applied to all
##' rows whenever possible, but only rows specified in \code{subset}
##' will be used to estimate the parameters of the dispersion model.}
##'
##' \item{pi.pre}{a m by n matrix, preliminary estimates of the
##' relative frequencies.}
##'
##' \item{s}{Basis matrix of the natural cubic spline evalued at the z = log(pi.pre)}
disp.nbs = function(counts, eff.lib.sizes, x,
  df = 6,
  phi.pre=0.1, mu.lower=1, mu.upper=Inf) {

  ## Fit NB regression models to the rows of the count matrix using a
  ## preliminary constant dispersion parameter, phi0.
  obj = irls.nb(counts, eff.lib.sizes, x, phi=phi.pre, beta0 = rep(NA, dim(x)[2]));

  ## Obtain preliminary estimates of the mean relative frequencies.
  mu.pre = obj$mu;
  pi.pre = t(t(mu.pre)/eff.lib.sizes);

  ## Replace extremely small pi values 
  ## 11-26-2013
  eps = 1/sum(eff.lib.sizes);
  pi.pre[pi.pre<eps] = eps;

  ## Specify the subset of rows to be used when estimating the dispersion model
  subset = rowSums(is.na(mu.pre) | mu.pre<mu.lower | mu.pre>mu.upper)==0;

  ## The predictor vector is the initial estimates of the log mean relative frequencies
  ## 2013-12-13
  z = as.vector(log(pi.pre));

  ## Specify the boundary knots
  Boundary.knots = range(z[subset]);
  ## Boundary.knots = c(1e-4, quantile(z, 0.99));

  ## Specify the knots
  ## It is not clear yet how to select the knots
  zs = sort(z[subset]);
  m = length(zs);
  l = quantile(z[subset], 0.05);
  r = quantile(z[subset], 0.95);

##  if (m > 300) {
##   l = max(l, zs[101]+0.01);
##   r = min(r, zs[m-100]);
##  }

  knots = seq(l, r, length=df-2);

  ## Specify the spline basis by the providing the knots
  ## s = ns(z, df = df);
  s = ns(z, knots=knots, Boundary.knots = Boundary.knots, intercept=TRUE);
  d = dim(pi.pre);

  ## Specify the dispersion function as a function of par
  ## Create a minimal enviroment for the function
  env = new.env(parent=baseenv());
  assign("s", s, envir=env)
  assign("d", d, envir=env)

  fun = function(par) {
    phi = exp(s %*% par);
    dim(phi) = d;
    phi
  }
  environment(fun) = env;

  ## Specify the lower and upper bounds for the components of par
  ## Not needed

  ## Sepcify the initial values of par
  ## This set seems to work
  par.init = rep(-1, dim(s)[2]);

  list(name="NBS", fun=fun, par.init=par.init, subset=subset, pi.pre=pi.pre, s=s);
}


##' Specify a piecewise constant (step) dispersion model. The
##' specified model are to be estimated using the function
##' \code{optim.disp.pl} or \code{optim.disp.apl}.
##'
##' Under this model, the dispersion is modeled as a step (piecewise constant) function.
##'
##' @title (private) Specify a piecewise constant dispersion model
##' @internal 
##' @param counts a mxn matrix of NB counts
##' @param eff.lib.sizes a n-vector of estimated effective library sizes
##' @param x a nxp matrix, design matrix (specifying the treatment structure)
##' @param df the number of steps
##' @param knots  a numerical vector of length df-1, giving the knots or jump locations.
##' @param phi.pre a number, a preliminary constant dispersion value
##' @param mu.lower a number, rows with mu.pre < mu.lower will not be used for estimating the dispersion model 
##' @param mu.upper a number, rows with mu.pre > mu.upper will not be used for estimating the dispersion model
##' @return a list \item{fun}{a function that takes par as input and
##' outputs a matrix (same dimension as \code{counts}) of dispersion
##' values}
##' \item{par.init}{a vector of length \code{df}, initial values of par}
##' \item{subset}{a logical vector of length m, specifying a subset of genes to be used when estimating the model parameters. Note that the estimated model will be applied to all rows whenever possible, but only rows specified in \code{subset} will be used to estimate the dispersion model parameters.}
##' \item{pi.pre}{a m by n matrix, preliminary estimates of the relative frequencies.}
##' \item{knots}{a vector, the break points of the step function}
disp.step = function(counts, eff.lib.sizes, x,
  df = 1,
  knots = NULL,
  phi.pre=0.1, mu.lower=1, mu.upper=Inf) {

  ## Fit NB regression models to the rows of the count matrix using a
  ## preliminary constant dispersion parameter, phi0.
  obj = irls.nb(counts, eff.lib.sizes, x, phi=phi.pre, beta0 = rep(NA, dim(x)[2]));

  ## Obtain preliminary estimates of the mean relative frequencies.
  mu.pre = obj$mu;
  pi.pre = t(t(mu.pre)/eff.lib.sizes);

  ## Specify the subset of rows to be used when estimating the dispersion model
  subset = rowSums(is.na(mu.pre) | mu.pre<mu.lower | mu.pre>mu.upper)==0;

  ## Replace extremely small pi values 
  ## 11-26-2013
  eps = 1/sum(eff.lib.sizes);
  pi.pre[pi.pre<eps] = eps;

  ## The predictor vector is the initial estimates of the log mean relative frequencies
  ## z = log(pi.pre*1000);
  ## 2013-12-13
  z = log(pi.pre);

  if (is.null(knots)) {
    p = seq(0, 1, length=df+1)[-c(1, df+1)];
    knots = quantile(z[subset], p);
  } else {
    knots = sort(knots);
  }

  ## Specify the starting and ending positions of each step
  l = c(min(z, knots)-0.01, knots);
  r = c(knots, max(z, knots)+0.01);

  ## Identify indices of z values belonging to each step
  ids = list(df);
  m = length(z);
  for (i in 1:df) {
    ids[[i]] = (1:m)[z>=l[i] & z< r[i]];
  }

  d = dim(pi.pre);

  ## Specify the dispersion function as a function of par
  ## Create a minimal enviroment for the function
  env = new.env(parent=baseenv());
  assign("df", df, envir=env)
  assign("ids", ids, envir=env)
  assign("d", d, envir=env)

  fun = function(par) {
    lphi = array(NA, d);

    for (i in 1:df) {
      lphi[ids[[i]]] = par[i];
    }

    exp(lphi)

  }
  environment(fun) = env;

  ## Specify the lower and upper bounds for the components of par
  ## Not needed

  ## Sepcify the initial values of par
  ## This set seems to work
  par.init = rep(-1, df);

  list(name="step", fun=fun, par.init=par.init, subset=subset, pi.pre=pi.pre, knots = knots);
}
