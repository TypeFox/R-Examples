##' Specify a NB2 dispersion model. The parameter of the specified
##' model are to be estimated from the data using the function
##' \code{optim.pcl}.
##'
##' Under this NB2 model, the log dispersion is a constant.
##'
##' @title (private) A NBP dispersion model
##'
##' @param par a number, log dispersion
##' @return a vector of length m, log dipserison.
log.phi.nb2=function(par, pi){
  rep(par, length(pi))
}

##' Specify a NBP dispersion model. The parameters of the specified
##' model are to be estimated from the data using the function
##' \code{optim.pcl}.
##'
##' Under this NBP model, the log dispersion is modeled as a linear
##' function of the log mean realtive
##' frequencies (\code{pi.pre}):
##'
##' log(phi) =  par[1] + par[2] * log(pi.pre/pi.offset),
##'
##' where the default value of \code{pi.offset} is 1e-4.
##'
##' @title (private) A NBP dispersion model
##'
##' @param par a vector of length 2, the intercept and the slope of the log linear model (see Details).
##' @param pi a vector of length m, estimated mean relative frequencies
##' @param pi.offset a number
##' @return a vector of length m, log dipserison.
log.phi.nbp=function(par, pi, pi.offset=1e-4){
  z = log(pi / pi.offset);
  par[1] + par[2] * z;
}

##' Specify a NBQ dispersion model. The parameters of the specified
##' model are to be estimated from the data using the function
##' \code{optim.pcl}.
##'
##' Under this NBQ model, the log dispersion is modeled as a quadratic
##' function of the log mean realtive
##' frequencies (\code{pi}):
##'
##' log(phi) =  par[1] + par[2] * log(pi/pi.offset) + par[3] * (log(pi/pi.offset))^2;
##'
##' where the (default) value of \code{pi.offset} is 1e-4.
##'
##' @title (private) A NBQ dispersion model
##' @noRd
##'
##' @param par a vector of length 3, see Details.
##' @param pi a vector of length m, estimated mean relative frequencies
##' @param pi.offset a number
##' @return a vector of length m, log dipserison.
log.phi.nbq=function(par, pi, pi.offset=1e-4){
  z = log(pi / pi.offset);
  par[1] + par[2] * z + par[3] * z^2;
}

##' Negative profile conditional likelihood of the dispersion model
##'
##' @title Negative profile conditional likelihood of the dispersion model
##' @param par  parameter of the disperesion model
##' @param log.phi.fun the disperison model
##' @param y counts
##' @param ls library sizes
##' @param n.grps number of groups
##' @param grps a boolean matrix of group membership
##' @param grp.sizes group sizes
##' @param mu.lower lower bound for mu
##' @param mu.upper upper bound for mu
##' @param print.level  print level
##' @return negative log likelihood of the dispersion function 
nll.log.phi.fun = function(par, log.phi.fun, y, ls, n.grps, grps, grp.sizes,
  mu.lower, mu.upper, print.level) {

  if (print.level>3)  print(par);
  
  l = double(n.grps);

  ## Compute the log conditional likelihood group by group
  for (i in 1:n.grps) {
    r = grp.sizes[i]; 

    ## We can only compute the conditional likelihood for groups with replicates
    if ( r > 1) {
      ## Compute the n vector of mean counts.
      mu = rowMeans(y[, grps[i,]]);
      pi = mu / ls;

      ## The size/shape parameter is the reciprical of the dispersion
      kappa = exp(-log.phi.fun(par, pi));

      ## thresholding
      ss = (mu > mu.lower) & (mu < mu.upper); # an n vector

      kappa = kappa[ss];
      ## print(range(1/kappa));

      y0 = y[,grps[i,]][ss,];
      
      l[i] = sum(lgamma(kappa + y0)); # kappa is n x 1, y0 is n x r 

      t = rowSums(y0);

      l[i] = l[i] + sum(lgamma(r * kappa) - r * lgamma(kappa) - lgamma(r * kappa + t));
    }
  }

  if (print.level > 3) cat(sprintf("- log likelihood of the dispersion model: %f.\n", -sum(l)));

  - sum(l);
}


##' Estimate parameters in a dispersion model by optimizing the
##' profile conditional likelihood.
##'
##'
##' This function only works with
##' object prepared \code{\link{prepare.nbp}}.
##' This function will use the pseudo counts in the \code{obj} which
##' are thinned (subsampled) raw counts that have equal library sizes.
##' 
##' In the input, the function \code{log.phi.fun} specifies a
##' parametric dispersion model that model the dispersion as
##' parametric function of the mean relative frequencies.  This
##' function takes \code{par} as the first parameter and a vector of
##' mean relative frequencies a second parameter. It returns a vector
##' of dispersion values according to the specified dispersion model.
##'
##' This function determines the parameters of the parametric dispersion model
##' (specified by \code{log.phi.fun}) by maximizing the conditional
##' likelihood.
##'
##' @title (private) Estimate parameters in a parametric dispersion model 
##' @noRd
##'
##' @param obj output from \code{\link{prepare.nbp}}.
##' @param log.phi.fun a disperison function, one of
##' \code{log.phi.nb2}, \code{log.phi.nbp} or \code{log.phi.nbq}.
##'
##' @param par.init  a vector, initial values of par 
##' @param par.lower a vector, lower bounds of par
##' @param par.upper a vector, upper bounds of par
##' @param mu.lower see mu.upper
##' @param mu.upper only genes with mean counts between  mu.lower and
##' mu.upper will be used to compute the conditional likilihood of the
##' dispersion parameters
##'
##' @return a list, output from \code{optim} with an additional component:
##' \item{log.phi.fun}{same as input}
optim.pcl = function(obj, log.phi.fun,
  par.init, par.lower, par.upper,
  mu.lower=1, mu.upper=Inf,
  method="L-BFGS-B",
  print.level=0, ...) {

  if (print.level > 0) cat("Estimate NB dispersion.\n");

  ## Use thinned counts
  y = obj$pseudo.counts;
  if (is.matrix(y)==FALSE) stop ("counts not a matrix");

  lib.sizes = obj$pseudo.lib.sizes;
  if (diff(range(lib.sizes)) / mean(lib.sizes) > 0.01)
    stop("The pseudo.lib.sizes are not the same!");

  ## Number of genes
  m = dim(y)[1];	

  ## Number of samples
  n = dim(y)[2];		

  ## The common library sizes
  ls = lib.sizes[1]; 

  ## Identify the treatment groups
  gids = unique(obj$grp.ids);  
  n.grps = length(gids);  
  grps = matrix(FALSE, n.grps, n);
  grp.sizes = integer(n.grps);
  for (i in  1:n.grps) {
    grps[i,] = (obj$grp.ids == gids[i]);
    grp.sizes[i] = sum(grps[i,]);
  }

  nll = function(par) {nll.log.phi.fun(par, log.phi.fun, y, ls, n.grps, grps, grp.sizes,
    mu.lower, mu.upper, print.level = print.level);}

  res = optim(par.init, nll, lower=par.lower, upper=par.upper, method = method, ...);

  c(list(log.phi.fun=log.phi.fun), res)
} 
