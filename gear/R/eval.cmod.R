#' Evaluate covariance or semivariance model.
#' 
#' \code{eval.cmod} evaluates the covariance or semivariance of a model based on the provided arguments.  See Details for what the function returns, as it changes depending on the class of \code{mod}.
#' 
#' If \code{mod} is of class \code{cmodStd} (from the \code{cmod.std} function), then the function returns an \eqn{n \times m} matrix with the evaluated standard covariance function. 
#' 
#' @param mod A covariance or semivariance model.
#' @param d An \eqn{n \times m} matrix of distances.
#' @param coords An numeric object with two columns.
#' 
#' @return Returns the evaluated model with necessary components needed for \code{mlefit}, \code{predict}, or \code{loglik} functions.
#' 
#' @author Joshua French
#' @examples 
#' n = 10
#' coords = matrix(runif(2*n), nrow = n, ncol = 2)
#' d = as.matrix(dist(coords))
#' cmod = cmod.std(model = "exponential", psill = 1, r = 1)
#' eval.cmod(cmod, d)

#' @rdname eval.cmod
#' @export
eval.cmod = function(mod, d = NULL, coords = NULL)
{
  UseMethod("eval.cmod")
}

#' @rdname eval.cmod
#' @export
eval.cmod.cmodStd = function(mod, d, coords = NULL)
{
  if(!is.numeric(d)) stop("d must be a numeric matrix")
  
  model = mod$model
  psill = mod$psill
  r = mod$r
  evar = mod$evar
  fvar = mod$fvar
  par3 = mod$par3
  
  if(model == "exponential")
  {
    V = psill*exp(-d/r)
    
  }else if(model == "gaussian")
  {
    V = psill*exp(-(d/r)^2)
    
  }else if(model == "matern")
  {
    sd = d/r
    V = (d > 0) * psill*(2^(1-par3)/gamma(par3)*sd^par3*besselK(sd, nu = par3))
    V[is.nan(V)] = psill
  }else if(model == "amatern")
  {
    sd = 2 * sqrt(par3) * d/r
    V = (d > 0) * psill*(2^(1-par3)/gamma(par3)*sd^par3*besselK(sd, nu = par3))
    V[is.nan(V)] = psill	
  }else if(model == "spherical")
  {
    sd = d/r
    V = psill*(1 - 1.5*sd + 0.5*(sd)^3)*(d < r)
  }else if(model == "wendland1")
  {
    sd = d/r
    V = psill*((1 - sd)^4 * (4*sd + 1))*(d < r)
  }else if(model == "wendland2")
  {
    sd = d/r
    V = psill*((1 - sd)^6 * (35*sd^2 + 18*sd + 3))/3*(d < r)
  }else if(model == "wu1")
  {
    sd = d/r
    V = psill*((1 - sd)^3 * (1 + 3*sd + sd^2))*(d < r)
  }else if(model == "wu2")
  {
    sd = d/r
    V = psill*((1 - sd)^4*(4 + 16*sd + 12*sd^2 + 3*sd^3))/4*(d < r)
  }else if(model == "wu3")
  {
    sd = d/r
    V = psill*((1 - sd)^6 * (1 + 6*sd + 41/3*sd^2 + 12*sd^3 + 5*sd^4 + 5/6*sd^5))*(d < r)
  }
  return(V + (d == 0)*(fvar + evar))
}


