#-----------------------------------------------------------------------------#
# Calculate matrix of log Likelihood derivatives  
#
# @param integrand function passed to \code{\link{log_likelihood}}. Defaults to
# \code{\link{logit_integrand}}
# @param X covariate matrix
# @param A vector of treatment assignments
# @param G vector of group assignments
# @param fixed.effects vector of fixed effect parameters
# @param random.effects OPTIONAL vector random effect parameters
# @param ... additional arguments passed to \code{integrand} or \code{\link{grad}}. 
# For example, one can change the \code{method} argument in \code{grad}.
# @return N X length(params) matrix of scores
# @export
#-----------------------------------------------------------------------------#

score_matrix <- function(integrand,
                         X, A, G, 
                         fixed.effects,
                         random.effects,
                         ...)
{
  ## Warnings ##
  if(length(fixed.effects) != ncol(X)){
    stop("The length of params is not equal to the number of predictors")
  }
  
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  XX <- cbind(X, A)
  pp <- length(fixed.effects)
  gg <- sort(unique(G))
  
  ## Compute score for each group and parameter ##
  int.args <- append(get_args(integrand, dots),
                     get_args(integrate, dots))
  fargs <- append(int.args, get_args(numDeriv::grad, dots))
  
  print("Calculating matrix of scores...")
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
               args <- append(fargs, 
                              list(integrand = integrand, 
                                   fixed.effects = fixed.effects,
                                   random.effects = random.effects,
                                   A = xx[ , (pp + 1)],
                                   X = xx[ , 1:pp]))
               return(do.call(score_calc, args = args))})

  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = pp + length(random.effects), 
                byrow = TRUE,
                dimnames = list(gg, names(c(fixed.effects, random.effects))))
  
  return(out)
}