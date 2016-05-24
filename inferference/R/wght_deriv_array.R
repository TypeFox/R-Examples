#-----------------------------------------------------------------------------#
# Create an array of group weight derivatives 
# 
# Uses \code{\link{wght_deriv_calc}} to compute the weight derivatives for each 
# group per coverage level
#
# @param integrand the function to passed to the argument 'f' of \code{\link{integrate}},
# which is part of \code{\link{wght_calc}}. 
# @param allocations coverage levels in (0, 1), possibly (probably) vector valued
# @param X covariate matrix
# @param A vector of treatment assignments
# @param G vector of group assignments
# @param fixed.effects vector of fixed effect parameters
# @param random.effects OPTIONAL vector random effect parameters
# @param ... additional arguments passed to integrand
# @return a length(unique(group)) X length(params) X length(alphas) array of 
# group weight derivatives
# @export
# 
#-----------------------------------------------------------------------------#

wght_deriv_array <- function(integrand, 
                             allocations, 
                             X, A, G,
                             fixed.effects,
                             random.effects,
                             ...)
{
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  XX <- cbind(X, A)
  p  <- length(fixed.effects) # number of predictors
  pp <- p + length(random.effects)
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  k  <- length(allocations) 
  N  <- length(unique(G))
  dots <- list(...)

  ## Warnings ##

  ## Compute weight (derivative) for each group, parameter, and alpha level ##
  print('Calculating array of IP weight derivatives...')
  
  w.list <- lapply(aa, function(allocation){
    w <- by(XX, INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              wght_deriv_calc(integrand = integrand, 
                              allocation = allocation, 
                              A = x[, p+1], 
                              X = x[, 1:p], 
                              fixed.effects = fixed.effects, 
                              random.effects = random.effects,
                              ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = pp, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, pp, k),
               dimnames = list(gg, names(c(fixed.effects, random.effects)), aa))
  
  return(out)
}