#-----------------------------------------------------------------------------#
# Create a matrix of group IP weights 
# 
# Creates a number of groups by number of allocation schemes matrix of group weights.
# Allocation schemes are selected by the user. 
# 
# Groups should be numbered 1, ..., N 
#  
# @param integrand the function used in the weight calculation. Defaults to 
# \code{\link{logit_integrand}}
# @param allocations coverage levels in [0, 1]. Can be vector.
# @param X covariate matrix
# @param A vector of treatment assignments
# @param G vector of group assignments
# @param fixed.effects vector of fixed effect parameters
# @param random.effects OPTIONAL vector random effect parameters
# @param ... additional arguments passed to \code{integrand}
# @return a length(unique(group)) X length(alphas) matrix of group weights 
# @export
#
#-----------------------------------------------------------------------------#

wght_matrix <- function(integrand, 
                        allocations, 
                        X, A, G,
                        fixed.effects,
                        random.effects = NULL,
                        ...)
{
  ## Gather necessary bits ##
  XX <- cbind(X, A)
  p  <- length(fixed.effects)
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  
  ## Compute weight for each group and allocation level ##
  print('Calculating matrix of IP weights...')
  
  w.list <- lapply(aa, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              wght_calc(integrand = integrand, 
                        allocation = allocation, 
                        A = x[, p+1], X = x[, 1:p], 
                        fixed.effects = fixed.effects, 
                        random.effects = random.effects, ...)})
    as.numeric(w)
  }) 
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, aa))
  
  return(w.matrix)
}
