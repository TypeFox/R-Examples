#' Update the Lambda parameters of clusters.
#' 
#' Updates the Lambda parameters to maximize the likelihood of the data under
#' Mallows' model.
#' 
#' 
#' @param r Matrix of partial rankings.
#' @param R Current modal sequences.
#' @param z Current probabilities of memberships in each cluster.
#' @param G Number of modal sequences.
#' @param dists.to.Rg Matrix of the distances between the data and the current
#' modal sequences.
#' @param dists.table Table of the distance distribution in N! space, under
#' Kendall's metric.
#' @param top.bound The maximum value for the lambda parameter.
#' @return Vector of new lambda parameters for the clusters.
#' @author Erik Gregory
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords lambda maximization
UpdateLambda <-
function(r, R, z, G, dists.to.Rg, 
                         dists.table, top.bound = 1000) {
  lambda <- 0*(1:G)
  # Lambda is the solution to an equation.
  # We update Lambda for each cluster.
  rhs <- 0*(1:G)
  flower1 <- C_lam(0, dists.table = dists.table)*sum(dists.table*as.numeric(names(dists.table)))
  fupper1 <- C_lam(top.bound, 
                   dists.table = dists.table)*sum(as.numeric(names(dists.table))*dists.table*exp(-top.bound*as.numeric(names(dists.table))))
  for (i in 1:G) {
    # Calculate the RHS of the equation needed to 
    # determine lambda values.
    rhs[i] <- sum(z[, i]*dists.to.Rg[, i])/sum(z[, i])
    # Find the root of the Lambda function, to 
    # determine lambda.
    flower  <- flower1 - rhs[i]
    fupper <- fupper1 - rhs[i]

    if (sign(flower) != sign(fupper)) {
      lambda[i] <- uniroot(Lambda, interval = c(0, top.bound), 
                           rhs = rhs[i], 
                           dists.table = dists.table, tol = 1E-5,
                           f.lower = flower, f.upper = fupper)$root
    }
    else {
      print("Solution Exceeded top bound, defaulting to top bound.")
      lambda[i] <- top.bound
    }
  }
  return(lambda)
}
