#' @title Affinity Propagation
#' @description Runs the Affinity Propagation clustering algorithm of Frey and Dueck, 2007.
#' 
#' @param sim.mat a similarity matrix between individuals to be clustered.
#' @param num.iter maximum number of iterations to attempt.
#' @param stable.iter number of sequential iterations for which consistent clustering is considered acceptable.
#' @param shared.pref type of shared preference to use. Can be one of "min", "median", or a numeric value.
#' @param lambda damping factor.
#' 
#' @return A matrix with one row per sample in 'sim.mat' and one column for each iteration. Values in columns indicate 
#' cluster assignment (arbitrary numbers) for each sample.
#' 
#' @references Frey, B.J., and D. Dueck. 2007. Clustering by passing messages between data points. Science 315:972-976
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(iris)
#' 
#' # Take 75 random iris rows for example
#' iris <- iris[sample(1:nrow(iris), 75), ]
#' iris <- droplevels(iris)
#' 
#' iris.sim <- -dist(iris[, -5])
#' 
#' iris.affin <- affin.prop(iris.sim, stable.iter = 5)
#' table(iris$Species, iris.affin[, ncol(iris.affin)])
#' 
#' @importFrom stats median
#' @export
#' 
affin.prop <- function(sim.mat, num.iter = 100, stable.iter = 10, shared.pref = "min", lambda = 0.5) {
  #
  # Affinity Propagation following :
  #   Frey, B.J., and D. Dueck.  2007.  Clustering by passing messges between data points.  Science 315:972-976.
  #
  # Eric Archer (eric.archer@noaa.gov)
  # 5/16/2007
  #

  # check similarity matrix structure
  if (!is.matrix(sim.mat)) sim.mat <- as.matrix(sim.mat)
  if (nrow(sim.mat) != ncol(sim.mat)) {
    warning("sim.mat is not square - NA returned")
    return(NA)
  } 

  # calculate shared preferences (resp[k,k])
  if (is.character(shared.pref)) {
    if (shared.pref == "median") diag(sim.mat) <- median(sim.mat[lower.tri(sim.mat)])
    if (shared.pref == "min") diag(sim.mat) <- min(sim.mat[lower.tri(sim.mat)])
  }
  if (is.numeric(shared.pref)) diag(sim.mat) <- shared.pref

  # setup matrices
  n <- nrow(sim.mat)
  avail <- resp <- matrix(0, nrow = n, ncol = n)
  max.k.mat <- matrix(NA, nrow = n, ncol = num.iter)

  for (iter in 1:num.iter) {
    # update responsibilities
    prop.resp <- sapply(1:n, function(k) {
      sapply(1:n, function(i) sim.mat[i, k] - max(avail[i, -k] + sim.mat[i, -k]))
    })
    resp <- (resp * lambda) + (prop.resp * (1 - lambda))

    # update availabilities    
    prop.avail <- sapply(1:n, function(k) {
      sapply(1:n, function(i) {
        if (k == i) {
          resp.other.points <- resp[-i, k]
          resp.other.points <- resp.other.points[resp.other.points > 0]
          sum(resp.other.points)
        } else {
          resp.other.points <- resp[-c(i, k), k]
          resp.other.points <- resp.other.points[resp.other.points > 0]
          min(0, resp[k, k] + sum(resp.other.points))
        }
      })
    })
    avail <- (avail * lambda) + (prop.avail * (1 - lambda))

    # find exemplars
    avail.resp.sum <- avail + resp
    max.k <- sapply(1:n, function(i) which.max(avail.resp.sum[i, ]))    

    # check if stable in last 'stable.iter' iterations
    max.k.mat[, iter] <- max.k
    if (iter >= stable.iter) {
      cols <- (iter - stable.iter):iter
      stable <- sapply(1:nrow(max.k.mat), function(row) length(unique(max.k.mat[row, cols])) == 1)
      if (all(stable)) {
        max.k.mat <- matrix(max.k.mat[, 1:iter], nrow = n, ncol = iter)
        break
      }
    }
  }

  rownames(max.k.mat) <- rownames(sim.mat)
  max.k.mat
}

