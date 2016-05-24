#' Phylogenetic weights matrix
#' 
#' This function can be used to compute a phylogenetic weights matrix
#' with different methods.
#' 
#' @param tree a \code{phylo}, \code{phylo4} or \code{phylo4d} object.
#' @param dist.phylo a character string specifying the method used to compute phylogenetic distances.
#' Available distances are "\code{patristic}","\code{nNodes}","\code{Abouheif}" and "\code{sumDD}".
#' See Details.
#' @param method a method to compute phylogenetic weights from phylogenetic distances.
#' Available methods are "\code{lag-norm}", "\code{clade}", "\code{inverse}" and "\code{exponential}".
#' See Details.
#' @param mu a numeric value giving the mean of the distribution if \code{method} is \code{lag-norm}.
#' This is a phylogenetic distance.
#' @param sigma a numeric value giving the standard deviation
#' of the distribution if \code{method} is \code{lag-norm}.
#' @param dmax the maximum phylogenetic distance to use to delineate clades.
#' @param alpha a numeric value giving the exponent to use if \code{method} is \code{inverse}.
#' @param beta a numeric value giving the factor to use if \code{method} is \code{exponential}.
#' 
#' @details
#' Method "\code{inverse}": \deqn{\frac{1}{d^{\alpha}}}{w = 1/d^alpha}
#' 
#' The phylogenetic distance matrix is computed internally
#' using the function \code{\link[adephylo]{distTips}} from the package \pkg{adephylo}.
#' See \code{\link[adephylo]{distTips}} for details about the methods.
#' 
#' @return A square matrix of phylogenetic weights whose sums of rows is 1.
#' 
#' @seealso \code{\link[adephylo]{proxTips}} in \pkg{adephylo}.
#' @export
phyloWeights <- function(tree, dist.phylo = "patristic", method = "lag-norm",
                         mu = 0, sigma = 5, dmax = 10, alpha = 1, beta = 1){
  if (inherits(tree, "phylo4")){
    tree <- as(tree, "phylo")
  }
  if (!inherits(tree, "phylo")){
    stop("x has to be a phylo, phylo4 or phylo4d object")
  }
  dist.phylo <- match.arg(dist.phylo, c("patristic", "nNodes", "Abouheif", "sumDD"))
  method <- match.arg(method, c("lag-norm", "clade", "inverse", "exponential"))
  
  d <- distTips(tree, method = dist.phylo)
  d <- as.matrix(d)
  
  if(method == "lag-norm"){
    w <- dnorm(d, mean = mu, sd = sigma)
  }
  
  if(method == "clade"){
    w <- ifelse(d < dmax, 1, 0)
    w[rowSums(w) == 1, ] <- NA
  }
  
  if(method == "inverse"){
    w <- 1 / (d ^ alpha)
  }
  
  if(method == "exponential"){
    w <- exp(-beta * d)
  }
  
  diag(w) <- 0
  w <- prop.table(w, 1)
  return(w)
}