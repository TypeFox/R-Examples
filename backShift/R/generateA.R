#' Generates a connectivity matrix A.
#' 
#' @description  Generates a connectivity matrix A with cycle product smaller 
#'  than 1.
#'
#' @details If \code{expNumNeigh} and \code{maxCoef} are large, function 
#' may fail to find a connectivity matrix with cycle product smaller one. In this
#' case, try to lower these parameters.
#'
#' @param p Number of variables.
#' @param expNumNeigh Expected number of neighbors, to be passed to 
#'  \code{\link[pcalg]{randDAG}}.
#' @param minCoef Minimal edge coefficient. The absolute magnitude of the 
#' coefficients will be sampled uniformly at random from the 
#' range \eqn{[minCoef, maxCoef]}.
#' @param maxCoef Maximal edge coefficient. The absolute magnitude of the 
#' coefficients will be sampled uniformly at random from the 
#' range \eqn{[minCoef, maxCoef]}.
#' @param cyclic If \code{TRUE}, connectivity matrix will contain at least one
#' cycle.
#' @param verbose If \code{TRUE}, comments will be printed.
#' 
#' @return A list with two elements
#' \itemize{
#' \item A Connectivity matrix
#' \item sizeCycle Size of the cycle, if \code{cyclic} was set to \code{TRUE}.
#'}
generateA <- function(p, expNumNeigh, minCoef, maxCoef, cyclic, verbose = FALSE){
  graph.obj <- randDAG(p, expNumNeigh, wFUN=list(runif, min=minCoef, max=maxCoef))
  A <- as(graph.obj, "matrix")
  
  while(sum(A) == 0){ # do not want to get empty graph
    graph.obj <- randDAG(p, expNumNeigh, wFUN=list(runif, min=minCoef, max=maxCoef))
    A <- as(graph.obj, "matrix")
  }
  
  # reverse sign of half the entries
  nz <- which(abs(A) > 0, arr.ind = F)
  reverse.ind <- sample(1:length(nz), ceiling(length(nz)/2))
  for(i in 1:length(reverse.ind)) A[nz[reverse.ind[i]]] <- - A[nz[reverse.ind[i]]]
  
  if(cyclic){
    candidate <- addCycle(A, p, minCoef, maxCoef)
    
    # check whether I-A is invertible and whether A has CP < 1
    Dhat <- t(diag(p) - candidate$A)
    while(is.singular.matrix(Dhat) && !hasCPsmallerOne(Dhat, FALSE)$success){
      candidate <- addCycle(A, p, minCoef, maxCoef)
      if(verbose) cat("I - A is singular or has CP >= 1, regenerating A...")
    }
    
    A <- candidate$A
    sizeCycle <- candidate$sizeCycle
    if(verbose) cat('Graph contains cycle of size', sizeCycle, '\n')
  }else{
    sizeCycle <- NULL
  }
  list(A = A, sizeCycle = sizeCycle)
}