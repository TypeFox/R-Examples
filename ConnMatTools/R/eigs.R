#' Compute some eigenvalues of a matrix
#' 
#' This function computes a limited number of eigenvalues and eigenvectors of a 
#' matrix. It uses \code{\link[igraph]{arpack}} function from the 
#' \link[igraph]{igraph} package. If this package is not available, it will use 
#' the standard \code{\link{eigen}} function to do the calculation, but will 
#' issue a warning.
#' 
#' @param M a matrix.
#' @param nev number of eigenvalues and eigenvectors to return
#' @param sym A boolean indicating if matrix is symmetric or not.  Defaults to 
#'   checking if this is the case or not.
#' @param which A character string indicating which eigenvalues to return. 
#'   Defaults to "LM", meaning largest magnitude eigenvalues. If not using 
#'   \code{\link[igraph]{arpack}}, then "SM" is also a possibility to return the
#'   smallest magnitude eigenvalues. If using \code{\link[igraph]{arpack}}, then
#'   a number of options are possible, though they are not all guaranteed to 
#'   work for all use cases. See that function for more details.
#' @param use.arpack Boolean determining if calculation is to be done with 
#'   \code{\link[igraph]{arpack}} function from the \link[igraph]{igraph} 
#'   package. This is much quicker for large matrices, but requires 
#'   \link[igraph]{igraph}. Defaults to TRUE, but will use eigen instead if 
#'   \link[igraph]{igraph} is not found.
#' @param options.arpack Additional options for \code{\link[igraph]{arpack}}. 
#'   See that function for details.  Not all options are compatible with this 
#'   function.
#'   
#' @return A list with at least the following two items:
#'   
#'   \item{values}{A set of eigenvalues}
#'   
#'   \item{vectors}{A matrix of eigenvectors}
#'   
#' @seealso See also \code{\link[igraph]{arpack}}
#'   
#' @author David M. Kaplan \email{dmkaplan2000@@gmail.com}
#' @export
#' @include eigs.R
eigs <- function( M, nev=min(dim(M)[1]-1,1),
                  sym=sum(abs(M-t(M)))/sum(abs(M))<1e-10,
                  which="LM",
                  use.arpack=TRUE,
                  options.arpack=NULL ) {
  if (class(M) != "matrix")
    stop("Input M must be a matrix.")
  
  n = dim(M)[1]
  
  if (nev > n) 
    stop("nev must be less than dimension of M")
  
  if (use.arpack && require(igraph)) {
    options.arpack$which = which
    options.arpack$n = n
    options.arpack$nev = nev
    options.arpack$ncv = nev+2
    
    #browser()
    ff = function(x,extra) { extra %*%x }
        
    v = arpack(ff, extra=M, sym=sym, options=options.arpack)
  } else {
    warning("Using eigen. May be slow.  Use arpack in igraph package to improve speed.")
    v = eigen(M,symmetric=sym)
    
    I = switch(which,
               LM = 1:nev,
               SM = n:(n-nev+1),
               stop("which can only be 'LM' or 'SM' when using eigen")
               )
    
    v$values = v$values[I]
    v$vectors = v$vectors[,I]
  }
  
  return(v)
}
