#### Wrapper functions for Doptimize and Defficiencies

#' Wrapper function for OApackage Doptimize function.
#'
#' This function generates
#' a single 2-level design with the specified number of runs and factors. The
#' design is optimized for the criterium that is specified by the
#' parameters \eqn{\alpha_i}{alpha_i}. 
#'
#' The criterium that is optimized is:
#'	\eqn{F = \alpha_1 D + \alpha_2 D_s + \alpha_3  D_1}{F = alpha1*D + alpha2*Ds + alpha3 * D1}
#'
#' Here D is the D-efficiency of the design and Ds and D1 other efficiency measures.
#' For the details on these efficiency measures see the \link{oapackage-package}.
#' The values \eqn{\alpha_i}{alpha_i} are scalar parameters.
#' 
#' 
#' @param N Number of runs
#' @param k Number of factors
#' @param nrestarts Number of restarts to generated an optimal design
#' @param alpha1 Parameter of the optimization function
#' @param alpha2 Parameter of the optimization function
#' @param alpha3 Parameter of the optimization function
#' @param verbose Integer that determines the amount of debug output
#' @param method Integer, default: 0. The method 0 uses updates of single elements of the design matrix. The method 1 uses swaps of 2 elements of the matrix.
#' @param niter Integer (maximum number if iteration steps in the optimization)
#' @param maxtime Float (maximum running time before aborting the optimization)
#' @return A matrix containing the generated design
Doptimize=function(N, k, nrestarts, alpha1=1, alpha2=0, alpha3=0, verbose=1, method=0, niter=100000, maxtime=500) {

nabort <- -1
nn <- N*k
#print('call')
tmp <- .C('DoptimizeR', as.integer(N), as.integer(k), as.integer(nrestarts), as.double(alpha1), as.double(alpha2), as.double(alpha3), as.integer(verbose), as.integer(method), as.integer(niter), as.double(maxtime), as.integer(nabort), result=double(nn) ) 
#print(tmp)
#message('Doptimize: done')
p = tmp[['result']]

A <- array(p, dim=c(N,k) )
A
 }

 #' Wrapper function for OApackage Defficiencies function.
#'
#' This function calculates the D, Ds- and D1-efficiency of a design. The definitions
#' of these efficiencies can be found on the main page of the package \link{oapackage}.
#' 
#' @param A An array containing the design (in 0, 1 format)
#' @return A list containing the calculated efficiencies
Defficiencies=function(A) {

sz <- dim(A)
ndim <- length(sz)
if ( ndim!=2 ) {
print('Defficiencies: input should be a 2-dimensional array')
return
}
N = sz[1]
k = sz[2]
if ( N > 5000 || k > 5000 ) {
  print('Defficiencies: input array should have dimensions smaller than 5000x5000')
  return
}

#message('Defficiencies: call')
tmp <- .C('DefficienciesR', as.integer(N), as.integer(k), as.double(A), D=double(1), Ds=double(1), D1=double(1) ) 
#message('Defficiencies: done')

dd=c( tmp['D'], tmp['Ds'], tmp['D1'])
dd
}

# To compile documentation from the source code:
#
# setwd(...)
# devtools::document() 

# Testing functions
#foo <- function(x){x*2}
#cfoo=function(a,b){.C('cfoo',as.double(a),as.double(b),c=as.double(0))$c}

