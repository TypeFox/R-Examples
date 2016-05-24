##################################################
## INTERFACES FOR C CODE - BASIC LINEAR ALGEBRA ##
##################################################
##Functions in this file:
## makeCholBlock  EX:ok
## invCholBlock   EX:with makeCholBlock
## solveTriBlock  EX:with solveTriBlock
## blockMult      EX:ok
## sumLogDiag     EX:ok
## sumLog         EX:with sumLogDiag
## norm2          EX:with dotProd
## dotProd        EX:ok


##' Provides block diagonal version of the base package functions
##' \code{\link[base:chol]{chol}}, \code{\link[base:chol2inv]{chol2inv}}, and
##' \code{\link[base:backsolve]{backsolve}}.
##' \cr\cr
##' Computes the Cholesky factor, the matrix inverse and solves matrix equation
##' systems for block diagonal matrices.
##' 
##' \code{makeCholBlock} computes the Cholesky factor of a block diagonal matrix
##' using the block diagonal structure to speed up computations.
##' 
##' \code{invCholBlock} uses the Cholesky factor from \code{makeCholBlock} to
##' compute the inverse of \code{mat}.
##' 
##' \code{solveTriBlock} solves equation systems based on the Cholesky factor,
##' using the block diagonal structure to speed up computations (c.f.
##' \code{\link{backsolve}}). The function solves equations of the form R*x = B,
##' and R'*x = B with respect to x, where the transpose is controlled by the
##' parameter \code{transpose}. Aplying the function twice solves mat*x=B, see
##' the examples.
##' 
##' For all three functions the block diagonal structure of the matrix is
##' defined by two input variables, the number of blocks \code{n.blocks}, and
##' the size of each block \code{block.sizes}. The size of the matrices must
##' match the total number of blocks, i.e. \code{sum(block.sizes)} \emph{must}
##' equal \code{dim(mat)}.
##' 
##' The functions can be used for full matrices by setting the number of blocks
##' to 1.
##' 
##' @title Computations for Block Diagonal Matrices
##' @param mat A block diagonal, square, positive definite matrix.
##' @param R Upper right block diagonal Cholesky factor. The output from
##'   \code{\link[base:chol]{chol}} or \cr \code{makeCholBlock}.
##' @param n.blocks Number of diagonal blocks in \code{mat} (or \code{R}).
##'   Defaults to 1 (i.e. a full matrix) if neither \code{n.blocks} nor
##'   \code{block.sizes} given, o.w. it defaults to  \code{length(block.sizes)}).
##' @param block.sizes A vector of length \code{n.blocks} with the size of each
##'   of the diagonal blocks. If not given it will assume equal size blocks.
##' @return \code{makeCholBlock} gives the Cholesky factor and
##'   \code{invCholBlock} gives the inverse of the matrix \code{mat}.
##'   \code{solveTriBlock} gives to answer to the equation system.
##' 
##' @example Rd_examples/Ex_makeCholBlock.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' 
##' @family basic linear algebra
##' @family block matrix functions
##' @export
##' @useDynLib SpatioTemporal C_make_chol_block
makeCholBlock <- function(mat, n.blocks=1,
                          block.sizes=rep(dim(mat)[1]/n.blocks,n.blocks)){
  ##special case for one block
  if( length(block.sizes)==1 ){
    return( chol(mat) )
  }
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##call the c-function
  tmp <- .Call(C_make_chol_block, as.integer(block.sizes), mat)
  ##check if matrix is pos.def
  if( tmp[1]==-1 ){
    stop("In 'makeCholBlock': Matrix not positive definite.")
  }
  return(tmp)
}##function makeCholBlock

##' @rdname makeCholBlock
##' @export
##' @useDynLib SpatioTemporal C_inv_chol_block
invCholBlock <- function(R, n.blocks=1,
                         block.sizes=rep(dim(R)[1]/n.blocks,n.blocks)){
  ##special case for one block
  if( length(block.sizes)==1 ){
    return( chol2inv(R) )
  }
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  
  ##call the c-function, error checking in C-code.
  .Call(C_inv_chol_block, as.integer(block.sizes), R)
}##function invCholBlock

##' @param B Vector or matrix containg the right hand side of the equations
##'   system to be solved; needs to be a multiple of \code{dim(R)[1]}.
##' @param transpose Transpose \code{R} before solving the equation system.
##'   Controls if we solve the equations system given by R*x = B or R'*x=B.
##' @rdname makeCholBlock
##' @export
##' @useDynLib SpatioTemporal C_solve_tri_block
solveTriBlock <- function(R, B, n.blocks=1,
                          block.sizes=rep(dim(R)[1]/n.blocks,n.blocks),
                          transpose=FALSE){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##run in R if length(block.sizes)==1
  if( length(block.sizes)==1 ){
    n.x <- length(B)/dim(R)[1]
    return( backsolve(R, as.matrix(B,block.sizes,n.x), transpose=transpose) )
  }else{
    ##call the c-function, error checking in C-code.
    .Call(C_solve_tri_block, as.integer(block.sizes),
          as.integer(transpose), R, B)
  }
}##function solveTriBlock

##' Computes the matrix product between a block diagonal square matrix and a
##' column vector (or matrix).
##'
##' @title Multiplication of Block Diagonal Matrix and Vector
##' @param mat A block diagonal, square matrix.
##' @param X Vector or matrix to multiply by \code{mat}; \code{length(X)}
##'   needs to be a multiple of \code{dim(mat)[1]}.
##' @param n.blocks Number of diagonal blocks in \code{mat} (or \code{R}).
##'   Defaults to 1 (i.e. a full matrix) if neither \code{n.blocks} nor
##'   \code{block.sizes} given, o.w. it defaults to  \code{length(block.sizes)}).
##' @param block.sizes A vector of length \code{n.blocks} with the size of each
##'   of the diagonal blocks. If not given it will assume equal-sized blocks.
##' @return Returns mat * X.
##' 
##' @example Rd_examples/Ex_blockMult.R
##' 
##' @author Johan Lindström
##' @family basic linear algebra
##' @family block matrix functions
##' @export
##' @useDynLib SpatioTemporal C_block_mult
blockMult <- function(mat, X, n.blocks=1,
                       block.sizes=rep(dim(mat)[1]/n.blocks,n.blocks)){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##call the c-function, error checking in C-code.
  .Call(C_block_mult, as.integer(block.sizes), mat, X)
}##function blockMult

##' Computes the sum of the logarithm of the diagonal elements in a matrix, or
##' of elements in a vector. This corresponds to the logarithm of the
##' determinant for a Cholesky factor.  Behaviour is undefined for any elements
##' that are <=0. 
##'
##' @title Sum the Logarithm of (Diagonal) Elements
##' @param mat A square matrix (preferably a Cholesky factor).
##' @return Sum of the logarithm of the (diagonal) elements.
##' 
##' @example Rd_examples/Ex_sumLogDiag.R
##' 
##' @author Johan Lindström
##' @family basic linear algebra
##' @export
##' @useDynLib SpatioTemporal C_sum_log_diag
sumLogDiag <- function(mat){
  ##call the c-function, error checking in C-code.
  .Call(C_sum_log_diag, mat)
}##function sumLogDiag

##' @rdname sumLogDiag
##' @param v A vector
##' @export
##' @useDynLib SpatioTemporal C_sum_log
sumLog <- function(v){
  ##call the c-function, error checking in C-code.
  .Call(C_sum_log, v)
}##function sumLog


##' @rdname dotProd
##' @export
##' @useDynLib SpatioTemporal C_norm2
norm2 <- function(v1){
  ##calculate sum of the squared elements of a vector (or matrix) (|x|^2)
  ##input vector and vector size
  .Call(C_norm2, v1)
}##function norm2

##' \code{dotProd} computes the inner (or dot/scalar) product between two
##' vectors.
##' \cr\cr
##' \code{norm2} computes the squared 2-norm of all the elements in a matrix or
##' vector.
##' \cr\cr
##' If the vectors are of unequal length \code{dotProd} will give a warning and
##' then truncates the longer vector, discarding any excess elements before the
##' computations.
##'
##' @title Computes Inner Product and Squared 2-norm
##' @param v1,v2 Two vectors
##' @return \code{dotProd} returns the inner product of v1 and v2. \code{norm2}
##'   returns the squared 2-norm of all elements in v1.
##' 
##' @example Rd_examples/Ex_dotProd.R
##' 
##' @author Johan Lindström
##' @family basic linear algebra
##' @export
##' @useDynLib SpatioTemporal C_dotProd
dotProd <- function(v1,v2){
  if(length(v1)!=length(v2)){
    warning("In 'dotProd': vectors of unequal length, truncating the longer vector")
  }
  ##no need to check dimensions
  .Call(C_dotProd, v1, v2)
}##function dotProd
