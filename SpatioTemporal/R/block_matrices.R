##################################################
## FUNCTIONS THAT COMPUTE BLOCK MATRIC PRODUCTS ##
##################################################
##Functions in this file:
## calc.mu.B        EX:ok
## calc.iS.X        EX:with calc.mu.B
## calc.X.iS.X      EX:with calc.mu.B

##' Computes either the product between a block diagonal, square matrix
##' \code{iS} and a block matrix \code{X}; the quadratic form of a block
##' diagonal, square matrix, \code{t(X)*iS*X}; or a block matrix multiplied by a
##' vector, \code{X*alpha}.
##' 
##' @title Matrix Multiplication with Block Matrices
##' 
##' @param X A list of \code{m} matrices with which to form the block matrix;
##'   each matrix should be \code{p[i]} - by - \code{n}.
##' @param alpha A list of \code{m} vectors, with the i:th vector being 
##'   of length \code{p[i]}.
##' 
##' @return matrix containing iS*X, X'*iS*X, or X*alpha.
##' 
##' @example Rd_examples/Ex_calc_mu_B.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' 
##' @family likelihood utility functions
##' @family block matrix functions
##' @export
calc.mu.B <- function(X, alpha){
  dim <- checkDimInternal(X)
  ##calculate the land use regression for the temporal trends
  mu.B <- matrix(0, dim$n, dim$m)
  for( i in c(1:dim$m) ) ##loop over diff. trends
    mu.B[,i] <- X[[i]] %*% alpha[[i]]
  ##add names
  colnames(mu.B) <- names(X)
  rownames(mu.B) <- rownames(X[[1]])

  return(mu.B)
}##function calc.mu.B


##' @param iS A block diagonal, square matrix, with \code{dm} blocks each of
##'   size \code{n} - by - \code{dn}.
##' @rdname calc.mu.B
##' @export 
calc.iS.X <- function(X, iS){
  dim <- checkDimInternal(X)
  iS.X <- matrix(0, dim(iS)[1], sum(dim$p))
  Ind <- c(0,0)
  ##populate the matrices, taking the block structure into account
  for(i in 1:dim$m){
    iS.X[Ind[1]+(1:dim$n), Ind[2]+(1:dim$p[i])] <-
      iS[Ind[1]+(1:dim$n), Ind[1]+(1:dim$n)] %*% X[[i]]
    ##increase the indecies
    Ind[1] <- Ind[1]+dim$n
    Ind[2] <- Ind[2]+dim$p[i]
  }
  return(iS.X)
}##function calc.iS.X

##' @param iS.X Matrix containing the product of \code{iS} and \code{X}. Output
##'   from \code{calc.iS.X}.
##' @rdname calc.mu.B
##' @export 
calc.X.iS.X <- function(X, iS.X){
  dim <- checkDimInternal(X)
  X.iS.X <- matrix(0, sum(dim$p), sum(dim$p))
  Ind <- c(0,0)
  ##populate the matrices, taking the block structure into account
  for(i in 1:dim$m){
    X.iS.X[Ind[2]+(1:dim$p[i]), Ind[2]+(1:dim$p[i])] <-
      t(X[[i]]) %*% iS.X[Ind[1]+(1:dim$n),Ind[2]+(1:dim$p[i])]
    ##increase the indecies
    Ind[1] <- Ind[1]+dim$n
    Ind[2] <- Ind[2]+dim$p[i]
  }
  return(X.iS.X)
}##function calc.X.iS.X
