######################################################
## INTERFACES FOR C CODE - SPARSE F TIMES SOMETHING ##
######################################################
##Functions in this file:
## calc.tFX    EX:OK
## calc.FX     EX:OK
## calc.tFXF   EX:OK
## calc.FXtF2  EX:OK
## expandF     EX:OK

##' Computes the matrix products between the transpose of a sparse matrix \code{F}
##' containing temporal trends the and a vector/matrix.
##' \cr\cr
##' See the examples for details.
##'
##' @title Compute Matrix Product Bewteen Temporal Trends and a Matrix
##' @param F A (number of obs.) - by - (number of temporal trends) matrix
##'   containing the temporal trends. Usually \code{\link{mesa.model}$F}, where
##'   \code{\link{mesa.model}} is obtained from
##'   \code{\link{createSTmodel}}.
##' @param X A vector or matrix; needs to be a multiple of \code{dim(F)[1]}.
##' @param loc.ind A vector indicating which location each row in \code{F}
##'   corresponds to, usually \cr \code{\link{mesa.model}$obs$idx}.
##' @param n.loc Number of locations.
##' @return Returns a matrix of size \code{n.loc*dim(F)[2]}-by-code{n.x}.
##' 
##' @example Rd_examples/Ex_calc_tFX.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' @family block matrix functions
##' @family temporal trend functions
##' @export
##' @useDynLib SpatioTemporal C_calc_tFX
calc.tFX <- function(F, X, loc.ind, n.loc=max(loc.ind)){
  ##call the c-function, error checking in C-code
  .Call(C_calc_tFX, X, F, as.integer(n.loc), as.integer(loc.ind))
}##function calc.tFX

##' Computes the matrix products between a sparse matrix \code{F}
##' containing the temporal trends and a list of land-use-regression components.
##' \cr\cr
##' See the examples for details.
##'
##' @title Compute Matrix Product Bewteen Temporal Trends and a LUR components
##' @param F A (number of obs.) - by - (number of temporal trends) matrix
##'   containing the temporal trends. Usually \code{\link{mesa.model}$F}, where
##'   \code{\link{mesa.model}} is obtained from
##'   \code{\link{createSTmodel}}.
##' @param LUR A list of matrices, usually \code{\link{mesa.model}$X}. Each matrix
##'   in the list should have the same number of rows, but the number of columns
##'   may vary.
##' @param loc.ind A vector indicating which location each row in \code{F}
##'   corresponds to, usually \cr \code{\link{mesa.model}$obs$idx}.
##' @return Returns a matrix
##' 
##' @example Rd_examples/Ex_calc_FX.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' @family block matrix functions
##' @family temporal trend functions
##' @export
##' @useDynLib SpatioTemporal C_calc_F_part_X
calc.FX <- function(F, LUR, loc.ind){
  ##number of LUR components
  m <- length(LUR)
  ##size of each LUR component
  p <- sapply(LUR, function(x){ dim(x)[2] })
  ##number of locations for each LUR component
  n.loc <- sapply(LUR, function(x){ dim(x)[1] })
  ##check dimensions
  if( any(n.loc != n.loc[1]) ){
    stop( sprintf("In 'calc.FX': All components of 'LUR' does not have the same no. of locations (1 not matching %s).", paste(which(n.loc!=n.loc[1]), collapse=", ") ))
  }
  ##all n-loc equal, use the first one
  n.loc <- n.loc[1]
  FX <- matrix(0, dim(F)[1], sum(p))
  Ind <- 0
  for(i in 1:m){
    ##call the c-function, some error checking in C-code
    tmp <- .Call(C_calc_F_part_X, LUR[[i]], F[,i], as.integer(loc.ind),
                 as.integer(p[i]))
    FX[,(Ind+1):(Ind+p[i])] <- tmp
    Ind <- Ind + p[i]
  }
  return(FX)
}##function calc.FX

##' Computes the quadratic form between a sparse matrix \code{F} containing the
##' temporal trends and the covariance matrix for the residual fields (Sigma_nu)
##' \cr\cr
##' See the examples for details.
##' 
##' @title Compute Quadratic Form Bewteen Temporal Trends and Sigma nu
##' @param F A (number of obs.) - by - (number of temporal trends) matrix
##'   containing the temporal trends. Usually \code{\link{mesa.model}$F}, where
##'   \code{\link{mesa.model}} is obtained from
##'   \code{\link{createSTmodel}}.
##' @param mat A block diagonal, square matrix.
##' @param loc.ind A vector indicating which location each row in \code{F}
##'   corresponds to, usually \cr \code{\link{mesa.model}$obs$idx}.
##' @param n.blocks Number of diagonal blocks in \code{mat} (or \code{R}).
##'   Defaults to 1 (i.e. a full matrix) if neither \code{n.blocks} nor
##'   \code{block.sizes} given, o.w. it defaults to  \code{length(block.sizes)}).
##' @param block.sizes A vector of length \code{n.blocks} with the size of each
##'   of the diagonal blocks. If not given it will assume equal size blocks.
##' @param n.loc Number of locations.
##' @return Returns a square matrix with side \code{dim(F)[2]*n.loc}
##' 
##' @example Rd_examples/Ex_calc_tFXF.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' @family block matrix functions
##' @family temporal trend functions
##' @export
##' @useDynLib SpatioTemporal C_calc_tFXF
calc.tFXF <- function(F, mat, loc.ind, n.blocks=1,
                      block.sizes=rep(dim(mat)[1]/n.blocks,n.blocks),
                      n.loc=max(loc.ind)){
  ##ensure that block sizes are integer valued
  block.sizes <- round(block.sizes)
  ##call the c-function, error checking in C-code
  .Call(C_calc_tFXF, mat, F, as.integer(n.loc),
        as.integer(loc.ind), as.integer(block.sizes))
}##function calc.tFXF

##' Computes the quadratic form  between a sparse matrix \code{F} containing the
##' temporal trends and the covariance matrix for the beta fields (Sigma_B). Or
##' possibly the product between two different \code{F}'s and a cross-covariance
##' matrix.
##' \cr\cr
##' See the examples for details.
##' 
##' @title Compute Quadratic Form Bewteen Temporal Trends and Sigma B
##' @param F,F2 (number of obs.) - by - (number of temporal trends) matrices
##'   containing the temporal trends. Usually \code{\link{mesa.model}$F}, where
##'   \code{\link{mesa.model}} is obtained from
##'   \code{\link{createSTmodel}}.
##' @param mat A block diagonal, with equal size blocks. The
##'   number of blocks need to equal \code{dim(F)[2]}
##' @param loc.ind,loc.ind2 A vector indicating which location each row in \code{F}
##'   corresponds to, usually \cr \code{\link{mesa.model}$obs$idx}.
##' @return Returns a square matrix with side \code{dim(F)[1]}
##' 
##' @example Rd_examples/Ex_calc_FXtF2.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' @family block matrix functions
##' @family temporal trend functions
##' @export
calc.FXtF2 <- function(F, mat, loc.ind, F2=F, loc.ind2=loc.ind){
  ##check dimensions
  if( dim(F)[2] != dim(F2)[2] ){
    stop("In 'calc.FXtF2': dim(F)[2] != dim(F2)[2]")
  }
  if( (dim(mat)[1] %% dim(F)[2]) != 0 ){
    stop("In 'calc.FXtF2': dim(mat)[1] needs to be multiple of dim(F)[2]")
  }
  if( (dim(mat)[2] %% dim(F2)[2]) != 0 ){
    stop("In 'calc.FXtF2': dim(mat)[2] needs to be multiple of dim(F2)[2]")
  }
  ##compute size of the blocks in mat
  block.size <- c(round(dim(mat)[1]/dim(F)[2]),
                  round(dim(mat)[2]/dim(F2)[2]))
  
  ##split mat into a list
  mat.list <- vector( "list", dim(F)[2])
  for(i in 1:length(mat.list)){
    mat.list[[i]] <- mat[(1:block.size[1])+(i-1)*block.size[1],
                         (1:block.size[2])+(i-1)*block.size[2],drop=FALSE]
  }
  ##compute F*mat (using block structure in mat)
  FX <- calc.FX(F, mat.list, loc.ind)

  ##compute (F*mat) %*% F2', using sparse matrices
  F2 <- expandF(F2, loc.ind2, n.loc=block.size[2])
  return( as.matrix(FX %*% t(F2)) )
}##function calc.FXtF2


##' Expands the temporal trends in F to a full matrix (with lots of zeros).
##' Mainly used for testing, and illustration in examples.
##'
##' @title Expand F
##' @param F A (number of obs.) - by - (number of temporal trends) matrix
##'   containing the temporal trends. Usually \code{\link{mesa.model}$F}, where
##'   \code{\link{mesa.model}} is obtained from
##'   \code{\link{createSTmodel}}.
##' @param loc.ind A vector indicating which location each row in \code{F}
##'   corresponds to, usually \cr \code{\link{mesa.model}$obs$idx}.
##' @param n.loc Number of locations.
##' @param sparse Should the returned matrix be sparse (uses the Matrix-package,
##'   see \code{\link[Matrix:sparseMatrix]{sparseMatrix}})
##' @return Returns the expanded F, a \code{dim(F)[1]}-by-\code{n.loc*dim(F)[2]}
##'   matrix
##' 
##' @example Rd_examples/Ex_expandF.R
##' 
##' @author Johan Lindström and Adam Szpiro
##' @family temporal trend functions
##' @export
##' @import Matrix
expandF <- function(F, loc.ind, n.loc=max(loc.ind), sparse=TRUE){
  ##check dimensions
  if( dim(F)[1] != length(loc.ind) ){
    stop("In 'expandF': dim(F)[1] != length(loc.ind)")
  }
  if( max(loc.ind) > n.loc ){
    stop("In 'expandF': max(loc.ind) > n.loc")
  }
  ##create large return matrix, using sparse
  i <- rep(1:dim(F)[1], dim(F)[2])
  j <- rep(loc.ind, dim(F)[2]) + rep((0:(dim(F)[2]-1))*n.loc,
                                     each=dim(F)[1])
  Fexp <- sparseMatrix(i, j, x=c(F), dims=c(dim(F)[1], dim(F)[2]*n.loc))
  if( !sparse ){
    ##cast to matrix if sparse return unwanted.
    Fexp <- as.matrix(Fexp)
  }##if(!sparse)
  ##return expanded matrix
  return(Fexp)
}##function expandF
