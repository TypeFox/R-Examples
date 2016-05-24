#####################################################
## INTERFACES FOR C CODE - COVARIANCE CONSTRUCITON ##
#####################################################
##Functions in this file:
## makeSigmaB     EX:ok
## makeSigmaNu    EX:ok
## parsCovFuns    EX:ok
## namesCovFuns   EX:ok
## evalCovFuns    EX:ok
## crossDist      EX:ok

##' Function that creates a block covariance matrix with equal sized blocks.
##' Used to construct the Sigma_B matrix.
##'
##' Any parameters given as scalars will be \code{rep}-ed to match
##' \code{length(pars)}.
##' 
##' @title Create Block Covariance Matrix (Equal Block Sizes)
##' @param pars List of parameters for each block; if not a list a single
##'   block matrix is assumed. Should match parameters suggested by
##'   \code{\link{parsCovFuns}}.
##' @param dist Distance matrix.
##' @param type Name(s) of covariance functions, see
##'   \code{\link{namesCovFuns}}.
##' @param nugget Vector of nugget(s) to add to the diagonal of each matrix.
##' @param symmetry \code{TRUE}/\code{FALSE} flag if the \code{dist} is
##'   symmetric, resulting in a symmetric covariance matrix.
##' @param ind2.to.1 Vectors, that for each index along the second dimension
##'   gives a first dimension index, used only if \code{symmetry=FALSE}
##'   to determine which covariances should have an added nugget (collocated
##'   sites).
##' @param sparse If \code{TRUE}, return a block diagonal sparse matrix, see
##'   \code{\link[Matrix:bdiag]{bdiag}}.
##' @return Block covariance matrix of size \code{dim(dist)*n.blocks}.
##' 
##' @example Rd_examples/Ex_makeSigmaB.R
##' 
##' @author Johan Lindström
##' @family block matrix functions
##' @family covariance functions
##' @export
##' @import Matrix
##' @useDynLib SpatioTemporal C_makeSigmaB
makeSigmaB <- function(pars, dist, type="exp", nugget=0,
                       symmetry=dim(dist)[1]==dim(dist)[2],
                       ind2.to.1=1:dim(dist)[2], sparse=FALSE){
  if( !is.list(pars) ){
    ##pars not a list, assuming that we only have one block
    pars <- list(pars)
  }  
  ##First repeat length-one variables to match number of blocks
  n.blocks <- length(pars);
  if( length(type)==1 ){
    type <- rep(type, n.blocks)
  }
  if( length(nugget)==1 ){
    nugget <- rep(nugget, n.blocks)
  }

  ##call C-code, internal error-checking
  tmp <- .Call(C_makeSigmaB, type, pars, nugget, dist,
               as.integer(symmetry), as.integer(ind2.to.1),
               as.integer(sparse))
  if( sparse ){
    if(symmetry){
      tmp <- lapply(tmp, forceSymmetric)
    }
    tmp <- Matrix::bdiag(tmp)
  }
  return( tmp )
}##function makeSigmaB


##' Function that creates a block covariance matrix with unequally
##' sized blocks. Used to construct the Sigma_nu matrix.
##'
##' @title Create Block Covariance Matrix (Unequal Block Sizes)
##' @param pars Vector of parameters, as suggested by
##'   \code{parsCovFuns}.
##' @param dist Distance matrix.
##' @param type Name of the covariance function to use, see
##'   \code{\link{namesCovFuns}}.
##' @param nugget A value of the nugget or a vector of length
##'   \code{dim(dist)[1]} giving (possibly) location specific nuggets.
##' @param random.effect A constant variance to add to the covariance matrix,
##'   can be interpereted as either and partial sill with infinite
##'   range or as a random effect with variance given by \code{random.effect}
##'   for the mean value.
##' @param symmetry \code{TRUE}/\code{FALSE} flag if the \code{dist} matrix is
##'   symmetric. If also \code{ind1==ind2} and \code{blocks1==blocks2} the
##'   resulting covariance matrix will be symmetric.
##' @param blocks1,blocks2 Vectors with the size(s) of each of the
##'   diagonal blocks, usually \code{\link{mesa.model}$nt}. If \code{symmetry=TRUE}
##'   and then \code{blocks2} defaults to \code{blocks1} if missing.
##' @param ind1,ind2 Vectors indicating the location of each element in the
##'   covariance matrix, used to index the \code{dist}-matrix to
##'   determine the distance between locations, usually
##'   \code{\link{mesa.model}$obs$idx}. If \code{symmetry=TRUE}
##'   and then \code{ind2} defaults to \code{ind1} if missing.
##' @param ind2.to.1 Vectors, that for each index along the second dimension,
##'   \code{ind2}, gives a first dimension index, \code{ind1}, used only if
##'   \code{symmetry=FALSE} to determine which covariances should have an
##'   added nugget (collocated sites).
##' @param sparse If \code{TRUE}, return a block diagonal sparse matrix, see
##'   \code{\link[Matrix:bdiag]{bdiag}}.
##' @return Block covariance matrix of size
##'   \code{length(ind1)}-by-\code{length(ind2)}.
##' 
##' @example Rd_examples/Ex_makeSigmaNu.R
##' 
##' @author Johan Lindström
##' @family block matrix functions
##' @family covariance functions
##' @export
##' @import Matrix
##' @useDynLib SpatioTemporal C_makeSigmaNu
makeSigmaNu <- function(pars, dist, type="exp", nugget=0, random.effect=0,
                        symmetry=dim(dist)[1]==dim(dist)[2],
                        blocks1=dim(dist)[1], blocks2=dim(dist)[2],
                        ind1=1:dim(dist)[1], ind2=1:dim(dist)[2], 
                        ind2.to.1=1:dim(dist)[2], sparse=FALSE){
  if( missing(blocks2) && symmetry ){
    blocks2 <- blocks1
  }
  if( missing(ind2) && symmetry ){
    ind2 <- ind1
  }
  ##call C-code, internal error-checking
  tmp <- .Call(C_makeSigmaNu, type, pars, nugget, random.effect, dist,
               as.integer(blocks1), as.integer(blocks2),
               as.integer(ind1), as.integer(ind2), as.integer(ind2.to.1),
               as.integer(symmetry), as.integer(sparse))
  if( sparse ){
    ##order differs agains makeSigmaB,
    ##sigmaNu often has more and smaller blocks -> this option is faster
    tmp <- Matrix::bdiag(tmp)
    if(symmetry && isTRUE(all.equal(blocks1,blocks2)) &&
       isTRUE(all.equal(ind1,ind2)) ){
      tmp <- forceSymmetric(tmp)
    }
  }
  return( tmp )
}##function makeSigmaNu


##' Provides a list of parameter names for the given covariance function(s),
##' excluding the nugget which is added elsewhere.
##'
##' @title Parameter Names for Covariance Function(s)
##' @param type Name(s) of covariance functions, see \code{\link{namesCovFuns}}.
##' @param list Always return a list (if FALSE returns a vector if possible)
##' @return Character vector with parameter names (excluding the nugget),
##'   \code{NULL} if the name is unknown. Returns a list if type contains
##'   more than one element.
##' 
##' @examples
##'   ##all possible parameters
##'   parsCovFuns()
##'   ##just one covariance function
##'   parsCovFuns("exp")
##'   ##non existant covariance function
##'   parsCovFuns("bad.name")
##' 
##' @author Johan Lindström
##' @family covariance functions
##' @export
##' @useDynLib SpatioTemporal C_cov_pars
parsCovFuns <- function(type = namesCovFuns(), list=FALSE){
  if( length(type)==0 ){
    stop("'type' has to be of length>0")
  }
  ##special case for length one type
  if( length(type)==1 && !list){
    return( .Call(C_cov_pars, type[1]) )
  }
  ##o.w. return a list of possible names
  par.names <- vector("list", length(type))
  for(i in 1:length(type)){
    tmp <- .Call(C_cov_pars, type[i])
    if( !is.null(tmp) ){
      par.names[[i]] <- tmp
    }
  }
  names(par.names) <- type
  return( par.names )
}##function parsCovFuns

##' Returns a list of possible covariance function names
##'
##' Available covariance functions (\code{d} is the distance
##' between points):
##' \describe{
##'   \item{\code{exp,exponential}}{Exponential covariance:
##'     \deqn{\sigma^2 \exp(-d/\rho)}{sill * exp( -d/range )} }
##'   \item{\code{exp2,exponential2,gaussian}}{Gaussian/double
##'     exponential covariance:
##'     \deqn{\sigma^2 \exp(-(d/\rho)^2)}{sill * exp( -(d/range)^2 )} }
##'   \item{\code{cubic}}{Cubic covariance: 
##'     \deqn{\sigma^2 (1 - 7 (d/\rho)^2 + 8.75 (d/\rho)^3 -
##'           3.5 (d/\rho)^5 + 0.75 (d/\rho)^7)}{sill*(1 -
##'           7*(d/range)^2 + 8.75*(d/range)^3 - 3.5*(d/range)^5 +
##'           0.75*(d/range)^7)}
##'     if \eqn{d<\rho}{d<range}.}
##'   \item{\code{spherical}}{Spherical covariance: 
##'     \deqn{\sigma ^2(1 - 1.5(d/\rho) + 0.5 (d/\rho)^3)}{sill * (1 -
##'           1.5*(d/range) + 0.5*(d/range)^3)}
##'     if \eqn{d<\rho}{d<range}.}
##'   \item{\code{matern}}{Matern covariance: 
##'     \deqn{\frac{\sigma^2}{\Gamma(\nu) 2^{\nu-1}}
##'           \left(\frac{d\sqrt{8\nu}}{\rho}\right)^\nu
##'           K_\nu\left(\frac{d\sqrt{8\nu}}{\rho}\right)}{sill /
##'           (gamma(shape)*2^(shape-1)) *
##'           (d*sqrt(8*shape)/range)^shape *
##'           besselK( (d*sqrt(8*shape)/range), shape ) } }
##'   \item{\code{cauchy}}{Cauchy covariance: 
##'     \deqn{ \frac{\sigma^2}{(1 + (d/\rho)^2)^{\nu}}}{sill *
##'           (1 + (d/range)^2) ^ -shape)} }
##'   \item{\code{iid}}{IID covariance, i.e. zero matrix 
##'     since nugget is added afterwards. \deqn{0}{0}}
##' }
##' 
##' @title Available covariance functions
##' @return Character vector with valid covariance function names.
##' 
##' @examples
##'   namesCovFuns()
##' 
##' @author Johan Lindström
##' @family covariance functions
##' @export
##' @useDynLib SpatioTemporal C_cov_names
namesCovFuns <- function(){
  ##call C-function that returns a list of available covariance functions
  return( .Call(C_cov_names) )
}##function namesCovFuns


##' Computes covariance functions (excluding nugget) for a given vector or
##' matrix of distances.
##'
##' @title Compute Covariance Function
##' @param type Name of covariance functions, see \code{\link{namesCovFuns}}.
##' @param pars Parameter for the covariance function, see
##'   \code{\link{parsCovFuns}}.
##' @param d Vector/matrix for which to compute the covariance function.
##' @return Covariance function computed for all elements in d.
##' 
##' @example Rd_examples/Ex_evalCovFuns.R
##' 
##' @author Johan Lindström
##' @family covariance functions
##' @export
##' @useDynLib SpatioTemporal C_cov_simple
evalCovFuns <- function(type="exp", pars=c(1,1), d=seq(0,10,length.out=100)){
  ##call C-code, internal error-checking
  cov <- .Call(C_cov_simple, type, pars, d)
  if( is.array(d) ){
    cov <- array(cov, dim(d))
  }else{
    cov <- as.vector(cov)
  }
  return( cov )
}##function evalCovFuns

##' Computed the Euclidian distance matrix between to sets of points.
##'
##' @title Computed the Euclidian Distance Matrix
##' @param coord1,coord2 Matrices with the coordinates of locations, between
##'   which distances are to be computed.
##' @return A \code{dim(coord1)[1]}-by-\code{dim(coord2)[1]} distance matrix.
##' 
##' @example Rd_examples/Ex_crossDist.R
##' 
##' @author Johan Lindström
##' @family covariance functions
##' @family basic linear algebra
##' @export
##' @useDynLib SpatioTemporal C_dist
crossDist <- function(coord1, coord2=coord1){
  if( missing(coord2) ){
    symmetric <- TRUE
  }else{
    symmetric <- FALSE
  }
  ##cast to double
  coord1 <- as.matrix(coord1)
  storage.mode(coord1) <- "double"
  coord2 <- as.matrix(coord2)
  storage.mode(coord2) <- "double"
  
  ##call C-code, internal error-checking
  return( .Call(C_dist, coord1, coord2, as.integer(symmetric)) )
}
