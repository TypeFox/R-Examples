#' Numerically integrate pdf of observed distances over specified ranges
#'
#' Computes integral of pdf of observed distances over x for each observation.
#' The method of computation depends on argument switches set and the type of
#' detection function.
#'
#' @param ddfobj distance detection function specification
#' @param select logical vector for selection of data values
#' @param width truncation width
#' @param int.range integration range matrix; vector is converted to matrix
#' @param standardize logical used to decide whether to divide through by the
#'   function evaluated at 0
#' @param point logical to determine if point count (\code{TRUE}) or line
#'   transect (\code{FALSE})
#' @return vector of integral values - one for each observation
#' @author Jeff Laake & Dave Miller
#' @keywords utility
# @importFrom mgcv uniquecombs
#' @importFrom stats integrate
integratepdf <- function(ddfobj, select, width, int.range,
                         standardize=TRUE, point=FALSE){
  # Make sure there is consistency between integration ranges and data
  # It is ok to have a single observation with multiple ranges or a single range
  # with multiple observations but otherwise the numbers must agree if both >1

  if(!is.matrix(int.range)){
    if(is.vector(int.range) && length(int.range)==2){
      int.range <- matrix(int.range,ncol=2,nrow=1)
    }else{
      stop("int.range is not a matrix and cannot be given the required matrix structure")
    }
  }

  ## select tells us which data we compute the integrals for,
  ## if it's null then we compute for all data in ddfobj$xmat
  if(is.null(select)){
    nobs <- nrow(ddfobj$xmat)
    index <- 1:nobs
  }else{
    nobs <- sum(select)
    index <- which(select)
  }

  # number of integration ranges must be 1 or match number of observations
  # OR only be 1 observation and many ranges
  if(nrow(int.range)>1 && nobs>1 && nrow(int.range)!=nobs){
    stop("\n Number of integration ranges (int.range) does not match number of observations\n")
  }

  ## Now compute the integrals

  # if there is only 1 integral to compute (no covariates/1 set of covariates
  # & only one set of integration ranges), that's easy
  if(nobs==1){
    return(gstdint(int.range[1,], ddfobj=ddfobj, index=1, select=NULL,
           width=width, standardize=standardize, point=point, stdint=FALSE))
  }else{
  # if there are multiple covariates or multiple ranges

    # make int.range have nintegrals rows if we want it to
    # already checked above that this is okay
    # this allows us to simplify the code below
    if(nrow(int.range)==1){
      int.range <- t(replicate(nobs,int.range,simplify=TRUE))
    }

    ### find unique observations
    # need unique model matrix-int.range combinations
    # want them within those rows we selected to compute for
    #   we know from above that int.range has either nrow(data) rows or
    #   length(index) rows.
    if(is.null(ddfobj$shape)){
      newdat <- cbind(ddfobj$scale$dm[select, , drop=FALSE], int.range)
    }else{
      if(ncol(ddfobj$shape$dm)>1){
        scale_dm <- ddfobj$shape$dm[select, , drop=FALSE]
        scale_dm[,"(Intercept)"] <- NULL
      }else{
        scale_dm <- NULL
      }
      newdat <- cbind(ddfobj$scale$dm[select, , drop=FALSE],
                      scale_dm, int.range)
    }
    u.rows <- mgcv::uniquecombs(newdat)
    uu.index <- sort(unique(attr(u.rows,"index")))
    u.index <- attr(u.rows,"index")

    # generate the indices that we want to calculate integrals for
    ind <- match(uu.index, u.index)

    # calculate the integrals
    ints <- gstdint(int.range[ind,,drop=FALSE], ddfobj=ddfobj,
                    index=index[ind], select=NULL, width=width,
                    standardize=standardize, point=point,
                    stdint=FALSE)

    ## now rebuild the integrals and populate the return vector
    integrals <- ints[attr(u.rows,"index")]
  }
  return(integrals)
}
