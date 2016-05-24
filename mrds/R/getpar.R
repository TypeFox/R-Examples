#' Extraction and assignment of parameters to vector
#'
#' Extracts parameters of a particular type (scale,
#' shape, adjustments or g0 (p(0))) from the vector of parameters in
#' \code{ddfobj}. All of the parameters are kept in a single vector for
#' optimization even though they have very different uses.  \code{assign.par}
#' parses them from the vector based on a known structure and assigns them into
#' \code{ddfobj}.  \code{getpar} extracts the requested types to be extracted
#' from \code{ddfobj}.
#'
#' @aliases getpar
#' @param ddfobj distance sampling object (see \code{\link{create.ddfobj}})
#' @param fitting character string which is either "all","key","adjust" which
#'   determines which parameters are retrieved
#' @param index logical that determines whether parameters are returned (FALSE)
#'   or starting indices in parameter vector for scale, shape, adjustment
#'   parameters
#' @return index==FALSE, vector of parameters that were requested or
#'   index==TRUE, vector of 3 indices for scale, shape, adjustment
#' @note Internal functions not intended to be called by user.
#' @author Jeff Laake
#' @seealso assign.par
#' @keywords utility
getpar <- function(ddfobj,fitting="all",index=FALSE){

  pars <- NULL

  if(!index){
    if(fitting=="key" | fitting=="all"){
      if(!is.null(ddfobj$shape)) pars <- c(pars,ddfobj$shape$parameters)
      if(!is.null(ddfobj$scale)) pars <- c(pars,ddfobj$scale$parameters)
    }else{
      if(!is.null(ddfobj$shape)){
        pars <- c(pars,rep(NA,length(ddfobj$shape$parameters)))
      }
      if(!is.null(ddfobj$scale)){
        pars <- c(pars,rep(NA,length(ddfobj$scale$parameters)))
      }
    }

    if(fitting=="adjust" | fitting=="all"){
      if(!is.null(ddfobj$adjustment)){
        pars <- c(pars,ddfobj$adjustment$parameters)
      }
    }else{
      if(!is.null(ddfobj$adjustment)){
        pars <- c(pars,rep(NA,length(ddfobj$adjustment$parameters)))
      }
    }

    return(pars)

  }else{
    indices <- rep(0,3)
    if(!is.null(ddfobj$shape)){
      indices[1] <- length(ddfobj$shape$parameters)
    }

    if(!is.null(ddfobj$scale)){
      indices[2] <- indices[1]+length(ddfobj$scale$parameters)
    }

    if(!is.null(ddfobj$adjustment)){
      if(indices[2]!=0){
        indices[3] <- indices[2]+length(ddfobj$adjustment$parameters)
      }else{
        indices[3] <- indices[1]+length(ddfobj$adjustment$parameters)
      }
    }
    return(indices)
  }
}
