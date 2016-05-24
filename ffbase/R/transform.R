#' Transform a ffdf data.frame 
#' 
#' Same functionality as \code{\link{transform}}, but on a \code{ffdf} object. Please note that you should write
#' your expression as if it is a normal \code{data.frame}. The resulting data.frame
#' however will be a \code{ffdf} data.frame.
#' @usage 
#' \method{transform}{ffdf}(`_data`, \dots)
#'
#' @export
#' @export transform.ffdf
#'
#' @example ../examples/transform.R
#' @param _data \code{\link{ffdf}} data object to be transformed.
#' @param ... named arguments that will be added to the \code{ffdf} data.frame
#' @return a modified clone of \code{`_data`}.
transform.ffdf <- function(`_data`, ...){
    expr <- substitute(list(...))
    parent <- parent.frame()
    
    #chunks <- chunk(`_data`, by=2) #debug chunking
    chunks <- chunk(`_data`)    
    cdat <- `_data`[chunks[[1]],,drop=FALSE]
    e <- eval(expr, cdat, parent)
    cdat[names(e)] <- e
    res <- as.ffdf(cdat)
    
    rownames(res) <- NULL
    nrow(res) <- nrow(`_data`)
    rownames(res) <- rownames(`_data`)
    
    for (i in chunks[-1]){
      Log$chunk(i)
      cdat <- `_data`[i,,drop=FALSE]
      e <- eval(expr, cdat, parent)
       
      cdat[names(e)] <- e
      res[i,] <- cdat
    }
    res
}