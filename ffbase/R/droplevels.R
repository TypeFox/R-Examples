#' The function \code{droplevels} is used to drop unused levels from a \code{ff}
#'  factor or , more commonly, from factors in a \code{ffdf}
#' 
#' @seealso \code{\link{droplevels}} \code{\link{droplevels.ffdf}}
#' @export
#' @export droplevels.ff
#' @method droplevels ff 
#' @param x \code{ff} object
#' @param ... not used
#' @param inplace if \code{TRUE} the columns will be physically changed, 
#' otherwise (default) a new \code{ff} vector will be created
#' @return \code{ff} object where levels of factors are dropped
droplevels.ff <- function(x, ..., inplace=FALSE){
   if (!is.factor.ff(x)){
		stop("droplevels can only applied to a factor")      
   }
   
   if (!inplace){
      x <- ff::clone.ff(x)
   }
   
   levs <- levels(x)
   used <- logical(length(levs))
   for (i in chunk(x)){
     Log$chunk(i)
     used <- used | (levs %in% x[i])
   }
   recodeLevels(x, levs[used])
}


#' The function \code{droplevels} is used to drop unused levels from factors 
#' in a \code{\link{ffdf}}
#' 
#' @seealso \code{\link{droplevels}} \code{\link{droplevels.ff}}
#' @export
#' @export droplevels.ffdf
#' @method droplevels ffdf
#' @param x \code{ffdf} object
#' @param except specify which columns will be excluded from dropping levels
#' @param ... further arguments passed to \code{\link{droplevels.ff}}
#' @param inplace if \code{TRUE} the columns will be physically changed, 
#' otherwise (default) new \code{ff} vectors will be created
#' @return \code{ffdf} object where levels of factors are dropped
droplevels.ffdf <- function(x, except=NULL, ..., inplace=FALSE){
   ffs <- bit::physical(x)
   names(ffs) <- names(x)
   ix <- sapply(ffs, is.factor.ff)
   if (!is.null(except))
        ix[except] <- FALSE
	
   ffs[ix] <- sapply(ffs[ix], droplevels.ff, ..., inplace, simplify=FALSE)
   do.call("ffdf", ffs)
}
