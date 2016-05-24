#' Value Matching for ff objects
#'
#' \code{ffmatch} returns an ff vector of the positions of (first) matches of its first argument in its second. 
#' Similar as \code{\link[base]{match}}. \cr \cr
#' \code{ffdfmatch} allows to match ffdf objects by \code{\link[base]{paste}}-ing together the columns of the ffdf and matching on the pasted column and
#' returns an ff vector of the positions of (first) matches of its first argument in its second.\cr \cr
#' \code{\%in\%} returns a logical ff vector indicating if there is a match or not for its left operand. 
#' ffdf objects are also allowed in the left and right operand of the \code{\%in\%} operator. See the examples.
#'
#' @export
#' @example ../examples/match.R
#' @param x a \code{ff} object for \code{ffmatch} or an \code{ffdf} object for \code{ffdfmatch}
#' @param table a \code{ff} object for \code{ffmatch} or an \code{ffdf} object for \code{ffdfmatch}
#' @param nomatch the value to be returned in the case when no match is found. Note that it is coerced to \code{integer}.
#' @param incomparables a vector of values that cannot be matched. Any value in \code{x} matching a value in this vector is assigned the nomatch value. For historical reasons, \code{FALSE} is equivalent to \code{NULL}.
#' @param trace logical indicating to show on which chunk the function is computing
#' @param ... other parameters passed on to chunk
#' @return An ff vector of the same length as \code{x}. An integer vector giving the position in table of the first match if there is a match, otherwise \code{nomatch}. 
#' @seealso \code{\link[base]{match}, \link[base]{paste}}
ffmatch <- function(x, table, nomatch = NA_integer_, incomparables = NULL, trace=FALSE, ...){
	stopifnot(inherits(x, "ff_vector") & inherits(table, "ff_vector"))
	
  nomatch <- as.integer(nomatch)
  xchunk <- chunk(x, ...)
  tablechunk <- chunk(table, ...)
  
  if(trace) {
    message(sprintf("%s, x has %s chunks, table has %s chunks", Sys.time(), length(xchunk), length(tablechunk)))
  }
  res <- ff(nomatch, length=length(x), vmode="integer") 
  ## First work on looping over x, then over the table
  for (i in xchunk){    
    Log$chunk(i)
   	xi <- x[i] 
    unmatched <- TRUE
    
    if (trace){
      message(sprintf("%s, working on x chunk %s:%s", Sys.time(), min(i), max(i)))
    }
    
    m <- rep(NA_integer_, sum(i))
    for (j in tablechunk){
      Log$chunk(j)
      if(trace) {
        message(sprintf("%s, working on table chunk %s:%s", Sys.time(), min(j), max(j)))      
      }
      m[unmatched] <- fmatch(x=xi[unmatched], table=table[j], incomparables=incomparables) +  min(j) - 1L      
      unmatched <- is.na(m)
      if (!any(unmatched)) break
    }
    m[unmatched] <- nomatch
    res[i] <- m
  }	
  res
}


#' @rdname ffmatch
#' @export
ffdfmatch <- function(x, table, nomatch = NA_integer_, incomparables = NULL, trace=FALSE, ...){
	stopifnot(inherits(x, "ffdf") & inherits(table, "ffdf"))
	
  nomatch <- as.integer(nomatch)
  xchunk <- chunk(x, ...)
  tablechunk <- chunk(table, ...)
  
  if(trace) {
    message(sprintf("%s, x has %s chunks, table has %s chunks", Sys.time(), length(xchunk), length(tablechunk)))
  }
  res <- ff(nomatch, length=nrow(x), vmode="integer") 
  ## First work on looping over x, then over the table
  for (i in xchunk){
    Log$chunk(i)
    xi <- x[i, , drop=FALSE]   
    unmatched <- TRUE
    
    if (trace){
      message(sprintf("%s, working on x chunk %s:%s", Sys.time(), min(i), max(i)))
    }
    
    m <- rep(NA_integer_, sum(i))
    for (j in tablechunk){
      if(trace) {
        message(sprintf("%s, working on table chunk %s:%s", Sys.time(), min(j), max(j)))      
      }
      m[unmatched] <- fmatch(x=do.call(paste, xi[unmatched, , drop=FALSE]), table=do.call(paste, table[j, , drop=FALSE]), incomparables=incomparables) +  min(j) - 1L     
      unmatched <- is.na(m)
      if (!any(unmatched)) break
    }
    m[unmatched] <- nomatch
    res[i] <- m
  }	
  res
}


in.default <- get(x="%in%")
#' @rdname ffmatch
#' @export
#' @usage x \%in\% table
"%in%" <- function(x, table){
	if(inherits(x, "ff_vector")){
		ffmatch(x=x, table=as.ff(table), nomatch = 0L) > 0L
	} else if(inherits(x, "ffdf") && inherits(table, "ffdf")){
		ffdfmatch(x=x, table=as.ffdf(table), nomatch = 0L) > 0L
	} else{
		in.default(x=x, table=table)
	}	
}


# quick testing
# ffmatch2(ff(factor(c("a", "c"))), ff(factor(c("b", "a"))), trace=TRUE, by=1)
