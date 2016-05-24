#' Duplicated for ff and ffdf objects
#'
#' Duplicated for ff and ffdf objects similar as in \code{\link[base]{duplicated}}.\cr
#' Remark that this duplicated function is slightly different from the duplicated method in the base package as it first orders
#' the ffdf or ff_vector object and then applies duplicated. This means you need to order the ffdf or ff_vector 
#' in case you want to have the exact same result as the result of the base package. See the example.
#'
#' @rdname duplicated.ff
#' @export
#' @export duplicated.ff
#' @method duplicated ff
#' @example ../examples/duplicated.R
#' @param x \code{ff} object or \code{ffdf} object
#' @param incomparables a vector of values that cannot be compared. 
#'   FALSE is a special value, meaning that all values can be compared, 
#'   and may be the only value accepted for methods other than the default. 
#'   It will be coerced internally to the same type as x.
#' @param fromLast logical indicating if duplication should 
#'   be considered from the last, i.e., the last (or rightmost) of identical elements will be kept
#' @param trace logical indicating to show on which chunk the function is computing
#' @param ... other parameters passed on to chunk
#' @return A logical ff vector of length \code{nrow(x)} or \code{length(x)} indicating if each row or element is duplicated.
#' @seealso \code{\link[base]{duplicated}, \link[ff]{ffdforder}, \link[ff]{fforder}}
duplicated.ff <- function(x, incomparables = FALSE, fromLast=FALSE, trace=FALSE, ...){
  if (!identical(incomparables, FALSE)){
    .NotYetUsed("incomparables != FALSE")
  }     
  
  res <- ff(vmode="logical", length= length(x))
  
  o <- fforder(x, decreasing = fromLast, na.last = TRUE)  
  
  i.last <- NULL
  for (i in chunk(x, ...)){
    if (trace){
      message(sprintf("%s, working on x chunk %s:%s", Sys.time(), min(i), max(i)))
    }
    Log$chunk(i)
    i.o <- o[i]
    i.x <- x[i.o]
    res[i.o] <- duplicated(i.x)
    if(min(i) > 1){
      res[i.o[1]] <- duplicated(c(i.last, utils::head(i.x, 1)))[2]
    }
    i.last <- tail(i.x, 1)
  }
  res
}


#' @rdname duplicated.ff
#' @method duplicated ffdf
#' @export
#' @export duplicated.ffdf
duplicated.ffdf <- function(x, incomparables = FALSE, fromLast=FALSE, trace=FALSE, ...){
  if (!identical(incomparables, FALSE)){
    .NotYetUsed("incomparables != FALSE")
  }     
  
  res <- ff(vmode="logical", length= nrow(x))
  
  o <- ffdforder(x, decreasing = fromLast, na.last = TRUE)  
  
  i.last <- NULL
  for (i in chunk(x, ...)){
    Log$chunk(i)
    if (trace){
      message(sprintf("%s, working on x chunk %s:%s", Sys.time(), min(i), max(i)))
    }
    i.o <- o[i]
    i.x <- x[i.o, , drop=FALSE]
    res[i.o] <- duplicated(i.x)
    if(min(i) > 1){
      res[i.o[1]] <- duplicated(rbind(i.last, utils::head(i.x, 1)))[2]
    }
    i.last <- tail(i.x, 1)
  }
  res
}
