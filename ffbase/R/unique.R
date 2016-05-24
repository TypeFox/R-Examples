#' Unique values for ff and ffdf objects
#'
#' @rdname unique.ff
#' @export
#' @export unique.ff
#' @method unique ff
#' @example ../examples/unique.R
#' @param x \code{ff} object or \code{ffdf} object
#' @param incomparables a vector of values that cannot be compared. 
#'   FALSE is a special value, meaning that all values can be compared, 
#'   and may be the only value accepted for methods other than the default. 
#'   It will be coerced internally to the same type as x.
#' @param fromLast logical indicating if duplication should 
#'   be considered from the last, i.e., the last (or rightmost) of identical elements will be kept
#' @param trace logical indicating to show on which chunk the function is computing
#' @param ... other parameters passed on to chunk
#' @return An ffdf with unique values in \code{x} or an ff vector with unique values in \code{x} if x is a ff vector. 
#' @seealso \code{\link[base]{unique}}
unique.ff <- function(x, incomparables = FALSE, fromLast = FALSE, trace=FALSE, ...){
  if (!identical(incomparables, FALSE)){
    .NotYetUsed("incomparables != FALSE")
  }
  if(vmode(x) == "integer" & length(res <- levels(x))>0){    
    ## Something strange is happening for factors with fforder, reported to ff maintainer, doing a workaround
    if(any(is.na(x))){
      res <- c(res, NA)
    }
    res <- ff(res, levels = res)
  }else{
    ## Order the ff    
    xorder <- fforder(x, decreasing = fromLast, na.last = TRUE)
    xchunk <- chunk(x, ...)
    ## Chunkwise adding of unique elements to the unique ff_vector called res
    res <- NULL
    lastel <- NULL
    for (i in xchunk){
      #if (trace){
      #  message(sprintf("%s, working on x chunk %s:%s", Sys.time(), min(i), max(i)))
      #}
      Log$chunk(i)
      iorder <- xorder[i]
      iorder <- as.integer(iorder) # make sure it is not a Date
      xi <- x[iorder]
      xi <- unique(xi)
      ## exclude the first row if it was already in the unique ffdf as this is the last one from the previous unique
      if(sum(duplicated(c(xi[1], lastel)))>0){
        xi <- xi[-1]  
      }   
      if(length(xi) > 0){   
        ## Add the result to an ff_vector
        lastel <- xi[length(xi)]
        res <- ffappend(x=res, xi)
      }
    }
  }  
  res
}

#' @rdname unique.ff
#' @method unique ffdf
#' @export
unique.ffdf <- function(x, incomparables = FALSE, fromLast=FALSE, trace=FALSE, ...){
  if (!identical(incomparables, FALSE)){
    .NotYetUsed("incomparables != FALSE")
  }     
  ## Order the ffdf    
  xorder <- ffdforder(x, decreasing = fromLast, na.last = TRUE)
  xchunk <- chunk(x, ...)
  ## Chunkwise adding of unique rows to the unique ffdf called res
  res <- NULL
  lastrow <- NULL
  for (i in xchunk){
    if (trace){
      message(sprintf("%s, working on x chunk %s:%s", Sys.time(), min(i), max(i)))
    }
    iorder <- xorder[i]
    iorder <- as.integer(iorder) # make sure it is not a Date
    xi <- x[iorder, ]
    xi <- unique(xi)
    rownames(xi) <- NULL
    ## exclude the first row if it was already in the unique ffdf as this is the last one from the previous unique
    if(sum(duplicated(rbind(xi[1, , drop=FALSE], lastrow)))>0){
      xi <- xi[-1, , drop=FALSE]  
    }   
    if(nrow(xi) > 0){   
      ## Add the result to an ffdf
      lastrow <- xi[nrow(xi), , drop=FALSE]
      if(!is.null(res)){
        rownames(xi) <- (nrow(res)+1):(nrow(res)+nrow(xi))
      }
      res <- ffdfappend(x=res, dat=xi, recode=FALSE)        
    }
  }
  res
}

