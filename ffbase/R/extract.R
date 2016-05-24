#' Reading and writing vectors extended to handle logical \code{ff} vectors as indexes
#'
#' Package \code{ff} does not allow to extract and set values of \code{ff} vectors based on logical \code{ff} vectors. For this reason
#' the extractor functions \code{[.ff} and \code{[<-.ff} defined in package ff are overloaded.\cr
#' If you supply a logical \code{ff} vector as an index to another \code{ff} vector, the overloaded function will convert it to an integer \code{ff}. 
#' index before using the \code{[.ff} and \code{[<-.ff} function from the ff package. \cr
#' This allows to do \code{ff(1:10)[ff(c(FALSE, TRUE, NA, TRUE))]}\cr\cr
#' Mark that all other functionality from the extractor functions \code{[.ff} and \code{[<-.ff} in package ff are retained. This is an extension 
#' to handle logical \code{ff} vectors as indexes.  
#'
#' @export
#' @export [.ff
#' @rdname ffextract
#' @example ../examples/extract.R
#' @param x an \code{ff} object 
#' @param i missing OR a single index expression OR a \code{\link[ff]{hi}} object 
#' @param add TRUE if the values should rather increment than overwrite at the target positions, see \code{\link[ff]{readwrite.ff}}
#' @param pack FALSE to prevent rle-packing in hybrid index preprocessing, see \code{\link[ff]{as.hi}}
#' @param value the values to be assigned, possibly recycled
#' @usage \method{[}{ff} (x, i, pack = FALSE)
#' @return See \code{\link[ff]{Extract.ff}}. Mark that if a logical \code{ff} vector is used for \code{i}, and if only \code{FALSE} or \code{NA} 
#' values are present, NULL is returned in case of the extractor function \code{[.ff} while for the setter function \code{[<-.ff}, if the length value
#' is zero, this is not allowed.
#' @seealso \code{\link[ff]{Extract.ff}}
`[.ff` <- function(x, i, pack = FALSE){
  if(!missing(i) && is.ff(i) && length(i) > 0 && is.logical(i[1])){
    idx <- ffwhich(i, i==TRUE)    
    if(length(idx) == 0){
      warning("you are indexing ff vectors which return zero length, do you really need to use ff?")
      return(NULL)
    }    
    finalizer(idx) <- "delete"
    ff::`[.ff`(x=x, i=idx, pack=pack)
  }else{
    ff::`[.ff`(x=x, i=i, pack=pack)
  }	
}

#' @rdname ffextract
#' @usage \method{[}{ff} (x, i, add = FALSE, pack = FALSE) <- value
#' @export
#' @export [<-.ff
`[<-.ff` <- function(x, i, add = FALSE, pack = FALSE, value){
  if(!missing(i) && is.ff(i) && length(i) > 0 && is.logical(i[1])){
    idx <- ffwhich(i, i==TRUE)    
    if(length(idx) == 0){
      stop("no value for replacement")
    }    
    finalizer(idx) <- "delete"
    ff::`[<-.ff`(x=x, i=idx, add=add, pack=pack, value=value)
  }else{
    ff::`[<-.ff`(x=x, i=i, add=add, pack=pack, value=value)
  }	
}


#' Reading and writing data.frames (ffdf)
#'
#' Package \code{ff} does not allow to extract and set values of \code{ffdf} objects based on logical \code{ff} vectors. For this reason
#' the extractor functions \code{[.ffdf} and \code{[<-.ffdf} defined in package ff are overloaded.\cr
#' If you supply a logical \code{ff} vector as an index to subset an \code{ffdf} object, the overloaded function will convert the logical \code{ff} vector
#' to an integer \code{ff} index before using the \code{[.ffdf} and \code{[<-.ffdf} functions from the ff package. \cr
#' This allows to do \code{as.ffdf(iris)[as.ff(iris$Sepal.Length > 5), ]}\cr\cr
#' This is an extension to handle logical \code{ff} vectors as indexes to \code{ffdf} objects.
#'
#' @export
#' @export [.ffdf
#' @rdname ffextractffdf
#' @example ../examples/extractffdf.R
#' @param x an \code{ff} object 
#' @param i a row subscript 
#' @param j a column subscript
#' @param drop logical. If TRUE the result is coerced to the lowest possible dimension. 
#' @param value A suitable replacement value
#' @usage \method{[}{ffdf} (x, i, j, drop = TRUE)
#' @return See \code{\link[ff]{Extract.ffdf}}. Mark that if a logical \code{ff} vector is used for \code{i}, and if only \code{FALSE} or \code{NA} 
#' values are present, this is not allowed as ffdf with zero rows do not exist.
#' @seealso \code{\link[ff]{Extract.ffdf}}
"[.ffdf" <- function(x, i, j, drop = TRUE){  
  if(!missing(i) && is.ff(i) && length(i) > 0 && is.logical(i[1])){
    idx <- ffwhich(i, i == TRUE)
    if(length(idx) == 0){
      stop('Your expression selected zero rows, but ffdf objects cannot have zero rows')
    }
    finalizer(idx) <- "delete"
    if(missing(j)){
      ff::"[.ffdf"(x, i=idx, j=colnames(x), drop = drop)
    }else{
      ff::"[.ffdf"(x, i=idx, j=j, drop = drop)
    }
  }else{
    if(missing(i) & missing(j)){
      ff::"[.ffdf"(x, i=i, j=j, drop = drop) 
    }else if(missing(j)){
      n.args <- nargs() - !missing(drop)
      if(n.args == 2){
        ff::"[.ffdf"(x, i=i)              
      }else{
        ff::"[.ffdf"(x, i=i, j=colnames(x), drop=drop)
      }      
    }else if(missing(i)){
      ff::"[.ffdf"(x, i=i, j=j, drop = drop)  
    }else{
      ff::"[.ffdf"(x, i=i, j=j, drop = drop)  
    }    
  }
}


#' @rdname ffextractffdf
#' @usage \method{[}{ffdf} (x, i, j) <- value
#' @export
#' @export [<-.ffdf
"[<-.ffdf" <- function(x, i, j, value){  
  if(!missing(i) && is.ff(i) && length(i) > 0 && is.logical(i[1])){
    idx <- ffwhich(i, i == TRUE)
    if(length(idx) == 0){
      stop('Your expression selected zero rows, but ffdf objects cannot have zero rows')
    }
    finalizer(idx) <- "delete"
    if(missing(j)){
      ff::"[<-.ffdf"(x, i=idx, j=colnames(x), value=value)
    }else{
      ff::"[<-.ffdf"(x, i=idx, j=j, value=value)
    }
  }else{
    if(missing(i) & missing(j)){
      ff::"[<-.ffdf"(x, value=value) 
    }else if(missing(j)){
      ff::"[<-.ffdf"(x, i=i, j = colnames(x), value=value)        
    }else if(missing(i)){
      ff::"[<-.ffdf"(x, i=i, j=j, value=value)  
    }else{
      ff::"[<-.ffdf"(x, i=i, j=j, value=value)
    }
  }
}