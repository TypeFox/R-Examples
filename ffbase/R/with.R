#' Evaluate an expression in a ffdf data environment 
#' 
#' Evaluate an R expression in an environment constructed from a ffdf data frame.
#' (see \code{\link{with}}). Please note that you should write
#' your expression as if it is a normal \code{data.frame}. The resulting return value
#' however will be a \code{ff} object.
#' 
#' @note `with.ffdf` assumes that the returned object is of equal length as 
#' `nrow(data)` and must be converted to a `ff` object
#' In case this is not true, the result won't be correct.
#' 
#' @seealso \code{\link{ffdfwith}}
#' @method with ffdf 
#' @export
#' @export with.ffdf
#'
#' @example ../examples/with.R
#' @param data \code{\link{ffdf}} data object used as an environment for evaluation.
#' @param expr expression to evaluate.
#' @param ... arguments to be passed to \code{\link{chunk}}.
#' @return if expression is a \code{vector} a newly created \code{ff} vector will be returned 
#' otherwise if the expression is a data.frame a newly created \code{ffdf} object will be returned.
with.ffdf <- function(data, expr, ...){
   e <- substitute(expr)
   #chunks <- chunk(data, by=2) #debug chunking
   chunks <- chunk(data, ...)
   
   cdat <- data[chunks[[1]],,drop=FALSE]
   res <- eval(e, cdat, enclos=parent.frame())
   if (NROW(res)!= nrow(cdat)){
     stop("'with.ffdf' only returns `ff` object of equal length of `nrow(data)`")          
   }
   fc <- FALSE
   
#    if (!is.atomic(res) && !is.data.frame(res)){
#      stop("'with.ffdf' only returns `ff` object of equal length of `nrow(data)`")
#    }
   
   if (is.character(res) || is.factor(res)){
     res <- as.factor(res)
     fc <- TRUE
   } else if (is.data.frame(res)){
     fc <- sapply(res, function(x) is.factor(x) || is.character(x))
     res[fc] <- lapply(res[fc], as.factor)
   }
   if (is.vector(res) || is.factor(res) || inherits(res, "Date") || inherits(res, "POSIXct")){
      res <- as.ff(res)
      length(res) <- nrow(data)
      for (i in chunks[-1]){
        Log$chunk(i)
        d_i <- data[i,,drop=FALSE]
        r <- eval(e, d_i, enclos=parent.frame())
        
        if (length(r)!= nrow(d_i)){
          stop("'with.ffdf' only returns `ff` object of equal length of `nrow(data)`")          
        }
        
        if (fc){ 
             r <- as.factor(r)
             levels(res) <- appendLevels(res, levels(r))
         }
         res[i] <- r
      }
   } else if (is.data.frame(res)){
      res <- as.ffdf(res)
      rownames(res) <- NULL
      nrow(res) <- nrow(data)
      for (i in chunks[-1]){
        Log$chunk(i)
        d_i <- data[i,,drop=FALSE]
        r <- eval(e, d_i, enclos=parent.frame())
        
        if (nrow(r)!= nrow(d_i)){
          stop("'with.ffdf' only returns `ff` object of equal length of `nrow(data)`")          
        }
        if (any(fc)){
           r[fc] <- lapply(which(fc), function(x) {
                r[[x]] <- as.factor(r[[x]])
                levels(res[[x]]) <<- appendLevels(res[[x]], r[[x]])
                r[[x]]
             })
        }
        res[i,] <- r
      }
   } else {
     stop("'with.ffdf' only returns `ff` object of equal length of `nrow(data)`")          
   }
   res
}

#' Evaluate an expression in a ffdf data environment 
#' 
#' Same functionality as \code{\link{within}}. Please note that you should write
#' your expression as if it is a normal \code{data.frame}. The resulting data.frame
#' however will be a new \code{ffdf} data.frame.
#' @method within ffdf 
#' @export
#'
#' @example ../examples/within.R
#' @param data \code{\link{ffdf}} data object used as an environment for evaluation.
#' @param expr expression to evaluate.
#' @param ... arguments to be passed to \code{\link{chunk}}.
#' @return a modified clone of \code{data}.
within.ffdf <- function(data, expr, ...){
    expr <- substitute(expr)
    parent <- parent.frame()
    
    #chunks <- chunk(data, by=2) debug chunking
    chunks <- chunk(data, ...)
    cdat <- data[chunks[[1]],,drop=FALSE]
   
    e <- evalq(environment(), cdat, parent)
    eval(expr, e)
    l <- as.list(e)
    
    l <- l[!sapply(l, is.null)]
    del <- setdiff(names(cdat), names(l))
    #delete 
    cdat[del] <- list()
    cdat[names(l)] <- l

    res <- as.ffdf(cdat)
    rownames(res) <- NULL
    nrow(res) <- nrow(data)
    rownames(res) <- rownames(data)
    for (i in chunks[-1]){
       cdat <- data[i,,drop=FALSE]
       e <- evalq(environment(), cdat, parent)
       eval(expr, e)
       l <- as.list(e)
       l <- l[!sapply(l, is.null)]
       cdat[names(l)] <- l
       cdat[del] <- list()

       for(f in names(cdat)[sapply(cdat, is.factor)]) {
           levels(res[f]) <- appendLevels(levels(res[f]), levels(cdat[f]))
       }

       res[i,] <- cdat
    }
    res
}
