#' Create an index from a filter statement
#'
#' \code{ffwhich} creates an \code{\link{ff}} integer index vector
#' from a filter expression. The resulting vector can be used to index or subset
#' a ffdf or ff vector.
#' @example ../examples/ffwhich.R
#' @seealso ffindexget ffindexset
#' @param x \code{ff} or \code{ffdf} object
#' @param expr R code that evaluates to a logical
#' @param ... not used
#' @export
ffwhich <- function(x, expr, ...){
  UseMethod("ffwhich")
}

#' @method ffwhich ff_vector
#' @export
ffwhich.ff_vector <- function(x, expr, ..., envir=parent.frame()){
  #chunkify expression
  es <- deparse(substitute(expr))
  xs <- deparse(substitute(x))
  
  varre <- paste("\\b(",xs,")\\b", sep="")
  es <- gsub(varre, ".x[.i]", es)
  e <- parse(text=es)
  ###
  
  nl <- list(.x=x)
  
  fltr <- NULL
  for (.i in chunk(x, ...)){
    Log$chunk(.i)
    nl$.i = .i
    idx  <- which(eval(e, nl, envir)) +  min(.i) - 1L
    fltr <- ffappend(fltr, idx, ...)
  }
  fltr
}

#' @method ffwhich ffdf
#' @export
ffwhich.ffdf <- function(x, expr, ..., envir=parent.frame()){
  nl <- list(._x = x)
  es <- substitute(expr)  
  try( { if (is.expression(expr)){
           es <- expr
         } else {
           nl$.filter <- as.ff(expr)
           e <- expression(.filter[.i])
         }
       }
      , silent=TRUE
      )

  if (is.null(nl$.filter)){
    #### chunkify expression
    e <- chunkexpr(es, names(x), prefix="._x$")
  }
  ####
  
  #print(list(e=e, es=es))
  fltr <- NULL
  for (.i in chunk(x, ...)){
    Log$chunk(.i)
    nl$.i <- .i
    a <- which(eval(e, nl, envir)) +  min(.i) - 1L
    if (length(a))
      fltr <- ffappend(fltr, a)
  }
  fltr
}

###### quick testing
# x <- ff(10:1)
# idx <- ffwhich(x, x < 5)
# x[idx][]
# 
# dat <- ffdf(x1=x, y1=x)
# idx <- ffwhich(dat, x1 < 5 & y1 > 2)
# dat[idx,][,]
# f <- dat$x1[] > 3
# ffwhich(dat, f)
# ffwhich(dat, x > 2)
# ffwhich(dat, expression(x1 > 2))