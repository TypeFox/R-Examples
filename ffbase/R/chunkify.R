#' Chunkify an element-wise function
#' 
#' Chunkify creates a new function that operates on a ff vector. 
#' It creates chunks from the ff vector and calls the orginal function \code{fun} on each chunk.
#' @export
#' @param fun function to be 'chunkified', the function must accept a vector and 
#'    return a vector of the same \code{length}
#' @return 'chunkified' function that accepts a \code{ff} vector as its first argument.
chunkify <- function (fun){
  cfun <- function(x, ..., inplace = FALSE) {
    chunks <- chunk(x, ...)
    i <- chunks[[1]]
    res <- fun(x[i], ...)
    if(inherits(res, "character")){
      res <- as.ff(as.factor(res))
      ret <- ffappend(x=NULL, y=res)
      for (i in chunks[-1]) {
        Log$chunk(i)
        res <- as.ff(as.factor(fun(x[i], ...)))
        ret <- ffappend(x=ret, y=res)
      }
      ret
    }else{
      ret <- as.ff(res)
      length(ret) <- length(x)
      for (i in chunks[-1]) {
        Log$chunk(i)
        ret[i] <- fun(x[i], ...)
      }
      ret
    }
    
  }
  cfun
}

#' Chunk an expression 
#' 
#' chunkexpr replaces variables in an expression with a indexed version.
#' It main use it to rewrite "normal" R expression into chunked versions
#' that can be evaluated in a chunked-for-loop.
#'@param expr \code{expression} vector or language object
#'@param x \code{character} with variables to be chunked
#'@param i name of index that will be used in the for loop, typically a \code{ri} or \code{hi}.
#'@param prefix prefix for variables to be replaced.
#'@keywords internal
chunkexpr <- function(expr, x = all.vars(expr), i=".i", prefix=""){
  es <- lapply(as.expression(expr), deparse)
  es <- lapply(es, paste0, collapse="\n")
  xs <- x
  for (var in xs){
    varre <- paste("\\b(",var,")\\b", sep="")
    varsub <- paste(prefix, "\\1[",i,"]", sep="")
    es <- gsub(varre, varsub, es)
  }
  parse(text=es)
}

#chunkexpr(c("x","y"), c("x>2 & y==1\nz==3", "y > 3"), i=".i", prefix="data$")
