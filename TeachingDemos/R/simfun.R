simfun <- function(expr, drop, ...) {
  dots <- list(...)
  expr <- substitute(expr)
  has.drop <- !missing(drop)
  char2seed <- TeachingDemos::char2seed
  force(char2seed)
  function(data,seed) {
    if(!missing(seed)) {
      if(is.character(seed)) {
        char2seed(seed)
      } else {
        set.seed(seed)
      }
    }
    data.is.df <- FALSE
    if(!missing(data) && is.data.frame(data)) {
      data.is.df <- TRUE
      df.rn <- row.names(data)
      dots <- c(as.list(data),dots)
    } else if(!missing(data)) {
      dots <- c(as.list(data),dots)
    }
    outlist <- within(dots,eval(expr))
    if(has.drop) outlist[drop] <- NULL
    out.df <- as.data.frame(outlist)
    if(data.is.df) {
      row.names(out.df) <- df.rn
    }
    out.df
  }
}



