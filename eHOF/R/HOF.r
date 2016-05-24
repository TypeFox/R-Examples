"HOF" <- function(...) UseMethod("HOF")


"[.HOF.list" <- function(x, s,...) {
    out <- NextMethod("[", drop=TRUE)
    class(out) <- if(length(s)>1) 'HOF.list' else 'HOF'
    return(out)
  }
