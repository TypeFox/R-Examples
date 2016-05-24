test <- function # simple function arguments
## this function does nothing in particular and does it very well
( first,         ##<< the first argument
  second=TRUE    ##<< the second argument
 ){
  ##<<alias Iolanthe
  ##<<details
  ## if second is TRUE then first is returned
  if ( second ){
    res <- first
  } else {
    ##<<details
    ## if second is not TRUE then first is negated
    res <- !first
  }
  invisible(res)
### invisible something not unrelated to first
}
