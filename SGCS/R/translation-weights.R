#' Translation weights
#' Compute the translation edge correction weights
#' 
#' @param x Some form of point pattern data.
#' @useDynLib SGCS
#' @export 

translation_weights <- function(x) {
  if(0){
    x <- internalise_pp(x)
    return(.External("translation_weights_c", x, PACKAGE="SGCS"))
  }
  if(class(x)=="ppp"){
    w <- 1/(edge.Trans(x)/area(x$window))
    w <- w[lower.tri(w)]
  }
  else{
    x <- internalise_pp(x)
    if(x$dim == 2){
      w <- translation_weights(internal_to_ppp(x))
    }
    else{
      w <- .External("translation_weights_c", x, PACKAGE="SGCS")
    }
  }
  w
}