#' Distances to observation window edge
#' 
#' Compute the distances from points to observation window edge.
#' 
#' @param x Point pattern
#' 
#' Uses bdist.points for ppp objects, otherwise relies on a cuboidal/rectangular window.
#' 
#' @export

edge_distance <- function(x) {
  if(class(x)=="ppp"){
    w <- bdist.points(x)
  }
  else{
    x <- internalise_pp(x)
    if(x$dim == 2){
      w <- edge_distance(internal_to_ppp(x))
    }
    else{ # 3 dimension, not spatstat format
      v<-sapply(1:3, function(i) pmin(x$bbox[2,i]-x$coord[,i], x$coord[,i]-x$bbox[1,i] ))
      w <- apply(v, 1, min)
    }
  }
  w
}