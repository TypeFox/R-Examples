listNbh <- function(object,y,val){
  stopifnot(class(object)=="neighborhood")
  if (missing(y)) y=object$neighbors
  else{
    stopifnot(length(y)==object$n)
    y=y[object$neighbors]
  }
  if (missing(val)) val <- object$values
  posVal <- match(val,object$values,nomatch=FALSE)
  stopifnot(all(posVal!=0))
  out <- split(y,rep(1:object$nu,object$size.nbh))[posVal]
  out
}
