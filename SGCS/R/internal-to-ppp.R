#' Data to spatstat format
#' 
#' @param x Pattern in internal format
#' 
#' @export
internal_to_ppp <-function(x) {
  x <- internalise_pp(x) # just to be sure
  win <- if(!is.null(x$owin)) x$owin else as.owin(c(x$bbox))
  coords <- x$coord
  if(x$dim==2)
    ppp(coords[,1], coords[,2], window=win, mark=cbind(mass=x$mass, type=x$type))
  else 
    pp3(coords[,1], coords[,2], coords[,3], c(x$bbox))
}