#' Define a collection of multivariate polynomials.
#' 
#' Combine a series of mpoly objects into a mpolyList.
#' 
#' @param ... a series of mpoly objects.
#' @return An object of class mpolyList.
#' @export
#' @examples
#' ( p1 <- mp("t^4 - x") )
#' ( p2 <- mp("t^3 - y") )
#' ( p3 <- mp("t^2 - z") )
#' ( ms <- mpolyList(p1, p2, p3) )
#' is.mpolyList( ms )
#' 
#' mpolyList(mp("x + 1"))
#' p <- mp("x + 1")
#' mpolyList(p)
#' 
#' ps <- mp(c("x + 1", "y + 2"))
#' is.mpolyList(ps)
#' 
#' 
#' f <- function(){
#'   a <- mp("1")
#'   mpolyList(a)
#' }
#' f()
#' 
#' 
mpolyList <- function(...){
	
  arguments <- as.list(match.call()[-1])  

  out <- lapply(arguments, eval, parent.frame(1))

  if(is.mpoly(out)) out <- list(out)
  
  if(!all(vapply(out, is.mpoly, logical(1)))){
  	stop("each argument must be of class mpoly.", call. = FALSE)
  }
  
  class(out) <- "mpolyList"
  out
}
