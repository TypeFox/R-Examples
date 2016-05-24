#' Test whether an mpoly object is linear.
#'
#' Test whether an mpoly object is linear.
#'
#' @param x an mpoly or mpolyList object
#' @return a logical vector
#' @examples
#' \dontrun{
#' 
#' is.linear(mp("0"))
#' is.linear(mp("x + 1"))
#' is.linear(mp("x + y"))
#' is.linear(mp(c("0", "x + y")))
#' 
#' is.linear(mp("x + x y"))
#' is.linear(mp(c("x + x y", "x")))
#' 
#' 
#' }
is.linear <- function(x){
  
  stopifnot(is.mpoly(x) || is.mpolyList(x))
  
  if(is.mpolyList(x)) return(sapply(x, is.linear))

  all(
    sapply(x, function(term){
      if(all(length(term) <= 2)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    })
  )  
}

