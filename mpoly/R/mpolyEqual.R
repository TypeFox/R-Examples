#' Determine whether two multivariate polynomials are equal.
#' 
#' Determine whether two multivariate polynomials are equal.
#' 
#' @param e1 an object of class mpoly
#' @param e2 an object of class mpoly
#' @method == mpoly
#' @aliases ==.mpoly ==
#' @return A logical value.
#' @name mpolyEqual
#' @export
#' @examples
#' 
#' p1 <- mp("x + y + 2 z")
#' p1 == p1
#' 
#' p2 <- reorder(p1, order = "lex", varorder = c("z","y","x"))
#' p1 == p2
#' p2 <- reorder(p1, order = "lex", varorder = c("z","w","y","x"))
#' p1 == p2
#' p1 == ( 2 * p2 )
#' 
#' p1 <- mp("x + 1") 
#' p2 <- mp("x + 1")
#' identical(p1, p2)
#' p1 == p2
#' 
#' mp("x + 1") == mp("y + 1")
#' mp("2") == mp("1")
#' mp("1") == mp("1")
#' mp("0") == mp("-0")
#' 
`==.mpoly` <- function(e1, e2){
	
  if(!is.mpoly(e1)  || !is.mpoly(e2)){
    stop("e1 and e2 must be of class mpoly.", call. = FALSE)
  }
  
  diff <- e1 - e2
  
  if(
    length(diff) == 1 && 
    length(diff[[1]]) == 1 && 
    diff[[1]] == 0
  ){
    return(TRUE)
  } else {
    return(FALSE)	
  }
  
}




