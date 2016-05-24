is.number <- function(char) str_detect(char, "[0-9]{1}")
is.period <- function(char) str_detect(char, "[\\.]{1}")
is.letter <- function(char) str_detect(char, "[a-zA-Z]{1}")







#' mpoly predicate functions
#' 
#' Various functions to deal with mpoly and mpolyList objects.
#' 
#' @param x object to be tested
#' @return Vector of logicals.
#' @name predicates
#' @examples
#' 
#' p <- mp("5")
#' is.mpoly(p)
#' is.constant(p)
#' 
#' is.constant(mp(c("x + 1", "7", "y - 2")))
#' 
#' p <- mp("x + y")
#' is.mpoly(p)
#' 
#' is.mpolyList(mp(c("x + 1", "y - 1")))
#'
#' 
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
#' (p <- bernstein(2, 5))
#' is.mpoly(p)
#' is.bernstein(p)
#' 
#' (p <- chebyshev(5))
#' is.mpoly(p)
#' is.chebyshev(p)
#' str(p)








#' @export
#' @rdname predicates
is.constant <- function(x){
  if(is.mpoly(x) && length(x) == 1 && length(x[[1]]) == 1) return(TRUE)
  if(is.mpolyList(x)) return(vapply(x, is.constant, logical(1)))
  FALSE
}






#' @export
#' @rdname predicates
is.mpoly <- function(x) any(class(x) == "mpoly")






#' @export
#' @rdname predicates
is.bernstein <- function(x) any(class(x) == "bernstein")




#' @export
#' @rdname predicates
is.bezier <- function(x) any(class(x) == "bezier")






#' @export
#' @rdname predicates
is.chebyshev <- function(x) any(class(x) == "chebyshev")









#' @export
#' @rdname predicates
is.mpolyList <- function(x){
  if(any(class(x) == "mpolyList")){
    return(TRUE)  
  } else {
    return(FALSE)	
  }
}









#' @export
#' @rdname predicates
is.linear <- function(x){
  
  stopifnot(is.mpoly(x) || is.mpolyList(x))
  
  if(is.mpolyList(x)) return(vapply(x, is.linear, logical(1)))
  
  all(
    vapply(x, function(term) all(length(term) <= 2), logical(1))
  )  
}














#' Test whether an object is a whole number
#'
#' Test whether an object is a whole number.
#'
#' @param x object to be tested
#' @param tol tolerance within which a number is said to be whole
#' @return Vector of logicals.
#' @export
#' @examples
#' is.wholenumber(seq(-3,3,.5))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

