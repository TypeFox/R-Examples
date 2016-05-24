##' Convert character or logical indices to numeric
##' @param x the object that is to be indexed
##' @param n names of the object
##' @param i the indices to be converted
##' @return numeric indices
##' @author Claudia Beleites
##' @export 
numericindex <- function (x, i, n = names (x)){
  if (is.character (i))
    match (i, n)
  else if (is.logical (i))
    seq_along (x) [i]
  else if (is.numeric (i))
    i
   else
    stop ("i must be numeric, logical, or character")
}
.test (numericindex) <- function (){
  checkEquals (numericindex (v, c("b", "a", "x")), c (2L, 1L, NA))
  checkEquals (numericindex (v, c(TRUE, FALSE, TRUE)), c(1L, 3L))
  checkEquals (numericindex (v, TRUE), 1 : 3)
  checkEquals (numericindex (v, FALSE), integer (0L))
  checkEquals (numericindex (v, c(TRUE, FALSE)), c(1L, 3L))
  checkEquals (numericindex (v, 1L), 1L)
}
