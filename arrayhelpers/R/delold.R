##' Strip the attributes keeping track of the former shape
##'
##' Convenient for printing
##' @param a the array
##' @return \code{a} stripped of the \code{old.*} attributes.
##' @author Claudia Beleites
##' @export
##' @examples
##'
##' a <- arrayhelpers:::a
##' makeNd (a, 2)
##' delold (makeNd (a, 2))
##' 
##' 
delold <- function (a){
   attr (a, "old") <- NULL
   a
}
.test (delold) <- function (){
  tmp <- makeNd (a, 0)
  old <- attributes (tmp)
  old <- old [grepl ("^old", names (old))]
  checkIdentical (old,
                  list (old = list (list (
                          names = NULL, 
                          dimnames = list(rows = c("a", "b", "c", "d"),
                            columns = c("A", "B", "C"), 
                            d3 = c("1", "2")),
                          dim = c(4L, 3L, 2L))))
                  )
  
  tmp <- delold (tmp)
  checkTrue (is.null (attr (old, "old")))
}
