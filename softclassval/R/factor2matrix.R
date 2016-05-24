##' Convert hard class labels to membership matrix
##'
##' Converts a factor with hard class memberships into a membership matrix
##' @param f factor with class labels
##' @return matrix of size \code{length (f)} x \code{nlevels (f)}
##' @author Claudia Beleites
##' @seealso \code{\link[softclassval]{hardclasses}} for the inverse 
##' @export
factor2matrix <- function (f){
	if (! is.factor (f))
		f <- as.factor (f)
	
	res <- matrix (0, nrow = length (f), ncol = nlevels (f))
	colnames (res) <- levels (f)
	
	res [cbind (seq_along (f), as.numeric (f))] <- 1
	
	res
}

.test (factor2matrix) <- function (){
  checkEquals (factor2matrix (c ("a", "b", "a", NA)),
               structure(c(1, 0, 1, 0, 0, 1, 0, 0), .Dim = c(4L, 2L),
                         .Dimnames = list(NULL, c("a", "b"))))
}
