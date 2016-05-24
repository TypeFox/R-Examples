#' print a multcompLetters object
#' 
#' print method for an object of class 'multcompLetters'.
#' 
#' Prints only the Letters component of the 'multcompLetters' list unless
#' all=TRUE.
#' 
#' @param x an object of class 'multcompLetters'
#' @param all FALSE to print only the character vector representations of the
#' 'multcompLetters' comparison summary; TRUE to print also the matrix
#' representation.
#' @param ...  Other optional print parameters as described on the
#' \code{\link{print}} help page.
#' @return x\$Letters = the named, character vector representation of the
#' 'multcompLetters' evaluation of the distance relationships
#' @author Spencer Graves
#' @seealso \code{\link{multcompLetters}}
#' @keywords dplot
#' @export
#' @examples
#' 
#' dif3 <- c(FALSE, FALSE, TRUE)
#' names(dif3) <- c("A-B", "A-C", "B-C")
#' dif3L <- multcompLetters(dif3)
#' dif3L
#' print(dif3L)
#' print(dif3L, TRUE)
#' 
"print.multcompLetters" <-
function(x, all=FALSE, ...){
  {
    if(all){
      class(x) <- NULL
      print(x, ...)
    }
    else 
      print(x$Letters, ...)
  }
  invisible(x$Letters)
}

