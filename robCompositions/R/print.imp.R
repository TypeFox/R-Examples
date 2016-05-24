#' Print method for objects of class imp
#' 
#' The function returns a few information about how many missing values are
#' imputed and possible other information about the amount of iterations, for
#' example.
#' 
#' 
#' @param x an object of class \sQuote{imp}
#' @param \dots additional arguments passed trough
#' @return None (invisible NULL).
#' @author Matthias Templ
#' @seealso \code{\link{impCoda}}, \code{\link{impKNNa}}
#' @keywords print
#' @export
#' @examples
#' 
#' data(expenditures)
#' expenditures[1,3]
#' expenditures[1,3] <- NA
#' \dontrun{
#' xi <- impCoda(expenditures)
#' xi
#' summary(xi)
#' plot(xi, which=1:2)
#' }
#' 
print.imp <- function(x, ...){
  cat("\n --------------------------------------- \n")
  if( x$w > 1 ){
    print(paste(x$w, "missing values were imputed"))
  } else{
    print(paste(x$w, "missing value was imputed"))
  }

  if( length(x$criteria) > 0 ){
    if( x$iter > 1 ){
      print(paste(x$iter, "iterations were needed"))
    } else{
      print(paste(x$iter, "iteration was needed"))
    }
    print(paste("the last change was", round(x$criteria,4)))
  }
  cat(" --------------------------------------- \n")
}
