#' Debugging function
#' 
#' A simple but useful debugging function. It first announces the object to
#' printed and then prints it.
#' 
#' 
#' @param x The object to be printed
#' @param method The method for printing \code{x}. Default is \code{"print"},
#' which uses \code{\link{print}} for printing; \code{"cat"} uses
#' \code{\link{cat}} for printing. The latter is useful for short objects
#' (scalar and vectors), the former for more structured objects (data frames,
#' matrices, lists etc).
#' @author Hein Putter \email{H.Putter@@lumc.nl}
#' @keywords univar
#' @examples
#' 
#' tm <- c(0.2,0.5,1,1.2,1.8,4)
#' ta <- 2*tm
#' dfr <- data.frame(time=tm, stepf=ta)
#' deb(dfr, method="print")
#' deb(tm, method="cat")
#' 
#' @export deb
deb <- function(x, method=c("print","cat")) {
  call <- match.call()
  method <- match.arg(method)
  if (method=="print") {
    cat(deparse(call$x), ":\n")
    print(x)
  }
  if (method=="cat") cat(deparse(call$x), ":", x, "\n")
  flush.console()
}
