#' Summary method for objects of class imp
#' 
#' A short comparison of the original data and the imputed data is given.
#' 
#' Note that this function will be enhanced with more sophisticated methods in
#' future versions of the package.  It is very rudimental in its present form.
#' 
#' @param object an object of class \sQuote{imp}
#' @param \dots additional arguments passed trough
#' @return None (invisible NULL).
#' @author Matthias Templ
#' @seealso \code{\link{impCoda}}, \code{\link{impKNNa}}
#' @keywords print
#' @export
#' @method summary imp
#' @examples
#' 
#' data(expenditures)
#' expenditures[1,3]
#' expenditures[1,3] <- NA
#' xi <- impKNNa(expenditures)
#' xi
#' summary(xi)
#' # plot(xi, which=1:2)
#' 
summary.imp <- function(object, ...){
  geometricmean <- function (x) {
    if (any(na.omit(x == 0)))
        0
    else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
  }
gm <- apply(object$xOrig, 2, function(x) {
  geometricmean(as.numeric(x[complete.cases(x)]))
})
gmI <- apply(object$xImp, 2, function(x) {
  geometricmean(as.numeric(x[complete.cases(x)])) ## gewichten!
})

  d <- data.frame(orig=gm,
                  imp=gmI)
  cat("\n geometric mean of the original data and the imputed data: \n")

  ## Einfluss der Imputation mittels bootstrap schaetzen
  ## laut script sim

  d
}
