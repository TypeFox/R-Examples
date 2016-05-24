#' Select Non Factor Variables
#'
#' This is a function to select non-factor variables for the \code{Rcmdr} package.
#'
#' @seealso \code{\link[Rcmdr:Rcmdr.Utilities]{Rcmdr.Utilities}}
#'
#' @rdname util-nonFactors
#' @keywords hplot
#' @export
nonFactors <- function() {

  setdiff(Variables(), Factors())

}



#' Check Non Factor Variables in The Active Dataset
#'
#' This is a function to check the active dataset has more than or equal to \code{n} non-factor variables.
#'
#' @param n Numeric: the target number of \code{n}.
#' @return
#' \describe{
#' \item{result}{Boolean; the active dataset has or does not have more than or equal to \code{n} non-factor variables.}
#' }
#' @seealso \code{\link[Rcmdr:Rcmdr.Utilities]{Rcmdr.Utilities}}
#'
#' @rdname util-nonFactorsP
#' @keywords hplot
#' @export
nonFactorsP <- function (n = 1) {

  activeDataSetP() && length(setdiff(listVariables(), listFactors())) >= n

}
