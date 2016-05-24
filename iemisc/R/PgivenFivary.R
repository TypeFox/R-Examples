#' "Present equivalent of a series of future cash flows subject to varying interest rates" (Engineering Economics)
#'
#' Compute P given F and i that varies
#'
#' P is expressed as
#'
#' 	\deqn{P = \frac{F_n}{\prod \limits_{k=1}^n{\left(1 + i_k\right)}}}
#'
#' \describe{
#'	\item{\emph{P}}{the "present equivalent"}
#'	\item{\emph{\eqn{F_n}}}{the "future cash flows subject to varying interest rates"}
#'	\item{\emph{\eqn{i_k}}}{the "interest rate for the kth period"}
#'	\item{\emph{k}}{the "number of interest periods"}
#' }
#'
#' @param Fn numeric vector that contains the future value(s) at the end of a
#'    period n
#' @param k numeric vector that contains the kth period values
#' @param ik numeric vector that contains the effective interest rate(s) per
#'    period as a percent for the kth period
#'
#' @return PgivenFivary numeric vector that contains the present value(s)
#'
#'
#' @source
#' \enumerate{
#'    \item r - Add a Column to a Dataframe From a List of Values - Stack Overflow answered by Matthew Plourde on Jun 21 2012. See \url{http://stackoverflow.com/questions/11130037/add-a-column-to-a-dataframe-from-a-list-of-values/11130178}.
#'    \item r - Why does is.vector() return TRUE for list? - Stack Overflow answered by Andrie on May 17 2011. See \url{http://stackoverflow.com/questions/6032772/why-does-is-vector-return-true-for-list/6032909}.
#' }
#'
#' @references
#' William G. Sullivan, Elin M. Wicks, and C. Patrick Koelling, \emph{Engineering Economy}, Fourteenth Edition, Upper Saddle River, New Jersey: Pearson/Prentice Hall, 2009, page 142, 162.
#'
#' @encoding UTF-8
#'
#'
#' @examples
#' library(iemisc)
#' # Example for equation 4-31 from the Reference text (page 162)
#' PgivenFivary(Fn = 1000, ik = c(10, 12, 13, 10), k = 1)
#' # i1 is 10%, i2 is 12%, i3 is 14%, and i4 is 10% & k = 1 year
#'
#'
#'
#'
#' @export
PgivenFivary <- function (Fn, ik, k) {

ik <- ik / 100

ika <- vector("list", length(ik)) # Source 1 and 2 / pre-allocate the list since it is being used in a for loop

for (t in 1:length(ik)) {

ika[[t]] <- (1 / (1 + ik[t])) ^ k
}

PgivenFivary <- Fn * prod(unlist(ika))

return(round(PgivenFivary, digits = 2))
}
