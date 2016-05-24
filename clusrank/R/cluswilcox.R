#' Wilcoxon Rank SUm and Signed Rank Tests for Clustered Data
#'
#' Performs one- and two-sample Wilcoxon tests on vectors of clustered data.
#'
#' @param x  A numeric vector of data values. Non-finite (e.g.,
#' infinite or missing) values will be omitted.
#' @param y An optional numeric vector of data values.
#' @param cluster numeric or charater vector, the id of clusters.
#'  If not specified, each observation will
#' be assigned a distinct cluster, i.e., no cluster in the data.
#'@param data A optional data frame.
#'@param alternative a character string specifying the
#' alternative hypothesis, must be one of "two.sided" (default),
#'  "greater" or "less". You can specify just the initial letter.
#'@param mu a number specifying an optional parameter used to
#'  form the null hypothesis. See 'Details'.
#'@param permutation A logical, whether to use permutation test.
#'@param n.rep number of samples generated for permutation test.
#'@param formula   an object of class \code{"formula"} in the
#'form of  \code{lhs \~ rhs}, where \code{lhs} is a numeric
#'variable giving the data values and \code{rhs} contains
#'the \code{cluster}, \code{group}, \code{stratum}, e.g.,
#'\code{z ~ cluster(a) + group(b) + stratum(c)}, where
#'\code{cluster}, \code{group}, \code{stratum} are special terms.
#'@param subset an optional vector specifyin.csize a
#'subset of observations to be used.
#'@param na.action a function which indicates what should happen
#'when the data contains NAs. The  default action is to omit them.
#'@param group.x a character or a number, indicates which group id
#'  is for treatment x.
#' @param ...  Further arguments to be passed to or from methods.
#' @details
#' THe formula interface is only applicable for the 2-sample tests,
#' and vice versa.
#'
#' If only \code{x} is given, or if both \code{x} and \code{y} are given,
#' then a Wilcoxon signed rank test of the null that the distribution of
#' \code{x} (in the one sample case) or of \code{x - y} (in the paried
#' two sample case) is symmetric about \code{mu} is performed.
#'
#' By default(if \code{permutation} is not specified),
#' a normal approximation is used. Otherwise, a permutation test
#' is used.
#'@return  a list with class "ctest" containing the following components:
#' \item{rstatistic}{the value of the signed rank statistic
#'  with a name describing it.}
#'\item{erstatistics}{Expected value clustered Wilcoxon ranksum statistic.}
#'\item{vrstatistics}{Variance of clustered Wilcoxon ranksum statistic.}
#'\item{statistics}{the value of the test statistic.}
#'\item{p.value}{the p-value for the test}
#'\item{data.name}{a character string giving the names of the data.}
#'\item{method}{the name of the method}
#'\item{balance}{a logical, indicating if the data is balanced.}
#' \item{n}{Total number of observations.}
#' \item{cn}{Total number of clusters.}
#' \item{adjusted}{indicator of whether adjusted signed rank statistic is used.}
#' @seealso \code{\link{cluswilcox.test.formula}}, \code{\link{cluswilcox.test.numeric}}
#' @examples
#' ## Formula interface, only for rank sum test.
#' data(crd)
#' cluswilcox.test(z ~ group(group) + cluster(id), data = crd)
#' data(crdStr)
#' cluswilcox.test(z ~ group(group) + cluster(id) + stratum(stratum), data = crdStr)
#' ## Numeric interface, only for signed rank test.
#' data(crsd)
#' cluswilcox.test(z, cluster = id, data = crsd)
#' data(crsdUnb)
#' cluswilcox.test(z, cluster = id, data = crsdUnb)
#' @export

cluswilcox.test <- function(x, ...) {
  pars <- as.list(match.call())
  if( !exists(as.character(pars$x)) && !is.null(pars$data)) {
    x <- eval(pars$data)[,deparse(substitute(x))]
  }
  UseMethod("cluswilcox.test", x)
}

