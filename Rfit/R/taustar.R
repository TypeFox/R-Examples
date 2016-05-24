#' Estimate of the Scale Parameter taustar
#' 
#' An estimate of the scale parameter taustar = 1/(2*f(0)) is needed for the
#' standard error of the intercept in rank-based regression.
#' 
#' Confidence interval estimate of taustar. See, for example, Hettmansperger
#' and McKean (1998) p.7-8 and p.25-26.
#' 
#' @param resid full model residuals
#' @param p is the number of regression coefficients (without the intercept)
#' @param conf confidence level of CI used
#' @return Length-one numeric object containing the estimated scale parameter
#' taustar.
#' @author Joseph McKean, John Kloke
#' @seealso \code{\link{rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' ##  This is an internal function.  See rfit for user-level examples.
#' 
#' ## The function is currently defined as
#' function (resid, p, conf = 0.95) 
#' {
#'     n = length(resid)
#'     zc = qnorm((1 + conf)/2)
#'     c1 = (n/2) - ((sqrt(n) * zc)/2) - 0.5
#'     ic1 = floor(c1)
#'     if (ic1 < 0) {
#'         ic1 = 0
#'     }
#'     z = sort(resid)
#'     l = z[ic1 + 1]
#'     u = z[n - ic1]
#'     df = sqrt(n)/sqrt(n - p - 1)
#'     taustar = df * ((sqrt(n) * (u - l))/(2 * zc))
#'     taustar
#'   }
#' 
#' @export taustar
taustar <-
function (resid, p, conf = 0.95) 
{
    n = length(resid)
    zc = qnorm((1 + conf)/2)
    c1 = (n/2) - ((sqrt(n) * zc)/2) - 0.5
    ic1 = floor(c1)
    if (ic1 < 0) {
        ic1 = 0
    }
    z = sort(resid)
    l = z[ic1 + 1]
    u = z[n - ic1]
    df = sqrt(n)/sqrt(n - p - 1)
    taustar = df * ((sqrt(n) * (u - l))/(2 * zc))
    taustar
}
