#' Walsh Averages
#' 
#' Given a list of n numbers, the Walsh averages are the \eqn{latex}{ n(n+1)/2
#' } pairwise averages.
#' 
#' 
#' @param x A numeric vector
#' @return The Walsh averages.
#' @author John Kloke \email{kloke@@biostat.wisc.edu}
#' @seealso \code{\link{signedrank}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Hollander, M. and Wolfe, D.A. (1999), \emph{Nonparametric Statistical
#' Methods}, New York: Wiley.
#' @examples
#' 
#' 
#' median(walsh(rnorm(100)))  # Hodges-Lehmann estimate of location
#' 
#' ## The function is currently defined as
#' function (x) 
#' {
#'     n <- length(x)
#'     w <- vector(n * (n + 1)/2, mode = "numeric")
#'     ind <- 0
#'     for (i in 1:n) {
#'         for (j in i:n) {
#'             ind <- ind + 1
#'             w[ind] <- 0.5 * (x[i] + x[j])
#'         }
#'     }
#'     return(w)
#'   }
#' 
#' @export walsh
walsh <- function (x) {
    n <- length(x)
    w <- vector(n * (n + 1)/2, mode = "numeric")
    ind <- 0
    for (i in 1:n) {
        for (j in i:n) {
            ind <- ind + 1
            w[ind] <- 0.5 * (x[i] + x[j])
        }
    }
    w
}
