#' Intermediate weighting function giving Non-Expected Utility theory weights.
#'
#' Computes cumulative probabilities and difference between consecutive
#' cumulative probabilities described in Vinod (2008) textbook.  This is a simpler version
#' of the version in the book without mapping to non-expected utility theory weights.
#'
#' @param n {A (usually small) integer.}
#' @return
#' \item{x}{sequence 1:n}
#' \item{p}{probabilities p= x[i]/n}
#' \item{pdif}{consecutive differences p[i] - p[i - 1]}
###' ##@note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
####@seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Vinod, H. D.", "Hands-On Intermediate Econometrics
#' Using R"  (2008) World Scientific Publishers: Hackensack, NJ.
#' \url{http://www.worldscibooks.com/economics/6895.html}
#' @keywords Prelec
#' @examples
#'
#'  \dontrun{prelec2(10)}
#'
#' ## The function is currently defined as
#' function (n)
#' {
#'     x = 1:n
#'     p = rep(0, length(x))
#'     for (i in 1:n) {
#'         p[i] = x[i]/n
#'     }
#'     pdif = p
#'     for (i in 2:n) {
#'         pdif[i] = p[i] - p[i - 1]
#'     }
#'     list(x = x, p = p, pdif = pdif)
#'   }
#' @export

prelec2 <-
function(n) #
#input: size n
{x=1:n; p=rep(0,length(x));
#loop to calculate cumulative probabilities
for (i in 1:n) {p[i]=x[i]/n}
pdif=p # to get the first value of wpdif[1] as wp[1]
for (i in 2:n) {
pdif[i]=p[i]-p[i-1]} #compute first differences
list(x=x,p=p,pdif=pdif)  }
