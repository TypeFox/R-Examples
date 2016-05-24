#' Compute vectors measuring stochastic dominance of four orders.
#' 
#' Stochastic dominance is a sophisticated comparison of two distributions of
#' stock market returns.  The dominating distribution is superior in terms of
#' mean, variance, skewness and kurtosis respectively, representing dominance
#' orders 1 to 4, without directly computing four moments.  Vinod(2008) sec. 4.3
#' explains the details.  This function directly uses the output of `wtdpapb'.
#' 
#' 
#' @param dj these are (unequal) distances of consecutive intervals
#' @param wpa weighted probabilities of first set of probabilities
#' @param wpb weighted probabilities of second set of probabilities
#' @return 
#' \item{sd1b}{vector measuring stochastic dominance of order 1, SD1} 
#' \item{sd2b}{vector measuring stochastic dominance of order 2, SD2} 
#' \item{sd3b}{vector measuring stochastic dominance of order 3, SD3} 
#' \item{sd4b}{vector measuring stochastic dominance of order 4, SD4} 
## @note %% ~~further notes~~
#' @author Prof. H. D. Vinod, Economics Dept., Fordham University, NY
#' @seealso See Also \code{\link{wtdpapb}}
#' @references Vinod, H. D.", "Hands-On Intermediate Econometrics 
#'  Using R"  (2008) World Scientific Publishers: Hackensack, NJ.
#'  \url{http://www.worldscibooks.com/economics/6895.html}
#'  
#'
#' @references Vinod, H. D. "Ranking Mutual Funds Using 
#' Unconventional Utility Theory and Stochastic Dominance,"
#' Journal of Empirical Finance Vol. 11(3) 2004, pp. 353-377.
#' 
#' @keywords mean variance skewness kurtosis
#' @examples
#' 
#'  \dontrun{
#'  set.seed(234);x=sample(1:30);y=sample(5:34)
#'  w1=wtdpapb(x,y)
#'  stochdom2(w1$dj, w1$wpa, w1$wpb) }
#' 
#' @export

stochdom2 <-
function( dj, wpa, wpb) {
# input weighted pa and pb  IsubF and I sub f matrices
if (length(wpa) != length(dj)) print("wrong rows dj")
if (length(wpa) != length(wpb)) print("wrong rows wpa wpb")
rhs=cumsum(wpa-wpb)
sd1b=bigfp(d=dj, p=rhs)
#sd1=I.bigf %*% I.smallf %*% (wpa-wpb)
sd2b=bigfp(d=dj, p=sd1b)
sd3b=bigfp(d=dj, p=sd2b)
sd4b=bigfp(d=dj, p=sd3b)
list( sd1b=sd1b,sd2b=sd2b, sd3b=sd3b,sd4b=sd4b) }
