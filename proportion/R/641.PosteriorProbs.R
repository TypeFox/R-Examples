#' Bayesian posterior Probabilities
#' @param n - Number of trials
#' @param a - Prior Parameters
#' @param b - Prior Parameters
#' @param th - Theta value seeking Pr(Theta/X < th)
#' @details  Computes probability of the event \eqn{p < p0} (p0 is specified in [0, 1]) based
#' on posterior distribution of Beta-Binomial model with given parameters for prior Beta
#' distribution for all \eqn{x = 0, 1, 2......n } where \code{n} is the number of trials
#' @return A dataframe with
#'  \item{x}{ Number of successes}
#'  \item{PosProb}{ Posterior probability}
#' @family Miscellaneous  functions for Bayesian method
#' @examples
#' n=5;  a=0.5; b=0.5; th=0.5;
#' probPOS(n,a,b,th)
#' @references
#' [1] 2002 Gelman A, Carlin  JB, Stern HS and Dunson DB
#' Bayesian Data Analysis, Chapman & Hall/CRC
#' [2] 2006  Ghosh M, Delampady M and Samanta T.
#' An introduction to Bayesian analysis: Theory and Methods. Springer, New York
#' @export
#####Bayesian posterior Probabilites Pr(Theta/X < th) Other cases can be obtained from this
probPOS<-function(n,a,b,th)
{
  if (missing(n)) stop("'n' is missing")
  if (missing(a)) stop("'a' is missing")
  if (missing(b)) stop("'b' is missing")
  if (missing(th)) stop("'th' is missing")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(a) != "integer") & (class(a) != "numeric") || length(a) >1|| a<=0 ) stop("'a' has to be greater than 0")
  if ((class(b) != "integer") & (class(b) != "numeric") || length(b) >1|| b<=0 ) stop("'b' has to be greater than 0")
  if ((class(th) != "integer") & (class(th) != "numeric") || length(th) >1|| th<0 ) stop("'theta' has to be greater than x")


####INPUT n
x=0:n
k=n+1
####INITIALIZATIONS
PosProb=0
##############

for(i in 1:k)
{
bet=function(p) dbeta(p,shape1=x[i]+a,shape2=n-x[i]+b)
PosProb[i]=integrate(bet,0,th)$value
}
return(data.frame(x,PosProb))
}
################################################################

