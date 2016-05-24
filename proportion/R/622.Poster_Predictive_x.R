#' The Predicted probability - Bayesian approach
#' @param x - Number of successes
#' @param n - Number of trials from data
#' @param xnew - Required size of number of success
#' @param m - Future :Number of trials
#' @param a1 - Beta Prior Parameters for Bayesian estimation
#' @param a2 - Beta Prior Parameters for Bayesian estimation
#' @details  Computes posterior predictive probability for the required size of number of
#' successes  for \code{xnew} of \code{m} trials from the given number of successes \code{x}
#'  of \code{n} trials for the given parameters for Beta prior distribution
#' @return A dataframe with x,n,xnew,m,preprb
#'  \item{x}{ Number of successes}
#'  \item{n}{ Number of trials}
#'  \item{xnew}{ Required size of number of success}
#'  \item{m}{ Future - success, trails}
#'  \item{preprb }{ The predicted probability}
#' @family Miscellaneous  functions for Bayesian method
#' @examples
#' x=0; n=1; xnew=10; m=10; a1=1; a2=1
#' probPREx(x,n,xnew,m,a1,a2)
#' @references
#' [1] 2002 Gelman A, Carlin  JB, Stern HS and Dunson DB
#' Bayesian Data Analysis, Chapman & Hall/CRC
#' @export
probPREx<-function(x,n,xnew,m,a1,a2)
{
  if (missing(x)) stop("'x' is missing")
  if (missing(n)) stop("'n' is missing")
  if (missing(xnew)) stop("'xnew' is missing")
  if (missing(m)) stop("'m' is missing")
  if (missing(a1)) stop("'a1' is missing")
  if (missing(a2)) stop("'a2' is missing")
  if ((class(x) != "integer") & (class(x) != "numeric") || length(x) >1|| x<0 || x>n) stop("'x' has to be between 0 and n")
  if ((class(n) != "integer") & (class(n) != "numeric") || length(n) >1|| n<0 ) stop("'n' has to be greater or equal to 0")
  if ((class(xnew) != "integer") & (class(xnew) != "numeric") || length(xnew) >1 || xnew>m ) stop("'xnew' has to be between 0 and m")
  if ((class(m) != "integer") & (class(m) != "numeric") || length(m) >1|| m<0 ) stop("'m' has to be greater or equal to 0")
  if ((class(a1) != "integer") & (class(a2) != "numeric") || length(a1) >1|| a1<0 ) stop("'a1' has to be greater or equal to 0")
  if ((class(a2) != "integer") & (class(a2) != "numeric") || length(a2) >1|| a2<0 ) stop("'a2' has to be greater or equal to 0")

    preprb=(choose(m,xnew))/(beta(x+a1,n-x+a2))*beta(xnew+x+a1,m+n-xnew-x+a2)
    return(data.frame(x,n,xnew,m,preprb))
  }
