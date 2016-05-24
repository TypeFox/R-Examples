#' Computes CDF of KI between case profile and profile with stated relationship
#' 
#' Computes the Cumulative Distribution Function of a Kinship Index (KI) comparing hypotheses \code{hyp.1} vs \code{hyp.2} for profiles with a given relationship (\code{hyp.true}) to the case profile (e.g. \code{"FS"} for full siblings).
#' 
#' @param x (optional) An integer matrix specifying a single profile.
#' @param hyp.1 A character string giving the hypothesis in the numerator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param hyp.2 A character string giving the hypothesis in the denominator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param hyp.true A character string specifying the true relationship between the case profile and the other profile. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs.ki A list specifying the allelic frequencies that are used when computing the \eqn{KI}.
#' @param freqs.true (optionally) A list specifying the allelic frequencies that are used for computing the probabily distribution of the \eqn{KI} under \code{hyp.true}. When not provided, the function will use \code{freqs}. One might use different allelic frequencies \code{freqs.rel} when for example the case profile and relative come from some population, while \eqn{KI}s are computed with frequencies from another population.
#' @param markers Character vector stating the markers for which the KI distribution is derived. By default equal to all markers of \code{x}, or, if \code{x} is missing, to the intersection of the markers of \code{freqs.ki} and \code{freqs.true}.
#' @param theta.ki numeric value specifying the amount of background relatedness.
#' @param theta.true numeric value specifying the amount of background relatedness.
#' @param n.max Maximum number of events stored in memory. See \code{dists.product.pair} for details.
#' @examples
#' # for one profile, obtain the CDF of the SI,
#' # both for true sibs and unrelated profiles
#' data(freqsNLsgmplus)
#'
#' # sample a profile
#' x <- sample.profiles(N=1,freqsNLsgmplus)
#'
#' cdf.fs <- ki.cdf(x,hyp.1="FS",hyp.true="FS")
#' cdf.un <- ki.cdf(x,hyp.1="FS",hyp.true="UN")
#'
#' # the cdf's are *functions*
#' cdf.fs(1)
#' cdf.un(1)
#'
#' # we also obtain an ROC curve easily
#' t <- 10^(seq(from=-10,to=10,length=100)) # some thresholds
#' fpr <- cdf.un(t,exc.prob=TRUE)
#' tpr <- cdf.fs(t,exc.prob=TRUE)
#'
#' plot(log10(fpr),tpr,type="l")
#' @export
ki.cdf <- function(x,hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki=get.freqs(x),freqs.true=freqs.ki,markers={if (missing(x)) intersect(names(freqs.ki),names(freqs.true)) else get.markers(x)},theta.ki=0,theta.true=theta.ki,n.max=1e6){      
  if (missing(x)){
    stop("unconditional ki.cdf not yet implemented. Use ki.dist instead.")
  }else{
    x <- Zassure.matrix(x) 
    # obtain the cond ki dist for all markers
    x.cond.ki.dist <- ki.dist(x=x,hyp.1=hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,freqs.ki=freqs.ki,freqs.true=freqs.true,markers = markers,theta.ki=theta.ki,theta.true=theta.true)
    # return a nice function
    return(dist.pair.cdf(dists.product.pair(x.cond.ki.dist,n.max=n.max)))
  }
}