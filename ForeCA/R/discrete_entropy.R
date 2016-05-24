#' @title Shannon entropy for discrete pmf
#' 
#' @description
#' Computes the Shannon entropy \eqn{\mathcal{H}(p) = -\sum_{i=1}^{n} p_i \log p_i}
#' of a discrete RV \eqn{X} taking
#' values in \eqn{\lbrace x_1, \ldots, x_n \rbrace} with probability
#' mass function (pmf) \eqn{P(X = x_i) = p_i} with
#' \eqn{p_i \geq 0} for all \eqn{i} and \eqn{\sum_{i=1}^{n} p_i = 1}.
#' 
#' @details
#' \code{discrete_entropy} uses a plug-in estimator (\code{method = "MLE"}): 
#' \deqn{ 
#' \widehat{\mathcal{H}}(p) = - \sum_{i=1}^{n} \widehat{p}_i \log \widehat{p}_i. 
#' }
#' 
#' If \code{prior.weight > 0}, then it mixes the observed proportions \eqn{\widehat{p}_i}
#'  with a prior distribution
#' \deqn{ 
#' \widehat{p}_i \leftarrow (1-\lambda) \cdot \widehat{p_i} + \lambda \cdot prior_i, \quad i=1, \ldots, n,
#' }
#' where \eqn{\lambda \in [0, 1]} is the \code{prior.weight} parameter.  By default
#' the prior is a uniform distribution, i.e., \eqn{prior_i = \frac{1}{n}} for all i.
#' 
#' Note that this plugin estimator is biased. See References for an overview of alternative
#' methods.
#' 
#' @param probs numeric; probabilities (empirical frequencies). Must be non-negative and add up to \eqn{1}.
#' @inheritParams common-arguments
#' @param method string; method to estimate entropy; see Details below.
#' @param threshold numeric; frequencies below \code{threshold} are set to \eqn{0};
#' default \code{threshold = 0}, i.e., no thresholding.  
#' If \code{prior.weight > 0} then thresholding will be done \emph{before} smoothing.
#' @param prior.probs optional; only used if \code{prior.weight > 0}. 
#' Add a prior probability distribution to \code{probs}. By default it uses a 
#' uniform distribution putting equal probability on each outcome.
#' @param prior.weight numeric; how much weight does the prior distribution get in a mixture
#' model between data and prior distribution? Must be between \code{0} and \code{1}.
#' Default: \code{0} (no prior).
#' @return 
#' numeric; non-negative real value.
#' @keywords math univar
#' @seealso \code{\link{continuous_entropy}}
#' @references
#' Archer E., Park I. M., Pillow J.W. (2014). \dQuote{Bayesian Entropy Estimation for
#' Countable Discrete Distributions}. Journal of Machine Learning Research (JMLR) 15, 
#' 2833-2868. Available at \url{jmlr.org/papers/v15/archer14a.html}.
#' 
#' @export
#' @examples
#' 
#' probs.tmp <- rexp(5)
#' probs.tmp <- sort(probs.tmp / sum(probs.tmp))
#' 
#' unif.distr <- rep(1/length(probs.tmp), length(probs.tmp))
#' 
#' matplot(cbind(probs.tmp, unif.distr), pch = 19, 
#'         ylab = "P(X = k)", xlab = "k")
#' matlines(cbind(probs.tmp, unif.distr))
#' legend("topleft", c("non-uniform", "uniform"), pch = 19, 
#'        lty = 1:2, col = 1:2, box.lty = 0)
#' 
#' discrete_entropy(probs.tmp)
#' # uniform has largest entropy among all bounded discrete pmfs 
#' # (here = log(5))
#' discrete_entropy(unif.distr)
#' # no uncertainty if one element occurs with probability 1
#' discrete_entropy(c(1, 0, 0)) 
#' 
discrete_entropy <- function(probs, base = 2, method = c("MLE"), 
                             threshold = 0, 
                             prior.probs = rep(1 / length(probs), length = length(probs)),
                             prior.weight = 0) {
  
  stopifnot(!any(is.na(probs)),
            prior.weight >= 0 && prior.weight <= 1)
  
  if (!all(round(probs, 6) >= 0)) {
    stop("Not all probabilities are non-negative.")
  }
  
  stopifnot(all.equal(target = 1, current = sum(probs)),
            all(round(prior.probs, 6) >= 0),
            all.equal(target = 1, current = sum(prior.probs)))
  
  method <- match.arg(method)
  
  # set probabilities to 0 that are below threshold (and renormalize)
  if (threshold > 0) {
    probs[probs < threshold] <- 0
    probs <- probs / sum(probs)
  }
  
  if (prior.weight > 0) {
    # add prior
    probs <- (1 - prior.weight) * probs + prior.weight * prior.probs
  }

  if (any(probs == 0)) {
    probs <- probs[probs != 0]
  }
  stopifnot(all.equal(target = 1, current = sum(probs)))
  switch(method,
         MLE = {
           entropy.eval <- -sum(probs * log(probs, base = base))
         })
  attr(entropy.eval, "base") <- as.character(base)
  return(entropy.eval)
} 