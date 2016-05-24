#' binom.confint
#'
#' Exact confidence intervals for a binomial proportion parameter.
#'
#' @param k #"heads"
#' @param n sample size
#' @param alpha Confidence level
#' @param side Sidedness of the hypothesis: c("two", "upper", "lower")

binom.confint =
  function(k, n, alpha=0.05,
           side=c("two", "upper", "lower")) {

    sideLetter = toupper(substring(side[1], 1, 1))
    if(sideLetter=="T") {
      alphaUpper = alphaLower = alpha/2
    } else if(sideLetter=="U") {
      alphaUpper = alpha; alphaLower = 0
    } else if(sideLetter=="L") {
      alphaLower = alpha; alphaUpper = 0
    }
    else stop("bad value for \"side\" ")
    huntForBoundaryUpper = function(p)
      alphaUpper - (1 - pbinom((k-1),n,p))
    lb = try(uniroot(huntForBoundaryUpper, interval=0:1)$root)
    if(class(lb) == "try-error") lb = 0
    huntForBoundaryLower = function(p)
      alphaLower - pbinom(k,n,p)
    ub = try(uniroot(huntForBoundaryLower, interval=0:1)$root)
    if(class(ub) == "try-error") ub = 1
    return(c(lb=lb, ub=ub))
  }
