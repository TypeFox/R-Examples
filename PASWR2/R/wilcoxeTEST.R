#' @title Wilcoxon Exact Test
#' 
#' @description Performs exact one sample and two sample Wilcoxon tests on vectors of data
#' 
#' @details If only \code{x} is given, or if both \code{x} and \code{y} are given and \code{paired = TRUE}, a Wilcoxon signed rank test of the null hypothesis that the distribution of \code{x} (in the one sample case) or of \code{x - y} (in the paired two sample case) is symmetric about \code{mu} is performed.  
#' 
#' Otherwise, if both \code{x} and \code{y} are given and \code{paired = FALSE}, a Wilcoxon rank sum test is done.  In this case, the null hypothesis is that the distribution of \code{x} and \code{y} differ by a location shift \code{mu}, and the alternative is that they differ by some other location shift (and the one-sided alternative \code{"greater"} is that \code{x} is shifted to the right of \code{y}).
#' 
#' @param x is a numeric vector of data values.  Non-finite (i.e. infinite or missing) values will be omitted.
#' @param y an optional numeric vector of data values
#' @param mu a number specifying an optional parameter used to form the null hypothesis
#' @param paired a logical indicating whether you want a paired test
#' @param alternative a character string specifying the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.  You can specify just the initial letter.
#' @param conf.level confidence level of the interval
#' 
#' @return A list of class \code{htest}, containing the following components:
#' \item{\code{statistic}}{the value of the test statistic with a name describing it}
#' \item{\code{p.value}}{the p-value for the test}
#' \item{\code{null.value}}{the location parameter \code{mu}}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis}
#' \item{\code{method}}{the type of test applied}
#' \item{\code{data.name}}{a character string giving the names of the data}
#' \item{\code{conf.int}}{a confidence interval for the location parameter}
#' \item{\code{estimate}}{an estimate of the location parameter}
#' 
#' @section Note:
#' The function is rather primitive and should only be used for problems with fewer than 19 observations as the memory requirements are rather large.
#'
#' @references \itemize{
#'  \item Gibbons, J.D. and Chakraborti, S. 1992. \emph{Nonparametric Statistical Inference}. Marcel Dekker Inc., New York. 
#'  \item Hollander, M. and Wolfe, D.A. 1999. \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. 
#'  }
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#'  
#' @seealso \code{\link{wilcox.test}}
#' @export
#'  
#' @examples
#'  
#' # Wilcoxon Signed Rank Test
#' PH <- c(7.2, 7.3, 7.3, 7.4)
#' wilcoxe.test(PH, mu = 7.25, alternative = "greater")
#' # Wilcoxon Signed Rank Test (Dependent Samples)
#' with(data = AGGRESSION, 
#' wilcoxe.test(violence, noviolence, paired = TRUE, alternative = "greater"))
#' # Wilcoxon Rank Sum Test
#' x <- c(7.2, 7.2, 7.3, 7.3)
#' y <- c(7.3, 7.3, 7.4, 7.4)
#' wilcoxe.test(x, y)
#' rm(PH, x, y)
#'  
#' @keywords htest                                          
#####################################################################################
wilcoxe.test <- function(x, y = NULL, mu = 0, paired = FALSE, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
# Computes Exact distribution of Wilcoxon Signed Rank Test and Wilcoxon Sum 
# Rank Test.  In the presence of ties, computes the exact conditional 
# distribution of the Wilcoxon Signed Rank and the Wilcoxon Sum Rank Test.  
# For the Wilcoxon Signed Rank test, if any Di = 0, the procedure stops.
# Author: Alan T. Arnholt -- Last Revised: 06/18/13  
{
  alternative <- match.arg(alternative)
  if (!missing(mu) && ((length(mu) > 1) || !is.finite(mu)))
    stop("mu must be a single number")
  if (!((length(conf.level) == 1) && is.finite(conf.level) &&  (conf.level > 0) && (conf.level < 1)))
    stop("conf.level must be a single number between 0 and 1")
  if (!is.numeric(x))
    stop("x must be numeric")
  if (!is.null(y))
  {
    if (!is.numeric(y))
      stop("y must be numeric")
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (paired==TRUE) {
      if (length(x) != length(y))
        stop("x and y must have the same length")
    }
    else {
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
    }
  }
  
  else {
    dname <- deparse(substitute(x))
    if (paired==TRUE)
      stop("y missing for paired test")
    x <- x[is.finite(x)]
  }
  if (length(x) < 1)
    stop("not enough (finite) x observations")
  
  if (is.null(y) && paired==FALSE)
  {
    
    if (length(x) > 19)
      stop(" Sample size too large (n > 19).  Use wilcox.test")
    
    values <- x
    D1 <- values - mu
    D <- D1[D1 != 0]
    
    n.pitch <- length(D1) - length(D)
    if (n.pitch > 0)
      stop(" Not possible to get exact values with zero differences. Use wilcox.test")
    
    absD <- abs(D)
    rankabsD <- rank(absD)
    signD <- sign(D)
    signrank <- rankabsD*signD
    stat <- sum(signrank[signrank > 0])
    SWA <- outer(values, values, "+")
    SWA <- sort(SWA[!lower.tri(SWA)]) / 2     # The n*(n+1)/2 Walsh Averages
    HLE <- median(SWA)                        # Hodges-Lehman Estimator
    
############## A more pedagogical approach ##############################
# n2means <- apply(srs(D,2),1,mean)   # Computing the n choose 2 means  #
# WalshAverages <- c(D,n2means)       # The n*(n+1)/2 Walsh Averages    #
# SWA <- sort(WalshAverages)                                            #
#########################################################################
    
    estimate <- HLE
    method <- c("Wilcoxon Signed Rank Test")
    names(mu) <- "median"
    names(estimate) <- c("(pseudo)median")
    names(stat) <- "t+"
    CIS <- "Conf Intervals"
  }
  else if (!is.null(y) && paired == TRUE)
  {
    
    if (length(x) > 19)
      stop(" Dependent sample size too large (n > 19).  Use wilcox.test")
    
    values <- x - y
    values <- values[is.finite(values)]
    D1 <- values - mu
    D <- D1[D1 != 0]
    
    n.pitch <- length(D1)-length(D)
    if (n.pitch > 0)
      stop(" Not possible to get exact values with zero differences. Use wilcox.test")
    
    absD <- abs(D)
    rankabsD <- rank(absD)
    signD <- sign(D)
    signrank <- rankabsD*signD
    stat <- sum(signrank[signrank > 0])
    SWA <- outer(values, values, "+")
    SWA <- sort(SWA[!lower.tri(SWA)]) / 2     # The n*(n+1)/2 Walsh Averages
    HLE <- median(SWA)                        # Hodges-Lehman Estimator
    estimate <- HLE
    method <- c("Wilcoxon Signed Rank Test (Dependent Samples)")
    names(mu) <- "median difference"
    names(estimate) <- c("(pseudo)median")
    names(stat) <- "t+"
    CIS <- "Conf Intervals"
  }
  
###########################################################################
  if ( (!is.null(y) && paired == TRUE) || (is.null(y) && paired == FALSE) )
  {
    rd <- rank(sort(abs(D)))             # Sorted Ranked abs(D_i)
    n <- length(D)
    r <- 2^n
    mat <- rep(1, r)%*%t(rd)                    # Matrix 2^n rows of the n ranks
    signs <- matrix(1, nrow = r, ncol = n)     # (2^n * n) matrix of 1s
    i <- 0
    while(i < n)
    {
      k <- 2^i
      signs[ , i+1] <- rep(c(rep(1, k), rep(0, k)), 2^(n-i-1))
      i <- i + 1
    }
    ms <- mat*signs
    Tp <- apply(ms, 1, sum)
    MS <-cbind(ms, Tp)
    Tp <- sort(Tp)    
###########################################################################
    
    if (alternative == "less")
    {
      pval <- sum(Tp <= stat)/2^n
      alpha <- 1 - conf.level
      k <- alpha * (2^n) + 1      # Note: k is used here as the order statistic of T+ not the
      k <- max(floor(k), 2)      # value of the distribution of T+ such that P(T+<k)<=alpha/2
      LV <- Tp[k]
      UV <- Tp[2^n - k + 1]
      ci <- c(-Inf, SWA[UV + 1])
      ECL <- 1 - (sum(Tp < Tp[k]) / 2^n)
    }
    else if (alternative == "greater")
    {
      pval <- sum(Tp >= stat) / 2^n
      alpha <- 1- conf.level
      k <- alpha * (2^n) + 1      # Note: k is used here as the order statistic of T+ not the
      k <- max(floor(k), 2)       # value of the distribution of T+ such that P(T+<k)<=alpha/2
      LV <- Tp[k]
      UV <- Tp[2^n - k + 1]
      ci <- c(SWA[LV], Inf)
      ECL <- 1 - (sum(Tp < Tp[k]) / 2^n)
    }
    else
    {
      pval <- 2 * min(0.5, sum(Tp <= stat) / 2^n, sum(Tp >= stat) / 2^n )
      alpha <- 1- conf.level
      k <- alpha/2 * (2^n) + 1       # Note: k is used here as the order statistic of T+ not the
      k <- max(floor(k), 2)          # value of the distribution of T+ such that P(T+<k)<=alpha/2
      LV <- Tp[k]
      UV <- Tp[2^n - k + 1]
      ci <- c(SWA[LV], SWA[UV + 1])
      ECL <- 1 - 2*(sum(Tp < Tp[k]) / 2^n)
    }
  }
  
  else {
    if (length(y) < 1)
      stop("not enough y observations")
    
    method <- "Wilcoxon Rank Sum Test"
    names(mu) <- "median"
    
    n.x <- as.double(length(x))
    n.y <- as.double(length(y))
    N <- n.x + n.y
    
    if (choose(N, n.x) > 1000000)
      stop(" Sample sizes too large (choose(n.x+n.y, n.x) > 1,000,000).  Use wilcox.test")
    
    r <- rank(c(x - mu, y))
    STAT <- sum(r[seq(along = x)]) - n.x*(n.x + 1) / 2
    stat <- sum(r[seq(along = x)])
    W <- apply(srs(r, n.x), 1, sum)
    W <- sort(W)
    U <- W - n.x*(n.x + 1) / 2
    diffs <- sort(outer(x, y, "-"))
    
###########################################################################
    
    rci <- rank(c(x, y))
    if(mu == 0)
    {
      UCI <- U
    }
    else
    {
      WCI <- apply(srs(rci, n.x), 1, sum)
      WCI <- sort(WCI)
      UCI <- WCI - n.x*(n.x + 1) / 2
    }
    
###########################################################################
    
    estimate <- median(diffs)
    names(estimate) <- "difference in location"
    names(stat) <- "w"
    CIS <- "Conf Intervals"
    
    if (alternative == "less")
    {
      pval <- sum(U <= STAT) / length(U)
      alpha <- 1 - conf.level
      k <- max(floor(alpha*choose(N, n.x) + 1), 2)
      o.s. <- rank(unique(UCI))[unique(UCI) == UCI[k]]
      ci <- c(-Inf, diffs[n.x*n.y - o.s. + 1])
      ECL <- 1 - 2*sum(U <= U[k])/choose(N, n.x)
    }
    
    else if (alternative == "greater")
    {
      pval <- sum(U >= STAT) / length(U)
      alpha <- 1 - conf.level
      k <- max(floor(alpha*choose(N, n.x) + 1), 2)
      o.s. <- rank(unique(UCI))[unique(UCI) == UCI[k]]
      ci <- c(diffs[o.s.], Inf)
      ECL <- 1 - 2*sum(U <= U[k]) / choose(N, n.x)
    }
    else
    {
      pval <- 2 * min( 0.5, sum(U <= STAT) / length(U), sum(U >= STAT) / length(U) )
      alpha <- 1 - conf.level
      k <- max(floor(alpha / 2*choose(N, n.x) + 1), 2)
      o.s. <- rank(unique(UCI))[unique(UCI) == UCI[k]]
      ci <- c(diffs[o.s.], diffs[n.x*n.y - o.s. + 1])
      ECL <- 1 - 2*sum(U <= U[k]) / choose(N, n.x)
    }
    
    
  }
  
  cint <- ci
  attr(cint, "conf.level") <- ECL
  ans <- structure(list(statistic = stat, p.value = pval, estimate = estimate,
                        null.value = mu, alternative = alternative, method = method, 
                        data.name = dname, conf.int=cint ) )
  oldClass(ans) <- "htest"
  return(ans)
}