#' @title Sign Test
#' 
#' @description This function will test a hypothesis based on the sign test and reports linearly interpolated confidence intervals for one sample problems.
#' 
#' @details Computes a \dQuote{Dependent-samples Sign-Test} if both \code{x} and \code{y} are provided.  If only \code{x} is provided,  computes the \dQuote{Sign-Test.} 
#' 
#' @param x \code{NA}s and \code{Inf}s are allowed but will be removed 
#' @param y optional numeric vector; \code{NA}s and \code{Inf}s are allowed but will be removed  
#' @param md a single number representing the value of the population median specified by the null hypothesis
#' @param alternative is a character string, one of \code{"greater"}, \code{"less"}, or \code{"two.sided"}, or the initial letter of each, indicating the specification of the alternative hypothesis. For one-sample tests, \code{alternative} refers to the true median of the parent population in relation to the hypothesized value of the median.
#' @param conf.level confidence level for the returned confidence interval, restricted to lie between zero and one
#' 
#' @return A list of class \code{htest}, containing the following components:
#' \item{\code{statistic}}{the S-statistic (the number of positive differences between the data and the hypothesized median), with names attribute \dQuote{S.}}
#' \item{\code{p.value}}{the p-value for the test}
#' \item{\code{conf.int}}{is a confidence interval (vector of length 2) for the true median based on linear interpolation. The confidence level is recorded in the attribute \code{conf.level}. When the alternative is not \code{"two.sided,"} the confidence interval will be half-infinite, to reflect the interpretation of a confidence interval as the set of all values \code{k} for which one would not reject the null hypothesis that the true mean or difference in means is \code{k}. Here infinity will be represented by \code{Inf}.}
#' \item{\code{estimate}}{is a vector of length 1, giving the sample median; this estimates the corresponding population parameter. Component \code{estimate} has a names attribute describing its elements.}
#' \item{\code{null.value}}{is the value of the median specified by the null hypothesis. This equals the input argument \code{md}. Component \code{null.value} has a names attribute describing its elements.}
#' \item{\code{alternative}}{records the value of the input argument alternative: \code{"greater"}, \code{"less"}, or \code{"two.sided"}}
#' \item{\code{data.name}}{a character string (vector of length 1) containing the actual name of the input vector \code{x}}
#' 
#' @section Null Hypothesis:
#' For the one-sample sign-test, the null hypothesis is that the median of the population from which \code{x} is drawn is \code{md}. For the two-sample dependent case, the null hypothesis is that the median for the differences of the populations from which \code{x} and \code{y} are drawn is \code{md}. The alternative hypothesis indicates the direction of divergence of the population median for \code{x} from \code{md} (i.e., \code{"greater"}, \code{"less"}, \code{"two.sided"}.)
#' 
#' @section Assumptions: 
#' The median test assumes the parent population is continuous.
#' 
#' @section Confidence Interval:
#' A linear interpolation is returned for the related confidence interval (returned component \code{conf.int}), which can be obtained by interpolating between the possible achieved confidence levels closest to the desired level. Note that, as explained under the description of \code{conf.int}, the confidence interval will be half-infinite when alternative is not \code{"two.sided"}; infinity will be represented by \code{Inf}.
#' 
#' @references \itemize{
#' \item Gibbons, J.D. and Chakraborti, S. 1992. \emph{Nonparametric Statistical Inference}. Marcel Dekker Inc., New York. 
#' \item Kitchens, L.J. 2003. \emph{Basic Statistics and Data Analysis}. Duxbury. 
#' \item Conover, W. J. 1980. \emph{Practical Nonparametric Statistics, 2nd ed}. Wiley, New York. 
#' \item Lehmann, E. L. 1975. \emph{Nonparametrics: Statistical Methods Based on Ranks}. Holden and Day, San Francisco.
#' }
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#' 
#' @section Note:
#' The reported confidence interval is based on linear interpolation. The lower and upper confidence levels are exact.
#' 
#' @seealso \code{\link{z.test}}, \code{\link{zsum.test}}, \code{\link{tsum.test}}
#' @export
#'  
#' @examples
#' with(data = PHONE, SIGN.test(call.time, md = 2.1))
#' # Computes two-sided sign-test for the null hypothesis
#' # that the population median is 2.1.  The alternative
#' # hypothesis is that the median is not 2.1.  An interpolated
#' # upper 95% upper bound for the population median will be computed. 
#'
#' @keywords htest
#' 
########################################################################### 
# Performs one sample and two sample (dependent) Sign-Test
# for the Median as well as computing a corresponding
# confidence iterval for the median.  Author:  Alan T. Arnholt
# Last Revision: 06/12/13
###########################################################################
SIGN.test <- function(x, y = NULL, md = 0, alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
{
  alternative <- match.arg(alternative)
  if(!missing(md) && (length(md) != 1 || is.na(md)))
    stop("median must be a single number")
  if(!missing(conf.level) && (length(conf.level) != 1 || is.na(conf.level) || conf.level < 0 || conf.level > 1))
    stop("conf.level must be a number between 0 and 1")
###########################################################################
# One-Sample Sign-Test Exact Test  
  if( is.null(y) )
  {
    dname <- paste(deparse(substitute(x)))
    # fix to agree with documentation 6/12/13 - ie remove Infs and NAs
    x <- x[is.finite(x)]
    x <- sort(x)
    diff <- (x - md)
    n <- length(x)
    nt <- length(x) - sum(diff == 0)
    s <- sum(diff > 0)
    estimate <- median(x)
    method <- c("One-sample Sign-Test")
    names(estimate) <- c("median of x")
    names(md) <- "median"
    names(s) <- "s"
    CIS <- "Conf Intervals"
    if(alternative == "less")
    {
      # zobs <- (s-0.5*n)/sqrt(n*0.25)
      pval <- sum(dbinom(0:s, nt, 0.5))
      # Note: Code uses linear interpolation to arrive at the confidence intervals.
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if(k < 1)
      {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- -Inf
        xu <- x[n]
        ici <- c(xl, xu)
      }
      else
      {
        ci1 <- c(-Inf, x[n - k + 1])
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(-Inf, x[n - k])
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- -Inf
        xu <- (((x[n - k + 1] - x[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + x[n - k]
        ici <- c(xl, xu)
      } 
    }
    else if(alternative == "greater")
    {
      pval <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if(k < 1)
      {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- x[1]
        xu <- Inf
        ici <- c(xl, xu)
      }
      else
      {
        ci1 <- c(x[k], Inf)
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(x[k + 1], Inf)
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- (((x[k] - x[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + x[k + 1]
        xu <- Inf
        ici <- c(xl, xu)
      }
    }
    else
    {
      p1 <- sum(dbinom(0:s, nt, 0.5))
      p2 <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      pval <- min(2 * p1, 2 * p2, 1)
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)/2][1]
      if(k < 1)
      {
        conf.level <- (1 - 2 * (sum(dbinom(k, n, 0.5))))
        xl <- x[1]
        xu <- x[n]
        ici <- c(xl, xu)
      }
      else
      {
        ci1 <- c(x[k], x[n - k + 1])
        acl1 <- (1 - 2 * (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(x[k + 1], x[n - k])
        acl2 <- (1 - 2 * (sum(dbinom(0:k, n, 0.5))))
        xl <- (((x[k] - x[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + x[k + 1]
        xu <- (((x[n - k + 1] - x[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + x[n - k]
        ici <- c(xl, xu)
      }
    }
  }
  else
  {
    #   Paired-Samples Sign Test
    if(length(x) != length(y))
      stop("Length of x must equal length of y")
    xy <- sort(x - y)
    # remove NAs and Infs - 6/12/13
    xy <- xy[is.finite(xy)]
    diff <- (xy - md)
    n <- length(xy)
    nt <- length(xy) - sum(diff == 0)
    s <- sum(diff > 0)
    dname <-  paste(deparse(substitute(x)), " and ", deparse(substitute(y)), sep = "")
    estimate <- median(xy)
    method <- c("Dependent-samples Sign-Test")
    names(estimate) <- c("median of x-y")
    names(md) <- "median difference"
    names(s) <- "S"
    CIS <- "Conf Intervals"
    if(alternative == "less")
    {
      pval <- sum(dbinom(0:s, nt, 0.5))
      # Note: Code uses linear interpolation to arrive at the confidence intervals.
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if(k < 1)
      {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- -Inf
        xu <- xy[n]
        ici <- c(xl, xu)
      }
      else
      {
        ci1 <- c(-Inf, xy[n - k + 1])
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(-Inf, xy[n - k])
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- -Inf
        xu <- (((xy[n - k + 1] - xy[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + xy[n - k]
        ici <- c(xl, xu)
      }
    }
    else if(alternative == "greater")
    {
      pval <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)][1]
      if(k < 1)
      {
        conf.level <- (1 - (sum(dbinom(k, n, 0.5))))
        xl <- xy[1]
        xu <- Inf
        ici <- c(xl, xu)
      }
      else
      {
        ci1 <- c(xy[k], Inf)
        acl1 <- (1 - (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(xy[k + 1], Inf)
        acl2 <- (1 - (sum(dbinom(0:k, n, 0.5))))
        xl <- (((xy[k] - xy[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + xy[k + 1]
        xu <- Inf
        ici <- c(xl, xu)
      }
    }
    else
    {
      p1 <- sum(dbinom(0:s, nt, 0.5))
      p2 <- (1 - sum(dbinom(0:s - 1, nt, 0.5)))
      pval <- min(2 * p1, 2 * p2, 1)
      loc <- c(0:n)
      prov <- (dbinom(loc, n, 0.5))
      k <- loc[cumsum(prov) > (1 - conf.level)/2][1]
      if(k < 1)
      {
        conf.level <- (1 - 2 * (sum(dbinom(k, n, 0.5))))
        xl <- xy[1]
        xu <- xy[n]
        ici <- c(xl, xu)
      }
      else
      {
        ci1 <- c(xy[k], xy[n - k + 1])
        acl1 <- (1 - 2 * (sum(dbinom(0:k - 1, n, 0.5))))
        ci2 <- c(xy[k + 1], xy[n - k])
        acl2 <- (1 - 2 * (sum(dbinom(0:k, n, 0.5))))
        xl <- (((xy[k] - xy[k + 1]) * (conf.level - acl2))/(acl1 - acl2)) + xy[k + 1]
        xu <- (((xy[n - k + 1] - xy[n - k]) * (conf.level - acl2))/(acl1 - acl2)) + xy[n - k]
        ici <- c(xl, xu)
      }
    }
  } 
  # Below is how I had the stuff print at first...showing the exact intervals as
  # well as the interpolated interval.  Changed format to standard htest class.
  # Test.Values <- t(as.matrix(c(s, pval)))
  # values <- c("S", "p-value")
  # dimnames(Test.Values) <- list(NULL, values)
  if(k < 1)
  {
    cint <- ici
    attr(cint, "conf.level") <- conf.level
    rval <- structure(list(statistic = s, p.value = pval, estimate = estimate, null.value = md,
                           alternative = alternative, method = method, data.name = dname, conf.int = cint ))
    oldClass(rval) <- "htest"
    return(rval)
  }
  else
  {
    result1 <- c(acl2, ci2)
    result2 <- c(conf.level, ici)
    result3 <- c(acl1, ci1)
    Confidence.Intervals <- round(as.matrix(rbind(result1, result2, result3)),4)
    cnames <- c("Conf.Level", "L.E.pt", "U.E.pt")
    rnames <- c("Lower Achieved CI", "Interpolated CI", "Upper Achieved CI")
    dimnames(Confidence.Intervals) <- list(rnames, cnames)
    # return(Test.Values, Confidence.Intervals)
    cint <- ici
    
    attr(cint, "conf.level") <- conf.level
    
    rval <- structure(list(statistic = s, parameter = NULL, p.value = pval,
                           conf.int = cint, estimate = estimate, null.value = md,
                           alternative = alternative, method = method, data.name = dname ))
    ## Returns both htest format and matrix of CIs
    ## attr(rval, "class") <- "htest"
    oldClass(rval) <- "htest"
    ##   return(rval, Confidence.Intervals)
    # return(rval)
    print(rval)
    return(Confidence.Intervals)
  }
}