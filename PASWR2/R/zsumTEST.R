#' @title Summarized z-test
#' 
#' @description This function is based on the standard normal distribution and creates confidence intervals and tests hypotheses for both one and two sample problems based on summarized information the user passes to the function. Output is identical to that produced with \code{z.test}.
#' 
#' @details If \code{y} is \code{NULL} , a one-sample z-test is carried out with \code{x} provided \code{sigma.x} is finite.  If y is not \code{NULL}, a standard two-sample z-test is performed provided both \code{sigma.x} and \code{sigma.y} are finite.
#' 
#' @param mean.x a single number representing the sample mean of \code{x}
#' @param sigma.x a single number representing the population standard deviation for \code{x}
#' @param n.x a single number representing the sample size for \code{y}
#' @param mean.y a single number representing the sample mean of \code{y}
#' @param sigma.y a single number representing the population standard deviation for \code{y}
#' @param n.y a single number representing the sample size for \code{y}
#' @param alternative is a character string, one of \code{"greater"}, \code{"less"}, or \code{"two.sided"}, or the initial letter of each, indicating the specification of the alternative hypothesis. For one-sample tests, \code{alternative} refers to the true mean of the parent population in relation to the hypothesized value \code{mu}. For the standard two-sample tests, \code{alternative} refers to the difference between the true population mean for \code{x} and that for \code{y}, in relation to \code{mu}.
#' @param mu a single number representing the value of the mean or difference in means specified by the null hypothesis
#' @param conf.level confidence level for the returned confidence interval, restricted to lie between zero and one
#' @param ... Other arguments passed onto \code{z.test()}
#' 
#' @return A list of class \code{htest}, containing the following components:
#' \item{\code{statistic}}{the z-statistic, with names attribute \code{z}}
#' \item{\code{p.value}}{the p-value for the test}
#' \item{\code{conf.int}}{is a confidence interval (vector of length 2) for the true mean or difference in means. The confidence level is recorded in the attribute \code{conf.level}. When \code{alternative} is not \code{"two.sided,"} the confidence interval will be half-infinite, to reflect the interpretation of a confidence interval as the set of all values \code{k} for which one would not reject the null hypothesis that the true mean or difference in means is \code{k} . Here, infinity will be represented by \code{Inf}.}
#' \item{\code{estimate}}{vector of length 1 or 2, giving the sample mean(s) or mean of differences; these estimate the corresponding population parameters. Component \code{estimate} has a names attribute describing its elements.} 
#' \item{\code{null.value}}{the value of the mean or difference in means specified by the null hypothesis. This equals the input argument \code{mu}. Component \code{null.value} has a names attribute describing its elements.}
#' \item{alternative}{records the value of the input argument alternative: \code{"greater"}, \code{"less"}, or \code{"two.sided"}.}
#' \item{data.name}{a character string (vector of length 1) containing the names x and y for the two summarized samples.}
#' 
#' @section Null Hypothesis:
#'  For the one-sample z-test, the null hypothesis is that the mean of the population from which \code{x} is drawn is \code{mu}. For the standard two-sample z-test, the null hypothesis is that the population mean for \code{x} less that for \code{y} is \code{mu}. 
#'   
#'  The alternative hypothesis in each case indicates the direction of divergence of the population mean for \code{x} (or difference of means for \code{x} and \code{y}) from \code{mu} (i.e., \code{"greater"}, \code{"less"}, or \code{"two.sided"}).
#'  
#' @section Test Assumptions:
#' The assumption of normality for the underlying distribution or a sufficiently large sample size is required along with the population standard deviation to use Z procedures.
#' 
#' @section Confidence Intervals:
#'  For each of the above tests, an expression for the related confidence interval (returned component \code{conf.int}) can be obtained in the usual way by inverting the expression for the test statistic. Note that, as explained under the description of \code{conf.int}, the confidence interval will be half-infinite when alternative is not \code{"two.sided"}; infinity will be represented by \code{Inf}.
#'  
#' @references \itemize{
#'  \item Kitchens, L.J. 2003. \emph{Basic Statistics and Data Analysis}. Duxbury. 
#'  \item Hogg, R. V. and Craig, A. T. 1970. \emph{Introduction to Mathematical Statistics, 3rd ed}. Toronto, Canada: Macmillan. 
#'  \item Mood, A. M., Graybill, F. A. and Boes, D. C. 1974. \emph{Introduction to the Theory of Statistics, 3rd ed}. New York: McGraw-Hill.
#'  \item Snedecor, G. W. and Cochran, W. G. 1980. \emph{Statistical Methods, 7th ed}. Ames, Iowa: Iowa State University Press. 
#'  }
#' 
#' @author Alan T. Arnholt <arnholtat@@appstate.edu> 
#'  
#' @seealso \code{\link{z.test}}, \code{\link{tsum.test}}
#' @export
#'  
#' @examples
#'
#' zsum.test(mean.x = 56/30,sigma.x = 2, n.x = 30, alternative="greater", mu = 1.8)
#' # Example 9.7 part a. from PASWR.
#' x <- rnorm(12)
#' zsum.test(mean(x), sigma.x = 1, n.x = 12)
#' # Two-sided one-sample z-test where the assumed value for 
#' # sigma.x is one. The null hypothesis is that the population 
#' # mean for 'x' is zero. The alternative hypothesis states 
#' # that it is either greater or less than zero. A confidence 
#' # interval for the population mean will be computed.
#' # Note: returns same answer as: 
#' z.test(x, sigma.x = 1)
#'
#' x <- c(7.8, 6.6, 6.5, 7.4, 7.3, 7.0, 6.4, 7.1, 6.7, 7.6, 6.8)
#' y <- c(4.5, 5.4, 6.1, 6.1, 5.4, 5.0, 4.1, 5.5)
#' zsum.test(mean(x), sigma.x = 0.5, n.x = 11 ,mean(y), sigma.y = 0.5, n.y = 8, mu = 2)
#' # Two-sided standard two-sample z-test where both sigma.x 
#' # and sigma.y are both assumed to equal 0.5. The null hypothesis 
#' # is that the population mean for 'x' less that for 'y' is 2. 
#' # The alternative hypothesis is that this difference is not 2. 
#' # A confidence interval for the true difference will be computed.
#' # Note: returns same answer as: 
#' z.test(x, sigma.x = 0.5, y, sigma.y = 0.5)
#' #
#' zsum.test(mean(x), sigma.x = 0.5, n.x = 11, mean(y), sigma.y = 0.5, n.y = 8, 
#' conf.level=0.90)
#' # Two-sided standard two-sample z-test where both sigma.x and 
#' # sigma.y are both assumed to equal 0.5. The null hypothesis 
#' # is that the population mean for 'x' less that for 'y' is zero.  
#' # The alternative hypothesis is that this difference is not
#' # zero.  A 90% confidence interval for the true difference will 
#' # be computed.  Note: returns same answer as:
#' z.test(x, sigma.x=0.5, y, sigma.y=0.5, conf.level=0.90)
#' rm(x, y) 
#'  
#' @keywords htest 
#'  
###########################################################################
zsum.test <- function (mean.x, sigma.x = NULL, n.x = NULL, mean.y = NULL, sigma.y = NULL, n.y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, conf.level = 0.95, ... ) 
{
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || 
                                 !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  ################################################################################
  if ( (!is.null(mean.y) || !is.null(sigma.y)  || !is.null(n.y) ) && 
         (is.null(sigma.x) || is.null(n.x) ) )
    stop("You must enter summarized values for both samples (x and y)")
  if (!is.null(mean.x) && is.null(mean.y) && is.null(n.x) && 
        is.null(sigma.x)) 
    stop("You must enter the value for both sigma.x and n.x")
  if (is.null(n.x) && !is.null(mean.x) && !is.null(sigma.x) && 
        is.null(mean.y)) 
    stop("You must enter the value for n.x")
  if (is.null(sigma.x) && !is.null(mean.x) && !is.null(n.x) && 
        is.null(mean.y)) 
    stop("You must enter the value for sigma.x")
  if (is.null(n.y) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(sigma.y) && !is.null(sigma.x) && !is.null(n.x)) 
    stop("You must enter the value for n.y")
  if (is.null(n.y) && is.null(n.x) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(sigma.y) && !is.null(sigma.x)) 
    stop("You must enter the value for both n.x and n.y")
  if (is.null(sigma.x) && is.null(sigma.y) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(n.x) && !is.null(n.y)) 
    stop("You must enter the value for both sigma.x and sigma.y")
  if (!is.null(sigma.x) && is.null(sigma.y) && !is.null(mean.x) && 
        !is.null(mean.y) && !is.null(n.x) && !is.null(n.y)) 
    stop("You must enter the value for sigma.y")
  if (is.null(n.y) && is.null(sigma.y) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(sigma.x) && !is.null(n.x)) 
    stop("You must enter the value for both sigma.y and n.y") 
  if ( (!is.null(mean.y) && !is.null(sigma.y) && !is.null(n.y)) && 
         (is.null(mean.x) || is.null(sigma.x) || is.null(n.x)) ) 
    stop("You must enter values for mean.x, sigma.x, and n.x")
  ################################################################################  
  if (!is.null(mean.y)) {
    dname <- "User input summarized values for x and y"
    yok <- !is.na(mean.y)
    xok <- !is.na(mean.x)
    mean.y <- mean.y[yok]
  }
  else {
    dname <- "User input summarized values for x"   
    xok <- !is.na(mean.x)
    yok <- NULL
  }
  mean.x <- mean.x[xok]
  nx <- n.x
  mx <- mean.x
  vx <- sigma.x^2
  if (is.null(mean.y)) {
    if (nx < 2) 
      stop("not enough 'x' observations")
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
    zstat <- (mx - mu)/stderr
    method <- "One Sample z-test"
    estimate <- setNames(mx, "mean of x")
  }
  else {
    ny <- n.y
    if (nx < 1) 
      stop("not enough 'x' observations")
    if (ny < 1) 
      stop("not enough 'y' observations")
    my <- mean.y
    vy <- sigma.y^2
    method <- "Two Sample z-test"
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    stderrx <- sqrt(vx/nx)
    stderry <- sqrt(vy/ny)
    stderr <- sqrt(stderrx^2 + stderry^2)
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
      stop("data are essentially constant")
    zstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pnorm(zstat)
    cint <- c(-Inf, zstat + qnorm(conf.level))
  }
  else if (alternative == "greater") {
    pval <- pnorm(zstat, lower.tail = FALSE)
    cint <- c(zstat - qnorm(conf.level), Inf)
  }
  else {
    pval <- 2 * pnorm(-abs(zstat))
    alpha <- 1 - conf.level
    cint <- qnorm(1 - alpha/2)
    cint <- zstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(zstat) <- "z"
  names(mu) <- "difference in means"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = zstat, p.value = pval, 
               conf.int = cint, estimate = estimate, null.value = mu, 
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}