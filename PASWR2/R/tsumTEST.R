#' @title Summarized t-test
#' 
#' @description Performs a one-sample, two-sample, or a Welch modified two-sample t-test based on user supplied summary information. Output is identical to that produced with \code{t.test}.
#' 
#' @details If \code{y} is \code{NULL}, a one-sample t-test is carried out with \code{x}. If \code{y} is not \code{NULL}, either a standard or Welch modified two-sample t-test is performed, depending on whether \code{var.equal} is \code{TRUE} or \code{FALSE}.
#' 
#' @param mean.x a single number representing the sample mean of \code{x}
#' @param s.x a single number representing the sample standard deviation of \code{x}
#' @param n.x a single number representing the sample size of \code{x}
#' @param mean.y a single number representing the sample mean of \code{y}
#' @param s.y a single number representing the sample standard deviation of \code{y}
#' @param n.y a single number representing the sample size of \code{y}
#' @param alternative is a character string, one of \code{"greater"}, \code{"less"}, or \code{"two.sided"}, or just the initial letter of each, indicating the specification of the alternative hypothesis. For one-sample tests, \code{alternative} refers to the true mean of the parent population in relation to the hypothesized value \code{mu}.  For the standard two-sample tests, \code{alternative} refers to the difference between the true population mean for \code{x} and that for \code{y}, in relation to \code{mu}. For the one-sample and paired t-tests, \code{alternative} refers to the true mean of the parent population in relation to the hypothesized value \code{mu}. For the standard and Welch modified two-sample t-tests, \code{alternative} refers to the difference between the true population mean for \code{x} and that for \code{y}, in relation to \code{mu}.  For the one-sample t-tests, alternative refers to the true mean of the parent population in relation to the hypothesized value \code{mu}. For the standard and Welch modified two-sample t-tests, alternative refers to the difference between the true population mean for \code{x} and that for \code{y}, in relation to \code{mu}.
#' @param mu is a single number representing the value of the mean or difference in means specified by the null hypothesis.
#' @param var.equal logical flag: if \code{TRUE}, the variances of the parent populations of \code{x} and \code{y} are assumed equal. Argument \code{var.equal} should be supplied only for the two-sample tests.
#' @param conf.level is the confidence level for the returned confidence interval; it must lie between zero and one.
#' @param ... Other arguments passed onto \code{tsum.test()}
#'
#' @return A list of class \code{htest}, containing the following components:
#' \item{\code{statistic}}{the t-statistic, with names attribute \code{"t"}}
#' \item{\code{parameters}}{is the degrees of freedom of the t-distribution associated with statistic. Component \code{parameters} has names attribute \code{"df"}.}
#' \item{\code{p.value}}{the p-value for the test}
#' \item{\code{conf.int}}{is a confidence interval (vector of length 2) for the true mean or difference in means. The confidence level is recorded in the attribute \code{conf.level}. When alternative is not \code{"two.sided"}, the confidence interval will be half-infinite, to reflect the interpretation of a confidence interval as the set of all values \code{k} for which one would not reject the null hypothesis that the true mean or difference in means is \code{k} . Here infinity will be represented by \code{Inf}.}
#' \item{\code{estimate}}{is a vector of length 1 or 2, giving the sample mean(s) or mean of differences; these estimate the corresponding population parameters. Component \code{estimate} has a names attribute describing its elements.} 
#' \item{\code{null.value}}{is the value of the mean or difference in means specified by the null hypothesis. This equals the input argument \code{mu}. Component \code{null.value} has a names attribute describing its elements.}
#' \item{alternative}{records the value of the input argument alternative: \code{"greater"} , \code{"less"} or \code{"two.sided"}.}
#' \item{data.name}{is a character string (vector of length 1) containing the names x and y for the two summarized samples.}
#' 
#' @section Null Hypothesis:
#'  For the one-sample t-test, the null hypothesis is that the mean of the population from which \code{x} is drawn is \code{mu}. For the standard and Welch modified two-sample t-tests, the null hypothesis is that the population mean for \code{x} less that for \code{y} is \code{mu}.
#'   
#'  The alternative hypothesis in each case indicates the direction of divergence of the population mean for \code{x} (or difference of means for \code{x} and \code{y}) from \code{mu} (i.e., \code{"greater"}, \code{"less"}, or \code{"two.sided"}).
#' 
#' @section Test Assumptions:
#'  The assumption of equal population variances is central to the standard two-sample t-test. This test can be misleading when population variances are not equal, as the null distribution of the test statistic is no longer a t-distribution. If the assumption of equal variances is doubtful with respect to a particular dataset, the Welch modification of the t-test should be used. 
#'  
#'   The t-test and the associated confidence interval are quite robust with respect to level toward heavy-tailed non-Gaussian distributions (e.g., data with outliers). However, the t-test is non-robust with respect to power, and the confidence interval is non-robust with respect to average length, toward these same types of distributions.
#'   
#' @section Confidence Intervals:
#'  For each of the above tests, an expression for the related confidence interval (returned component \code{conf.int}) can be obtained in the usual way by inverting the expression for the test statistic. Note that, as explained under the description of \code{conf.int}, the confidence interval will be half-infinite when alternative is not \code{"two.sided"} ; infinity will be represented by \code{Inf}.
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
#' @seealso \code{\link{z.test}}, \code{\link{zsum.test}}
#' @export
#'  
#' @examples
#'  
#' # 95% Confidence Interval for mu1 - mu2, assuming equal variances
#' round(tsum.test(mean.x = 53/15, mean.y = 77/11, s.x=sqrt((222 - 15*(53/15)^2)/14),
#' s.y = sqrt((560 - 11*(77/11)^2)/10), n.x = 15, n.y = 11, var.equal = TRUE)$conf, 2)
#' # One Sample t-test
#' tsum.test(mean.x = 4, s.x = 2.89, n.x = 25, mu = 2.5)  
#'  
#' @keywords htest                                          
##########################  Updated 6/11/13  ###################################
tsum.test <- function (mean.x, s.x = NULL, n.x = NULL, mean.y = NULL, s.y = NULL, n.y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, var.equal = FALSE, conf.level = 0.95, ... ) 
{
  alternative <- match.arg(alternative)
  if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
    stop("'mu' must be a single number")
  if (!missing(conf.level) && (length(conf.level) != 1 || 
                                 !is.finite(conf.level) || 
                                 conf.level < 0 || conf.level > 1)) 
    stop("'conf.level' must be a single number between 0 and 1")
  ################################################################################
  if ( (!is.null(mean.y) || !is.null(s.y)  || !is.null(n.y) ) && 
        ( is.null(s.x) || is.null(n.x) ) )
    stop("You must enter summarized values for both samples (x and y)")
  if (!is.null(mean.x) && is.null(mean.y) && is.null(n.x) && 
        is.null(s.x)) 
    stop("You must enter the value for both s.x and n.x")
  if (is.null(n.x) && !is.null(mean.x) && !is.null(s.x) && 
        is.null(mean.y)) 
    stop("You must enter the value for n.x")
  if (is.null(s.x) && !is.null(mean.x) && !is.null(n.x) && 
        is.null(mean.y)) 
    stop("You must enter the value for s.x")
  if (is.null(n.y) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(s.y) && !is.null(s.x) && !is.null(n.x)) 
    stop("You must enter the value for n.y")
  if (is.null(n.y) && is.null(n.x) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(s.y) && !is.null(s.x)) 
    stop("You must enter the value for both n.x and n.y")
  if (is.null(s.x) && is.null(s.y) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(n.x) && !is.null(n.y)) 
    stop("You must enter the value for both s.x and s.y")
  if (!is.null(s.x) && is.null(s.y) && !is.null(mean.x) && 
        !is.null(mean.y) && !is.null(n.x) && !is.null(n.y)) 
    stop("You must enter the value for s.y")
  if (is.null(n.y) && is.null(s.y) && !is.null(mean.x) && !is.null(mean.y) && 
        !is.null(s.x) && !is.null(n.x)) 
    stop("You must enter the value for both s.y and n.y") 
  if ( (!is.null(mean.y) && !is.null(s.y) && !is.null(n.y)) && 
         (is.null(mean.x) || is.null(s.x) || is.null(n.x)) ) 
    stop("You must enter values for mean.x, s.x, and n.x")
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
  vx <- s.x^2
  if (is.null(mean.y)) {
    if (nx < 2) 
      stop("not enough 'x' observations")
    df <- nx - 1
    stderr <- sqrt(vx/nx)
    if (stderr < 10 * .Machine$double.eps * abs(mx)) 
      stop("data are essentially constant")
    tstat <- (mx - mu)/stderr
    method <- "One Sample t-test"
    estimate <- setNames(mx, "mean of x")
  }
  else {
    ny <- n.y
    if (nx < 1 || (!var.equal && nx < 2)) 
      stop("not enough 'x' observations")
    if (ny < 1 || (!var.equal && ny < 2)) 
      stop("not enough 'y' observations")
    if (var.equal && nx + ny < 3) 
      stop("not enough observations")
    my <- mean.y
    vy <- s.y^2
    method <- paste(if (!var.equal) 
      "Welch", "Two Sample t-test")
    estimate <- c(mx, my)
    names(estimate) <- c("mean of x", "mean of y")
    if (var.equal) {
      df <- nx + ny - 2
      v <- 0
      if (nx > 1) 
        v <- v + (nx - 1) * vx
      if (ny > 1) 
        v <- v + (ny - 1) * vy
      v <- v/df
      stderr <- sqrt(v * (1/nx + 1/ny))
    }
    else {
      stderrx <- sqrt(vx/nx)
      stderry <- sqrt(vy/ny)
      stderr <- sqrt(stderrx^2 + stderry^2)
      df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1))
    }
    if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) 
      stop("data are essentially constant")
    tstat <- (mx - my - mu)/stderr
  }
  if (alternative == "less") {
    pval <- pt(tstat, df)
    cint <- c(-Inf, tstat + qt(conf.level, df))
  }
  else if (alternative == "greater") {
    pval <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df), Inf)
  }
  else {
    pval <- 2 * pt(-abs(tstat), df)
    alpha <- 1 - conf.level
    cint <- qt(1 - alpha/2, df)
    cint <- tstat + c(-cint, cint)
  }
  cint <- mu + cint * stderr
  names(tstat) <- "t"
  names(df) <- "df"
  names(mu) <- "difference in means"
  attr(cint, "conf.level") <- conf.level
  rval <- list(statistic = tstat, parameter = df, p.value = pval, 
               conf.int = cint, estimate = estimate, null.value = mu, 
               alternative = alternative, method = method, data.name = dname)
  class(rval) <- "htest"
  return(rval)
}
