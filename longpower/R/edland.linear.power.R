#' Linear mixed model sample size calculations.
#' 
#' This function performs the sample size calculation for a linear mixed model
#' with random slope.
#' 
#' This function will also provide sample size estimates for linear mixed
#' models with random intercept only simply by setting \code{sig2.s = 0}
#' 
#' @param n sample size per group
#' @param delta group difference in slopes
#' @param t the observation times
#' @param sig2.s variance of random slope
#' @param sig2.e residual variance
#' @param sig.level type one error
#' @param power power
#' @param alternative one- or two-sided test
#' @return The number of subject required per arm to attain the specified
#' \code{power} given \code{sig.level} and the other parameter estimates.
#' @author Michael C. Donohue, Steven D. Edland
#' @seealso \code{\link{lmmpower}}, \code{\link{diggle.linear.power}},
#' \code{\link{liu.liang.linear.power}}
#' @references Edland, S.D. (2009) Which MRI measure is best for Alzheimer's
#' disease prevention trials: Statistical considerations of power and sample
#' size. \emph{Joint Stat Meeting Proceedings}. 4996-4999.
#' @keywords power sample size mixed effects random effects
#' @examples
#' 
#' \dontrun{
#' browseVignettes(package = "longpower")
#' }
#' # Reproduces the table on page 29 of Diggle et al
#' n = 3
#' t = c(0,2,5)
#' rho = c(0.2, 0.5, 0.8)
#' sigma2 = c(100, 200, 300)
#' tab = outer(rho, sigma2, 
#'       Vectorize(function(rho, sigma2){
#'         ceiling(edland.linear.power(
#'           delta=0.5,
#'           t=t,
#'           sig2.e=sigma2*(1-rho),
#'           alternative="one.sided",
#'           power=0.80)$n)}))
#' colnames(tab) = paste("sigma2 =", sigma2)
#' rownames(tab) = paste("rho =", rho)
#' tab
#' 
#' # An Alzheimer's Disease example using ADAS-cog pilot estimates
#' t = seq(0,1.5,0.25)
#' n = length(t)
#' 
#' edland.linear.power(delta=1.5, t=t, sig2.s = 24, sig2.e = 10, sig.level=0.05, power = 0.80)
#' 
edland.linear.power <- function(n = NULL, delta = NULL, t = NULL, sig2.s = 0, sig2.e = 1, 
         sig.level=0.05, power=NULL,
         alternative = c("two.sided", "one.sided"))
{
  if (sum(sapply(list(n, delta, sig2.s, sig2.e, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'n', 'delta', 'sig2.s', 'sig2.e', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  n.body <- quote({
    2*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) + qnorm(1-power))^2 * 
      (sig2.s + sig2.e/sum((t - mean(t))^2)) / delta^2
  })
  
  if (is.null(sig.level)) 
      sig.level <- uniroot(function(sig.level) eval(n.body) - 
          n, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power)) 
      power <- uniroot(function(power) eval(n.body) - 
          n, c(1e-3, 1 - 1e-10))$root
  else if (is.null(delta)) 
      delta <- uniroot(function(delta) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  else if (is.null(sig2.s)) 
      sig2.s <- uniroot(function(sig2.s) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  else if (is.null(sig2.e)) 
      sig2.e <- uniroot(function(sig2.e) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  n <- eval(n.body)

  METHOD <- "Power for longitudinal linear model with random slope (Edland, 2009)"
  structure(list(n = n, delta = delta, sig2.s = sig2.s, sig2.e = sig2.e, sig.level = sig.level, 
        t=t, power = power, alternative = alternative, 
        note = "n is number in *each* group",
        method = METHOD), class = "power.longtest")
}
