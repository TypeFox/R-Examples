#' Linear mixed model sample size calculations.
#' 
#' This function performs the sample size calculation for a linear mixed model.
#' See Diggle et al (2002) for parameter definitions and other details.
#' 
#' The parameters \code{u}, \code{v}, and \code{Pi} are expected to be the same
#' length and sorted with respect to each other. See Diggle, et al (1997) and
#' package vignette for more details.
#' 
#' @param n sample size per group
#' @param delta group difference in slopes
#' @param t the observation times
#' @param sigma2 the marginal model (GEE) scale parameter
#' @param R the working correlation matrix (or variance-covariance matrix if
#' \code{sigma2} is 1). If \code{R} is a scalar, an exchangeable working
#' correlation matrix will be assumed.
#' @param sig.level Type I error
#' @param power power
#' @param alternative one- or two-sided test
#' @return The number of subject required per arm to attain the specified
#' \code{power} given \code{sig.level} and the other parameter estimates.
#' @author Michael C. Donohue, Steven D. Edland
#' @seealso \code{\link{lmmpower}}, \code{\link{diggle.linear.power}}
#' @references Diggle P.J., Heagerty P.J., Liang K., Zeger S.L. (2002)
#' \emph{Analysis of longitudinal data}. Second Edition. Oxford Statistical
#' Science Series.
#' @keywords power sample size mixed effects random effects
#' @examples
#' 
#' \dontrun{
#' browseVignettes(package = "longpower")
#' }
#' 
#' # Reproduces the table on page 29 of Diggle et al
#' n = 3
#' t = c(0,2,5)
#' rho = c(0.2, 0.5, 0.8)
#' sigma2 = c(100, 200, 300)
#' tab = outer(rho, sigma2, 
#'       Vectorize(function(rho, sigma2){
#'         ceiling(diggle.linear.power(
#'           delta=0.5,
#'           t=t,
#'           sigma2=sigma2,
#'           R=rho,
#'           alternative="one.sided",
#'           power = 0.80)$n)}))
#' colnames(tab) = paste("sigma2 =", sigma2)
#' rownames(tab) = paste("rho =", rho)
#' tab
#' 
#' # An Alzheimer's Disease example using ADAS-cog pilot estimates
#' # var of random intercept
#' sig2.i = 55
#' # var of random slope
#' sig2.s = 24
#' # residual var
#' sig2.e = 10
#' # covariance of slope and intercep
#' cov.s.i <- 0.8*sqrt(sig2.i)*sqrt(sig2.s)
#' 
#' cov.t <- function(t1, t2, sig2.i, sig2.s, cov.s.i){
#'         sig2.i + t1*t2*sig2.s + (t1+t2)*cov.s.i 
#' }
#' 
#' t = seq(0,1.5,0.25)
#' n = length(t)
#' R = outer(t, t, function(x,y){cov.t(x,y, sig2.i, sig2.s, cov.s.i)})
#' R = R + diag(sig2.e, n, n)
#' 
#' diggle.linear.power(d=1.5, t=t, R=R, sig.level=0.05, power=0.80)
#' 
#' @export diggle.linear.power
diggle.linear.power <-
function(n=NULL, delta=NULL, t=NULL, sigma2=1, R=NULL, 
         sig.level=0.05, power=NULL,
         alternative=c("two.sided", "one.sided"))
{
  if (sum(sapply(list(n, delta, sigma2, power, sig.level), is.null)) != 1) 
      stop("exactly one of 'delta', 'sigma2', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
      sig.level | sig.level > 1)) 
      stop("'sig.level' must be numeric in [0, 1]")
  alternative <- match.arg(alternative)

  if(is.null(dim(R))) R = matrix(R, length(t), length(t)) + diag(1-R,length(t))
  
  n.body <- quote({
    V = sigma2*R
    xi = solve(rbind(1,t)%*%solve(V)%*%cbind(1,t))[2,2]
    2*(qnorm(ifelse(alternative=="two.sided", sig.level/2, sig.level)) +
       qnorm(1-power))^2*xi/delta^2
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
  else if (is.null(sigma2)) 
      sigma2 <- uniroot(function(sigma2) eval(n.body) - 
          n, c(1e-10, 1e5))$root
  n <- eval(n.body)

  METHOD <- "Longitudinal linear model slope power calculation (Diggle et al 2002, page 29)"
  structure(list(n = n, delta = delta, sigma2 = sigma2, R = R, sig.level = sig.level, 
        power = power, alternative = alternative, 
        note = "n is number in *each* group",
        method = METHOD), class = "power.longtest")
}
