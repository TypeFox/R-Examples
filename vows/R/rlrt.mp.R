#' Massively parallel restricted likelihood ratio tests
#' 
#' Conducts a possibly very large number of restricted likelihood ratio tests
#' (Crainiceanu and Ruppert, 2004), with common design matrix, for a polynomial
#' null against a smooth alternative.
#' 
#' The \pkg{RLRsim} package of Scheipl et al. (2008) is used to simulate the
#' common null distribution of the RLRT statistics.
#' 
#' @param Y ordinarily, an \eqn{n \times V} outcome matrix, where \eqn{V} is
#' the number of hypotheses (in brain imaging applications, the number of
#' voxels).  Can also be given by an object of class "\code{\link[fda]{fd}}".
#' @param x a vector or matrix of covariates.
#' @param nbasis number of B-spline basis functions.
#' @param norder order of B-splines.
#' @param nulldim dimension of the null space of the penalty.
#' @param loginvsp a grid of candidate values of the log inverse smoothing
#' parameter.
#' @param evalarg if \code{Y} is of class "\code{\link[fda]{fd}}", the argument
#' values at which the functions are evaluated.
#' @param get.df logical: Should the effective df of the smooth at each point
#' be obtained?
#' @param B evaluation matrix of the B-spline basis functions.
#' @param P penalty matrix.
#' @return A list with components \item{table}{matrix of log restricted
#' likelihood ratio values at each grid point, for each test.} \item{stat}{RLRT
#' statistics, i.e., the supremum of the values in \code{table} for each test.}
#' \item{logsp}{log smoothing parameter at which the supremum of the restricted
#' likelihood ratio is attained for each test.} \item{df}{if \code{get.df =
#' TRUE}, the effective degrees of freedom corresponding to the log smoothing
#' parameter values in \code{logsp}.} \item{sim}{values simulated from the null
#' distribution of the restricted likelihood ratio statistic.}
#' \item{pvalue}{p-values for the RLRT statistics.}
#' \item{fdr}{Benjamini-Hochberg false discovery rates corresponding to the
#' above p-values.} \item{call}{the call to the function.}
#' @author Lei Huang \email{huangracer@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{rlrt4d}}, and \code{\link{Fdr.rlrt}} for a more
#' sophisticated false discovery rate procedure.
#' @references Crainiceanu, C. M., and Ruppert, D. (2004).  Likelihood ratio
#' tests in linear mixed models with one variance component.  \emph{Journal of
#' the Royal Statistical Society, Series B}, 66(1), 165--185.
#' 
#' Reiss, P. T., Huang, L., Chen, Y.-H., Huo, L., Tarpey, T., and Mennes, M.
#' (2014). Massively parallel nonparametric regression, with an application to
#' developmental brain mapping. \emph{Journal of Computational and Graphical
#' Statistics}, \emph{Journal of Computational and Graphical Statistics},
#' 23(1), 232--248.
#' 
#' Scheipl, F., Greven, S. and Kuechenhoff, H. (2008). Size and power of tests
#' for a zero random effect variance or polynomial regression in additive and
#' linear mixed models. \emph{Computational Statistics & Data Analysis}, 52(7),
#' 3283--3299.
#' @examples
#' 
#' Y = matrix(rnorm(6000), nrow=20)
#' x = rnorm(20)
#' t4 = rlrt.mp(Y, x, loginvsp = -22:0)
#' f4 = Fdr.rlrt(t4, 6)
#' @export
rlrt.mp <-
function(Y, x=NULL, loginvsp, nbasis = 15, norder=4, nulldim=NULL, evalarg=NULL, get.df=FALSE, B=NULL, P=NULL) {
    if (is.null(x) & is.null(B))   stop("You must enter either x or B")
    if (is.null(x) & is.null(P))   stop("You must enter P")
    
    n = if (!is.null(x)) NROW(x) else NROW(B)
  	
    xz = rlr.xz(x, nbasis=nbasis, norder=norder, nulldim=nulldim, B=B, P=P)
    rlobj = rlrt.mp.fit(Y=Y, X=xz$X, Z=xz$Z, loginvsp=loginvsp, evalarg=evalarg, get.df=get.df)
    rlobj
}

