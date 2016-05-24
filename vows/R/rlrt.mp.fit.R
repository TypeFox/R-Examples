#' Massively parallel restricted likelihood ratio tests (internal)
#' 
#' Conducts a possibly very large number of restricted likelihood ratio tests
#' (Crainiceanu and Ruppert, 2004), with specified random-effects design matrix
#' and fixed-effects design matrix, for a polynomial null against a smooth
#' alternative.
#' 
#' The \pkg{RLRsim} package of Scheipl et al. (2008) is used to simulate the
#' common null distribution of the RLRT statistics.
#' 
#' @param Y ordinarily, an \eqn{n \times V} outcome matrix, where \eqn{V} is
#' the number of hypotheses (in brain imaging applications, the number of
#' voxels
#' @param X the fixed-effects design matrix.
#' @param Z the random-effects design matrix.
#' @param loginvsp a grid of candidate values of the log inverse smoothing
#' parameter.
#' @param evalarg if \code{Y} is of class "fd", the argument values at which
#' the functions are evaluated.
#' @param get.df logical: Should the effective df of the smooth at each point
#' be obtained?
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
#' @references Crainiceanu, C. M., and Ruppert, D. (2004).  Likelihood ratio
#' tests in linear mixed models with one variance component.  \emph{Journal of
#' the Royal Statistical Society, Series B}, 66(1), 165--185.
#' 
#' Scheipl, F., Greven, S. and Kuechenhoff, H. (2008). Size and power of tests
#' for a zero random effect variance or polynomial regression in additive and
#' linear mixed models. \emph{Computational Statistics & Data Analysis}, 52(7),
#' 3283--3299.
#' @examples
#' 
#' Y = matrix(rnorm(6000), nrow=20)
#' x = rnorm(20)
#' z = rep(1:5, each = 4)
#' t4. = rlrt.mp.fit(Y, x, z, loginvsp = -22:0)
#' @export
rlrt.mp.fit <-
function(Y, X, Z, loginvsp, evalarg=NULL, get.df=FALSE) {
    sim = RLRsim::RLRTSim(X, Z, qrX=qr(X), sqrt.Sigma=diag(NCOL(Z)))
	if (!exists("is.fd")) is.fd = function(mm) FALSE
	n = NROW(X); p = NCOL(X); q = NCOL(Z)
	svdZZ = svd(tcrossprod(Z))
	d = svdZZ$d; d[-(1:q)] = 0
	X. = crossprod(svdZZ$u, X)
	X.X. = crossprod(X.); detX.X. = det(X.X.)
	Y. = crossprod(svdZZ$u, if (is.fd(Y)) t(Y$coef) else Y)
	if (is.fd(Y)) {
		if (is.null(evalarg)) stop("'evalarg' must be specified when 'Y' is of class 'fd'")
		B = eval.basis(evalarg, Y$basis) 
	}
	m1 = Y. - X. %*% solve(X.X., crossprod(X., Y.))
	if (is.fd(Y)) term1 = (n-p) * log(colSums(tcrossprod(Y.,B) * tcrossprod(m1,B)))
	else term1 = (n-p) * log(colSums(Y. * m1))
	rldiff2 = Vectorize(function(log.inv.sp) {
		ev.v = 1 + exp(log.inv.sp)*d
		di = diag(1/sqrt(ev.v))
		X.. = di %*% X.
		X..X.. = crossprod(X..)
		Y.. = di %*% Y.
		m2 = Y.. - X.. %*% solve(X..X.., crossprod(X.., Y..))
		if (is.fd(Y)) term2 = (n-p) * log(colSums(tcrossprod(Y..,B) * tcrossprod(m2,B)))
	    else term2 = (n-p) * log(colSums(Y.. * m2))
		term1 - term2 - sum(log(ev.v)) - log(det(X..X..)) + log(detX.X.)
	})
	if (get.df) {  # see RWC, p. 336
		Rinv = solve(chol(crossprod(cbind(X, Z))))
		Rinv22 = Rinv[-(1:p), -(1:p)]
		svec = c(svd(Rinv22)$d^2, rep(0, p))
	}
	tabl = matrix(rldiff2(loginvsp), ncol = length(loginvsp))
	colnames(tabl) = loginvsp
	logsp = -loginvsp[apply(tabl, 1, which.max)]
    stat = apply(tabl, 1, max)
    pvalue = stat
    for (i in 1:length(stat)) pvalue[i] = mean(stat[i] < sim)
    fdr = p.adjust(pvalue, method="BH")
	list(table = tabl, stat = stat, logsp = logsp, df = if (get.df) rowSums(1/(1+exp(logsp) %o% svec)) else NULL, sim = sim, pvalue=pvalue, fdr=fdr)
}

