#' Nonparametric Test of Symmetry Using the Periodogram
#'
#' This function performs the nonparametric tests of reflection and complete symmetry using the periodogram from Lu and Zimmerman (2005) for spatial data with sampling locations on the integer grid. See Lu and Zimmerman (2005) for more details.
#'
#' @export
#' @keywords external
#'
#' @param spdata An \eqn{n \times 3} matrix. The first two columns provide \eqn{(x,y)} spatial coordinates. The third column provides data values at the coordinates.
#'
#' @param nrows	The number of rows of observed data.
#' @param ncols	The number of columns of observed data.
#' @param test	A string taking the value \code{reflection} or \code{complete} for a test of reflection or complete symmetry. If \code{test = "complete"}, a test for complete symmetry is performed after a test of reflection symmetry is performed and both p-values are returned.
#' @param nsim	The number simulations used to approximate the sampling distribution of CvM and CvM* from Lu ad Zimmerman (2005).
#'
#' @details The function assumes data are on the integer grid, \eqn{Z^2}. It uses the (unsmoothed) periodogram, the Fourier transform of the sample covariance function, to test symmetry properties.
#'
#' @return \item{pvalue.refl}{The p-value for the test of reflection symmetry computed by the CvM GoF test.}
#'  \item{pvalue.comp}{If \code{test = "complete"}, the p-value for the test of complete symmetry computed by using the CvM* GoF test.}
#'
#' @references Lu, N., & Zimmerman, D. L. (2005). Testing for directional symmetry in spatial dependence using the periodogram. \emph{Journal of Statistical Planning and Inference}, 129(1), 369-385.
#'
#' @seealso \code{\link{GuanTestGrid}}
#' @examples
#' library(mvtnorm)
#' set.seed(1)
#' #Number of rows and columns
#' nr <- 15
#' nc <- 15
#' n <- nr*nc
#' #Set up the coordinates
#' coords <- expand.grid(0:(nr-1), 0:(nc-1))
#' coords <- cbind(coords[,2], coords[,1])
#' #Compute the distance between sampling locations
#' D <- as.matrix(dist(coords))
#' #Set parameter values for exponential covariance function
#' sigma.sq <- 1
#' tau.sq <- 0.0
#' phi <- 1/4
#' #Simulate data using isotropic covariance function
#' D <- as.matrix(dist(coords))
#' R <- sigma.sq * exp(-phi*D)
#' R <- R + diag(tau.sq, nrow = n, ncol = n)
#' z <- rmvnorm(1,rep(0,n), R, method = c("chol"))
#' z <-  z-mean(z)
#' z <- t(z)
#' mydata <- cbind(coords, z)
#' #Run the test on the data from an isotropic (symmetric) covariance function
#' tr <- LuTest(mydata, nr, nc, test = "complete", nsim = 1000)
#' tr
#'
#' #Simulate data from anisotropic (non-symmetric) covariance function
#' aniso.angle <- pi/4
#' aniso.ratio <- 2
#' coordsA <- coords.aniso(coords, c(aniso.angle, aniso.ratio))
#' Da <- as.matrix(dist(coordsA))
#' R <- sigma.sq * exp(-phi*Da)
#' R <- R + diag(tau.sq, nrow = n, ncol = n)
#' z <- rmvnorm(1,rep(0,n), R, method = c("chol"))
#' z <-  z-mean(z)
#' z <- t(z)
#' mydata <- cbind(coords, z)
#' #Run the test on the data generated from an anisotropic 
#' #(and non reflection and non completely symmetric) covariance function
#' tr <- LuTest(mydata, nr, nc, test = "complete", nsim = 1000)
#' tr
LuTest = function(spdata, nrows, ncols, test = "complete", nsim = 5000)
{
	dname <- deparse(substitute(spdata))

	if(!is.matrix(spdata))
	{stop("spdata must be a matrix")}
	if(dim(spdata)[2] != 3)
	{stop("matrix of spatial data must have 3 columns")}
	if(dim(spdata)[1] < 10)
	{stop("matrix of spatial data must have at least 10 rows")}
	if(nrows <= 1 )
	{stop("nrows must be greater than 1")}
	if(ncols <= 1 )
	{stop("ncols must be greater than 1")}
	if(test != "reflection" & test != "complete")
	{stop("'test' must be either 'reflection' or 'complete'")}
	
	lags <- cov_lags(nrows, ncols)
	spdata <- scale_coords(spdata)
	chat <- est_cov(lags, spdata, nrows, ncols)
	cvec <- cov_complete(chat, lags, ret.mat = F)
	ffs <- get_Fourier_freqs(nrows, ncols)
	pe <- periodogram(cvec, ffs)
	pvalue.refl <- test_reflection_sym(pe,nrows,ncols, nsim)
	
	if(test == "reflection")
	{
		rv <- list("pvalue.refl" = pvalue.refl, pvalue.comp = NULL)
		return(rv)
	}

	if(test == "complete")
	{
		pvalue.comp <- test_complete_sym(pe, nsim)
		rv <- list("pvalue.refl" = pvalue.refl, "pvalue.comp" = pvalue.comp)
	}
	
	htestIso.LZ <- make_htestIso_LZ(rv, df = NULL)
	htestIso.LZ$data.name <- dname
	return(htestIso.LZ)
}


