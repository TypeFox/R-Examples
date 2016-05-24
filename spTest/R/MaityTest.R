#' Nonparametric Test of Isotropy Using the Sample Covariogram
#'
#' This function performs the nonparametric test of isotropy using the sample covariogram from Maity and Sherman (2012) for spatial data with sampling locations following any general spatial sampling design. It uses the Epanechnikov kernel function with an empirical bandwidth parameter to smooth over spatial lags. The asymptotic variance-covariance matrix is estimated using the grid based block bootstrap (GBBB) from Lahiri and Zhu (2006). See Maity and Sherman (2012) for more details.
#'
#' @export
#' @keywords external
#'
#' @param spdata	An \eqn{n \times 3} matrix. The first two columns provide \eqn{(x,y)} spatial coordinates. The third column provides data values at the coordinates.
#' @param lagmat A \eqn{k \times 2} matrix of spatial lags. Each row corresponds to a lag of the form \eqn{(x.lag, y.lag)} for which the covariogram value will be estimated.
#' @param A	A \eqn{d \times k} contrast matrix. The contrasts correspond to contrasts of the estimated covariogram at the lags given in 'lagmat'.
#' @param df A scalar indicating the row rank of the matrix \code{A}. This value gives the degrees of freedom for the asymptotic \eqn{\chi^2} distribution used to compute the p-value.
#' @param xlims A vector of length two providing the lower and upper x-limits of the sampling region. To ensure all sampling locations are included in the subsampling procedure, the x-limits should be \strong{slightly} wider than than the minimum and maximum observed x-coordinates of sampling locations.
#' @param ylims A vector of length two providing the lower and upper y-limits of the sampling region. To ensure all sampling locations are included in the subsampling procedure, the y-limits should be \strong{slightly} wider than than the minimum and maximum observed y-coordinates of sampling locations.
#' @param grid.spacing A vector of length two providing the x (width) and y (height) spacing, respectively, of the underlying grid laid on the sampling region to create spatial blocks. If the grid spacing width does not evenly divide the width of the sampling region, some data will be ommited during subsampling, i.e., the function does not handle partial windows. Same applies to grid spacing height and height of sampling region. See details for an example.
#' @param block.dims A vector of length two corresponding to the width and height of the blocks used in the GBBB. The width and height are given in terms of the spacing of the grid laid on the sampling region. See details for an example.
#' @param nBoot A scalar indicating the number of grid based block bootstrap (GBBB) samples to compute (see Lahiri and Zhu (2006) for details on GBBB). Recommended to be at least 50.
#' @param kappa A scalar corresponding to the tuning parameter to adjust empirical bandwidth.
#' @param user.bandwidth Logical. Defaults to \code{FALSE}. Set to \code{TRUE} to manually override the empirical bandwidth.
#' @param bandwidth When \code{user.bandwidth = TRUE}, a vector of length two providing the user-definted bandwidths for smoothing over spatial lags in the x and y directions.
#'
#' @details This function currently only supports square and rectangular sampling regions and does not support partial blocks. For example, suppose the sampling region runs from 0 to 20 in the x-direction (\code{xlims = c(0, 20)}) and from 0 to 30 in the y-direction (\code{ylims = c(0, 30)}) and an underlying grid of 1 by 1 is laid over the sampling region. Then an ideal value of \code{block.dims} would be (2,3) since  this would imply the blocks created have dimension 2 by 3 and these values evenly divide the width (20) and height (30), respectively, of the sampling region. Using the vector (3, 4.5) would imply that some data will not be used in the GBBB since these values would create partial blocks in the sampling region. 
#'
#'The value block.dims provides the width and height of the blocks in terms of the underlying grid laid on the sampling region. For example, if a grid with dimensions of grid.spacing = c(0.1, 0.2) is laid on the sampling region and block.dims = c(2,3) then the dimensions of the blocks created by the moving windows are (0.2, 0.6). The block dimensions of 0.2 and 0.6 should evenly divide the width and height, respectively, of the entire sampling region. Thus, the user must take care to ensure that the values of \code{grid.spacing} and \code{block.dims} are compatible with the dimensions of the sampling region. The easiest way to meet this constrain is to make the \code{grid.spacing} values a function of the \code{xlims} and \code{ylims} values. For example, to put down a \eqn{10 \times 10} grid on the domain, use \code{grid.spacing = (xlims[2]-xlims[1], ylim[2]-ylims[1])/10}. Then, setting \code{block.dims = c(2,2)} ensures that no data will be omitted during thes subsampling.
#'
#'To preserve the spatial dependence structure, the spatial blocks should have the same shape (i.e. square or rectangle) and orientation as the entire sampling domain.
#'
#' @return \item{C.hat}{A matrix of the spatial lags provided and the covariogram, \eqn{\hat{C}()}, point estimates at those lags used to construct the test statistic.}
#' \item{V.hat}{The estimate of the asymptotic variance-covariance matrix, \eqn{\hat{V}}, used to construct the test statistic.}
#' \item{n.boot}{The number of bootstrap samples used to estimate \eqn{V}.}
#' \item{test.stat}{The calculated test statistic.}
#' \item{pvalue.chisq}{The p-value computed using the asymptotic \eqn{\chi^2} distribution.}
#'
#' @references Maity, A., & Sherman, M. (2012). Testing for spatial isotropy under general designs. \emph{Journal of Statistical Planning and Inference}, 142(5), 1081-1091.
#' @references Lahiri, S. N., & Zhu, J. (2006). Resampling methods for spatial regression models under a class of stochastic designs. \emph{The Annals of Statistics}, 34(4), 1774-1813.
#'
#' @seealso \code{\link{GuanTestUnif}}
#'
#' @examples
#' library(mvtnorm)
#' set.seed(1)
#' #Sample Size
#' N <- 300
#' #Set parameter values for exponential covariance function
#' sigma.sq <- 1
#' tau.sq <- 0.0
#' phi <- 1/4
#' #Generate sampling locations
#' coords <-  cbind(runif(N,0,16), runif(N,0,16))
#' D <-  as.matrix(dist(coords))
#' R <- sigma.sq * exp(-phi*D)
#' R <- R + diag(tau.sq, nrow = N, ncol = N)
#' #Simulate Gaussian spatial data
#' z <- rmvnorm(1,rep(0,N), R, method = "chol")
#' z <- z - mean(z)
#' z <- t(z)
#' mydata <- cbind(coords, z)
#' mylags <- rbind(c(1,0), c(0, 1), c(1, 1), c(-1,1))
#' myA <- rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
#' my.xlims <- c(0, 16)
#' my.ylims <- c(0, 16)
#' my.grid <- c(1,1)
#' my.blockdims <- c(4,4)
#' #Number of bootstraps for demonstration only.
#' #In practice, use nBoot > 50
#' my.nBoot <- 3
#' tr <- MaityTest(mydata, mylags, myA, df = 2, my.xlims, my.ylims, my.grid, 
#' my.blockdims, nBoot = my.nBoot)
#' tr
#'
#' ####NOT RUN####
#' # #Simulate data from anisotropic covariance function
#' # aniso.angle <- pi/4
#' # aniso.ratio <- 2
#' # coordsA <- coords.aniso(coords, c(aniso.angle, aniso.ratio))
#' # Da <- as.matrix(dist(coordsA))
#' # R <- sigma.sq * exp(-phi*Da)
#' # R <- R + diag(tau.sq, nrow = N, ncol = N)
#' # z <- rmvnorm(1,rep(0,N), R, method = c("chol"))
#' # z <-  z-mean(z)
#' # z <- t(z)
#' # mydata <- cbind(coords, z)
#' # #Run the test on the data generated from an anisotropic covariance function
#' # #Increase the number of bootstraps
#' # my.nBoot = 100
#' # tr <- MaityTest(mydata, mylags, myA, df = 2, my.xlims, my.ylims, 
#' # my.grid, my.blockdims, nBoot = my.nBoot)
#' # tr
MaityTest = function(spdata, lagmat, A, df, xlims, ylims, grid.spacing = c(1,1), block.dims, nBoot = 100, kappa = 1, user.bandwidth = F, bandwidth = c(1,1))
{
	dname <- deparse(substitute(spdata))
	
	if(!is.matrix(spdata))
	{stop("spdata must be a matrix")}
	if(dim(spdata)[2] != 3)
	{stop("matrix of spatial data must have 3 columns")}
	if(dim(spdata)[1] <= 3)
	{stop("matrix of spatial data must have at least 4 rows")}
	if(dim(lagmat)[2] != 2)
	{stop("matrix of spatial lags must have 2 columns")}
	if(dim(lagmat)[1] != dim(A)[2])
	{stop("non-conformable `A` and `lagmat` matrix")}
	if(df <= 0)
	{stop("df must be greater than 0")}
	if(xlims[1] >= xlims[2])
	{stop("invalid x limits of sampling region")}
	if(ylims[1] >= ylims[2])
	{stop("invalid y limits of sampling region")}
	if(length(block.dims) != 2)
	{stop("block.dims must be length 2")}
	if(block.dims[1] <= 0 | block.dims[2] <= 0)
	{stop("block.dims must have positive entries")}
	if( (block.dims[1] %% 1) != 0 | (block.dims[2] %% 1) != 0)
	{warning("block.dims should be integer values")}
	block.w <- grid.spacing[1]*block.dims[1]
	block.h <- grid.spacing[2]*block.dims[2]
	if( (xlims[2]-xlims[1]) <= block.w)
	{stop("block width must be less than the width of sampling region, check grid.spacing[1] and block.dims[1]")}
	if( (ylims[2]-ylims[1]) <= block.h)
	{stop("block height must be less than the height of sampling region, check grid.spacing[2] and block.dims[2]")}
	if( ((xlims[2]-xlims[1])%%block.w) != 0 )
	{warning("width of blocks does not divide width of sampling region evenly, some data will be ommited during subsampling (function does not handle incomplete subblocks)")}
	if( ((ylims[2]-ylims[1])%%block.h) != 0 )
	{warning("height of blocks does not divide height of sampling region evenly, some data will be ommited during subsampling (function does not handle incomplete subblocks)")}
	if(nBoot <= 0)
	{stop("invalid number of bootstraps")}
	if(nBoot < 50)
	{warning("at least 50 bootstrap samples are recommended")}
	if(kappa <= 0)
	{stop("kappa must be positive")}
	if(user.bandwidth == T & length(bandwidth != 2))
	{stop("bandwidth must have length 2")}

	rawdata <- lag_dist_prod(spdata)
	bad <- which(rawdata[,1] == 0 & rawdata[,2] == 0)
	rawdata <- rawdata[-bad,]
	chat.mat <- est_chat_MS(rawdata, lagmat, kappa = 1, user.bandwidth, bandwidth)
	chat <- chat.mat[,3]
	
	blk.chats <- est_block_chats(lagmat, spdata, nBoot, block.dims, xlims, ylims, grid.spacing, kappa, user.bandwidth, bandwidth)
	Vhat <- cov(blk.chats)
	Mt <- A %*% Vhat %*% t(A)
	if(det(Mt) == 0)
	{stop("The matrix A %*% Vhat %*% t(A) is not invertible. Try using a different lag set for the test.")}
	Tn <-  t(A %*% chat) %*% solve(A %*% Vhat %*% t(A)) %*% (A %*% chat)
	Tn <- c(Tn)

	pval.chisq <-  pchisq(Tn, df, lower.tail = F)
	pval.chisq <- c(pval.chisq)

	rv <- list("C.hat" = chat.mat, "V.hat" = Vhat, "n.boot" = nBoot, "test.stat" = Tn, "pvalue.chisq" = pval.chisq )
	
	htestIso.MS <- make_htestIso_MS(rv, df)
	htestIso.MS$data.name <- dname
	return(htestIso.MS)
}
