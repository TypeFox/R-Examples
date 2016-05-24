#' Kernel Smoothing for Bivariate Copula Densities
#' 
#' This package provides fast implementations of kernel estimators for the
#' copula density. Due to its several plotting options it is particularly
#' useful for the exploratory analysis of dependece structures. It can be
#' further used for flexible nonparametric estimation of copula densities
#' and resampling.
#' 
#' The function \code{\link{kdecop}} can be used to estimate a copula density 
#' with a number of popular kernel estimators. The density estimate can be
#' evaluated on arbitrary points with \code{\link[kdecopula:dkdecop]{dkdecop}}; 
#' the cdf with \code{\link[kdecopula:pkdecop]{pkdecop}}. Furthermore, synthetic 
#' data can be simulated with \code{\link[kdecopula:rkdecop]{rkdecop}}, and
#' several plot options are provided by 
#' \code{\link[kdecopula:plot.kdecopula]{plot.kdecopula}}.
#' 
#' @name kdecopula
#' @docType package
#' @author Thomas Nagler
#' 
#' @references
#' Gijbels, I. and Mielniczuk, J. (1990).
#' Estimating the density of a copula function.
#' Communications in Statistics - Theory and Methods, 19(2):445-464. 
#' \cr \cr 
#' Charpentier, A., Fermanian, J.-D., and Scaillet, O. (2006).
#' The estimation of copulas: Theory and practice. 
#' In Rank, J., editor, Copulas: From theory to application in finance. Risk Books.
#' \cr \cr
#' Geenens, G., Charpentier, A., and Paindaveine, D. (2014). 
#' Probit transformation for nonparametric kernel estimation of the copula density.
#' arXiv:1404.4414 [stat.ME]. 
#' \cr \cr 
#' Nagler, T. (2014).
#' Kernel Methods for Vine Copula Estimation.
#' Master's Thesis, Technische Universitaet Muenchen,
#' \url{https://mediatum.ub.tum.de/node?id=1231221}
#' \cr \cr 
#' Wen, K. and Wu, X. (2015).
#' Transformation-Kernel Estimation of the Copula Density,
#' Working paper,
#' \url{http://agecon2.tamu.edu/people/faculty/wu-ximing/agecon2/public/copula.pdf}
#' 
#' @keywords package
#' 
#' @examples
#' 
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' 
#' ## estimation of copula density of variables 5 and 6
#' dens.est <- kdecop(udat[, 5:6])
#' plot(dens.est) 
#' 
#' ## evaluate density estimate at (u1,u2)=(0.123,0.321)
#' dkdecop(c(0.123, 0.321), dens.est) 
#' 
#' ## evaluate cdf estimate at (u1,u2)=(0.123,0.321)
#' pkdecop(c(0.123, 0.321), dens.est) 
#' 
#' ## simulate 500 samples from density estimate
#' rkdecop(500, dens.est)
#' 
NULL


#' Wisconsin Diagnostic Breast Cancer (WDBC)
#' 
#' The data contain measurements on cells in suspicious lumps in a women's
#' breast. Features are computed from a digitized image of a fine needle
#' aspirate (FNA) of a breast mass. They describe characteristics of the cell
#' nuclei present in the image. All samples are classsified as either
#' \emph{benign} or \emph{malignant}.
#' 
#' Ten real-valued features are computed for each cell nucleus: \cr
#' 
#' a) radius (mean of distances from center to points on the perimeter) \cr b)
#' texture (standard deviation of gray-scale values) \cr c) perimeter \cr d)
#' area \cr e) smoothness (local variation in radius lengths) \cr f)
#' compactness (perimeter^2 / area - 1.0) \cr g) concavity (severity of concave
#' portions of the contour) \cr h) concave points (number of concave portions
#' of the contour) \cr i) symmetry \cr j) fractal dimension ("coastline
#' approximation" - 1) \cr
#' 
#' The references listed below contain detailed descriptions of how these
#' features are computed.
#' 
#' The mean, standard error, and "worst" or largest (mean of the three largest
#' values) of these features were computed for each image, resulting in 30
#' features.
#' 
#' @name wdbc
#' @docType data
#' 
#' 
#' @format \code{wdbc} is a \code{data.frame} with 31 columns. The first column
#' indicates wether the sample is classified as benign (\code{B}) or malignant
#' (\code{M}). The remaining columns contain measurements for 30 features.
#' @note This breast cancer database was obtained from the University of
#' Wisconsin Hospitals, Madison from Dr. William H. Wolberg.
#' @references O. L. Mangasarian and W. H. Wolberg: "Cancer diagnosis via
#' linear programming",\cr SIAM News, Volume 23, Number 5, September 1990, pp 1
#' & 18.
#' 
#' William H. Wolberg and O.L. Mangasarian: "Multisurface method of pattern
#' separation for medical diagnosis applied to breast cytology", \cr
#' Proceedings of the National Academy of Sciences, U.S.A., Volume 87, December
#' 1990, pp 9193-9196.
#' 
#' K. P. Bennett & O. L. Mangasarian: "Robust linear programming discrimination
#' of two linearly inseparable sets",\cr Optimization Methods and Software 1,
#' 1992, 23-34 (Gordon & Breach Science Publishers).
#' 
#' 
#' @source
#' \url{https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)}
#' \cr \cr 
#' Bache, K. & Lichman, M. (2013).
#' UCI Machine Learning Repository.
#' Irvine, CA: University of California, School of Information and Computer
#' Science.
#' 
#' @keywords datasets
#' 
#' @examples
#' 
#' data(wdbc)
#' str(wdbc)
#' 
NULL



