#--------------------------------------------------------------------
#   npsp-package.R
#--------------------------------------------------------------------
#   npsp-package
#   earthquakes
#   aquifer
#
#   (c) R. Fernandez-Casal         Last revision: Aug 2014
#--------------------------------------------------------------------


#' npsp: Nonparametric spatial (geo)statistics
#'
#' This package implements nonparametric methods which may be useful
#' in geostatistical practice.
#' @section Main functions:
#' \code{\link{locpol}}, \code{\link{np.den}} and \code{\link{np.svar}}
#' use local polynomial kernel methods to compute
#' nonparametric estimates of a multidimensional regression function,
#' a probability density function or a semivariogram (or their first
#' derivatives), respectively.
#' Estimates of these functions can be constructed for any dimension
#' (the amount of available memory is the only limitation).
#' To speed up computations, linear binning is used to discretize the data.
#' A full bandwidth matrix and a multiplicative triweight kernel is used
#' to compute the weights. Main calculations are performed in FORTRAN
#' using the LAPACK library.
#'
#' \code{\link{np.svariso.corr}} computes a bias-corrected nonparametric semivariogram 
#' estimate using an iterative algorithm similar to that described in  
#' Fernandez-Casal and Francisco-Fernandez (2014). This procedure tries to correct
#' the bias due to the direct use of residuals, obtained from a 
#' nonparametric estimation of the trend function, in semivariogram estimation.
#'
#' \code{\link{fitsvar.sb.iso}} fits a `nonparametric' isotropic Shapiro-Botha variogram
#' model by WLS. Currently, only isotropic semivariogram estimation is supported.
#'
#' There are also functions for plotting data joint with a legend representing a
#' continuous color scale. \code{\link{splot}} allows to combine a standard R plot 
#' with a legend. \code{\link{spoints}}, \code{\link{simage}} and \code{\link{spersp}} 
#' draw the corresponding high-level plot with a legend strip for the color scale.
#' These functions are based on \code{\link[fields]{image.plot}} of package \pkg{fields}.
#'
#' Among the other functions intended for direct access by the user,
#' the following could be emphasized: \code{\link{binning}}, \code{\link{bin.den}},
#' \code{\link{svar.bin}}, \code{\link{h.cv}} and \code{\link{interp}}.
#' There are also some functions which can be used to interact with other packages.
#' For instance, \code{\link{as.variogram}} (\pkg{geoR}) or \code{\link{as.vgm}} (\pkg{gstat}).
#'
#' Kriging is not yet implemented in this package. Users are encouraged to use 
#' \code{\link[gstat]{krige}} (or \code{\link[gstat]{krige.cv}}) utilities in \pkg{gstat}
#' package together with \code{\link{as.vgm}}.
#' @author Ruben Fernandez-Casal (Dep. Mathematics, University of A Coru\~{n}a, Spain).
#' Please send comments, error reports or suggestions to \url{rubenfcasal@@gmail.com}.
#' @section Acknowledgments:
#' Important suggestions and contributions to some techniques included here were
#' made by Tomas Cotos-Ya\~{n}ez (Dep. Statistics, University of Vigo, Spain).
#' @name npsp-package
#' @aliases npsp
#' @docType package
#' @useDynLib npsp
#' @import graphics
#' @importFrom quadprog solve.QP
#' @keywords nonparametric smooth
#' @references
#' Fernandez-Casal R. and Francisco-Fernandez M. (2014) 
#' Nonparametric bias-corrected variogram estimation under non-constant trend, 
#' \emph{Stoch. Environ. Res. Ris. Assess}, \bold{28}, 1247-1259.
#'
#' Fernandez-Casal R., Gonzalez-Manteiga W. and  Febrero-Bande M. (2003) 
#' Flexible Spatio-Temporal Stationary Variogram Models, 
#' \emph{Statistics and Computing}, \bold{13}, 127-136.
#'
#' Rupert D. and Wand M.P. (1994) Multivariate locally weighted least squares regression.
#'   \emph{The Annals of Statistics}, \bold{22}, 1346-1370.
#'
#' Shapiro A. and Botha J.D. (1991) Variogram fitting with a general class of
#'   conditionally non-negative definite functions. \emph{Computational Statistics
#'   and Data Analysis}, \bold{11}, 87-96.
#'
#' Wand M.P. (1994) Fast Computation of Multivariate Kernel Estimators.
#'   \emph{Journal of Computational and Graphical Statistics}, \bold{3}, 433-445.
#'
#' Wand M.P. and Jones M.C. (1995) \emph{Kernel Smoothing}. Chapman and Hall, London.
NULL



#' Earthquake data
#'
#' The data set consists of 1859 earthquakes (with magnitude above or equal to
#' 2.0 in Richter's scale), which occurred from 25 November 1944 to 16 October
#' 2013 in the northwest (NW) part of the Iberian Peninsula.
#' The area considered is limited by the coordinates 41N-44N and 6W-10W,
#' which contains the autonomic region of Galicia (Spain) and northern Portugal.
#' @name earthquakes
#' @docType data
#' @format A data frame with 1859 observations on the following 6 variables:
#' \describe{
#'   \item{date}{Date and time (POSIXct format).}
#'   \item{time}{Time (years since first event).}
#'   \item{lon}{Longitude.}
#'   \item{lat}{Latitude.}
#'   \item{depth}{Depth (km).}
#'   \item{mag}{Magnitude (Richter's scale).}
#' }
#' @source National Geographic Institute (IGN) of Spain: \cr
#' \url{http://www.ign.es/ign/layout/sismologiaObtencionDatosSismiscos.do}.
#' @references
#' Francisco-Fernandez M., Quintela-del-Rio A. and Fernandez-Casal R. (2012)
#'    Nonparametric methods for spatial regression. An application to seismic
#'    events, \emph{Environmetrics}, \bold{23}, 85-93.
#' @keywords datasets
#' @examples
#' str(earthquakes)
#' summary(earthquakes)
#' with(earthquakes, spoints(lon, lat, mag, main = "Earthquake data"))
NULL



#' Wolfcamp aquifer data
#'
#' @description 
#' The Deaf Smith County (Texas, bordering New Mexico) was selected as an alternate
#' site for a possible nuclear waste disposal repository in the 1980s.
#' This site was later dropped on grounds of contamination of the aquifer,
#' the source of much of the water supply for west Texas.
#' In a study conducted by the U.S. Department of Energy, piezometric-head data
#' were obtained at 85 locations (irregularly scattered over the Texas panhandle)
#' by drilling a narrow pipe through the aquifer.
#'
#' This data set has been used in numerous papers.
#' For instance, Cressie (1989) lists the data and uses it to illustrate kriging,
#' and Cressie (1993, section 4.1) gives a detailed description of the data
#' and results of different geostatistical analyses.
#' @name aquifer
#' @docType data
#' @format A data frame with 85 observations on the following 3 variables:
#' \describe{
#'   \item{lon}{relative longitude position (miles).}
#'   \item{lat}{relative latitude position (miles).}
#'   \item{head}{piezometric-head levels (feet above sea level).}
#' }
#' @source Harper, W.V. and Furr, J.M. (1986) Geostatistical analysis of
#' potentiometric data in the Wolfcamp Aquifer of the Palo Duro Basin, Texas.
#' \emph{Technical Report BMI/ONWI-587}, Bettelle Memorial Institute, Columbus, OH.
#' @references
#' Cressie, N. (1989) Geostatistics. 
#'   \emph{The American Statistician}, \bold{43}, 197-202. 
#' 
#' Cressie, N. (1993) \emph{Statistics for Spatial Data}. New York. Wiley.
#' @keywords datasets
#' @examples
#' str(aquifer)
#' summary(aquifer)
#' with(aquifer, spoints(lon, lat, head, main = "Wolfcamp aquifer"))
NULL