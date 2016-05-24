
##' Estimates of material properties by Time-Temperature Superposition (TTS)
##' analysis
##' 
##' TTS analysis is often applied to frequency modulated data obtained by
##' Dynamic Mechanic Analysis (DMA) and Rheometry in the analytical chemistry
##' and physics areas. These techniques provide estimates of material
##' mechanical properties (such as moduli) at different temperatures in a wider
##' range of time. This package provides the Time-Temperature superposition
##' Master Curve at a referred temperature by the three methods: the two wider
##' used methods, Arrhenius based methods and WLF, and the newer methodology
##' based on derivatives procedure. The Master Curve is smoothed by B-splines
##' basis. The package output is composed of plots of experimental data,
##' horizontal and vertical shifts, TTS data, and TTS data fitted using
##' B-splines with bootstrap confidence intervals.
##' 
##' \tabular{ll}{ Package: \tab TTS\cr Type: \tab Package\cr Version: \tab
##' 1.0\cr Date: \tab 2015-09-14\cr License: \tab GPL >= 2\cr } The main
##' functions and data frame are \code{TTS}, \code{PLOT.TTS} and \code{PC}
##' 
##' @name TTS-package
##' @docType package
##' @author Antonio Meneses \email{antoniomenesesfreire@@hotmail.com}, Salvador
##' Naya \email{salva@@udc.es} and Javier Tarrio-Saavedra
##' \email{jtarrio@@udc.es}\cr
##' 
##' Maintainer: Antonio Meneses \email{antoniomenesesfreire@@hotmail.com}
#' @import mgcv
#' @import sfsmisc
#' @import splines
##' @references Naya, S., Meneses A., Tarrio-Saavedra, J., Artiaga R.,
##' Lopez-Beceiro, J. and Gracia-Fernandez C. (2013) New method for estimating
##' shift factors in time-temperatura superposition models. Journal of Thermal
##' Analysis and Calorimetry. ISSN 1388-6150. DOI 10.1007/s10973-013-3193-1.\cr
##' 
##' Williams, M. L. (1964) Structural analysis of Viscoelastic materials. AIAA
##' Journal, 785-808.\cr
##' 
##' Zou, J., You F., Su L., Yang Z., Chen G. and Guo S. (2011). Failure
##' Mechanism of Time-Temperature Superposition for Poly(vinyl
##' chloride)/Dioctylphthalate (100/70) System. DOI 10.1002/app.35113.\cr
##' 
##' Ferry J.D. (1980) Viscoelastic Properties of Polymers, Wiley: New York.\cr
##' 
##' Artiaga R., Garcia A. Fundamentals of DMA. In: 'Thermal analysis.
##' Fundamentals and applications to material characterization' (ed.: Artiaga
##' R.) Publicaciones de la Universidade da Coruna, A Coruna, Spain, 183-206
##' (2005).\cr
##' 
##' Chartoff R.P., Menczel J.D., Dillman S.H. Dynamic mechanical analysis
##' (DMA).  In: 'Thermal analysis of polymers. Fundamentals and applications'
##' (eds.: Menczel J.D., Prime R.B.) Wiley, San Jose, 387-496 (2009).\cr
##' @keywords package
NULL




##' Dataset obtained from polycarbonate (polymer) tests using Dynamic
##' Mechanical Analysis (DMA)
##'
##' PC contains 49 rows and 3 columns.
##'
##' The dataset corresponds to the storage modulus viscoelastic property of
##' different specimens of polycarbonate (PC) and obtained by DMA using TA
##' Instruments Q800 (Naya et al., 2013).
##'
##' @name PC
##' @docType data
##' @format This data frame is composed of the following columns:
#' @format A data frame with XXX observations on the following 3 variables:
#' \describe{
#'   \item{log10.frequency}{It accounts for seven different frequencies
##' (rad/s) in logarithmic scale for each temperature (overall 49).}
##' \item{log10.module}{It accounts for seven different storage
##' modulus, E' (Pa), in base-ten logarithmic scale for each temperature
##' (overall 49).}
##' \item{temperature}{Seven different temperatures:
##' 147, 148, 149, 150, 151, 152, 153 degrees celsius, each one with the
##' corresponding seven values of frequency and storage modulus (overall 49).}
##' }



##' @source Naya, S., Meneses A,. Tarrio-Saavedra, J,. Artiaga R.,
##' Lopez-Beceiro, J. and Gracia-Fernandez C. (2013) New method for estimating
##' shift factors in time-temperatura superposition models. Journal of Thermal
##' Analysis and Calorimetry. ISSN 1388-6150. DOI 10.1007/s10973-013-3193-1.
##' @keywords datasets
##' @examples
##'
## library(TTS)
##' data(PC)
##'
NULL

