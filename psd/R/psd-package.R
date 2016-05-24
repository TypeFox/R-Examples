#' @title Adaptive power spectral density estimation using optimal sine multitapers
#' 
#' @description
#' Estimate the power spectral density (PSD) 
#' of a timeseries using the sine multitapers, adaptively; the number of tapers 
#' (and hence the resolution and uncertainty) vary according to 
#' spectral shape. The main function to be used is \code{\link{pspectrum}}.
#'
#' @details
#' In frequency ranges where the spectrum  (\eqn{S})
#' is relatively flat, more tapers are taken and so a higher accuracy is 
#' attained at the expense of lower frequency resolution. 
#' The program makes a pilot estimate of the spectrum, then uses
#' Riedel and Sidorenko's (1995) estimate of the MSE (minimum square error), 
#' which is based on an estimate of the second derivative of the PSD (\eqn{S''}). 
#' The process is repeated \code{niter} times; further iteration may be necessary 
#' to reach convergence, or an acceptably low spectral variance. 
#' In this context the term "acceptable" is rather subjective: one can 
#' usually detect an unconverged state by a rather jagged appearence of the spectrum,
#' but this is uncommon in our experience.
#'
#' \subsection{Adaptive estimation}{
#' The adaptive process used is as follows. A quadratic fit to the logarithm of the
#' PSD within an adaptively determined frequency band is used to find an estimate of the local second 
#' derivative of the spectrum. This is used in an equation like R-S equation (13) for 
#' the MSE taper number, with the difference that a parabolic weighting is applied with 
#' increasing taper order. Because the FFTs of the tapered series can be found by 
#' resampling the FFT of the original time series (doubled in length and padded with zeros) 
#' only one FFT is required per series, no matter how many tapers are used. 
#' The spectra associated with the sine tapers are weighted before averaging with a 
#' parabolically varying weight. The expression for the optimal number of tapers 
#' given by R-S must be modified since it gives an unbounded result near points 
#' where \eqn{S''} vanishes, which happens at many points in most spectra. 
#' This program restricts the rate of growth of the number of tapers so that a 
#' neighboring covering interval estimate is never completely contained in the next 
#' such interval.
#' }
#'
#' \subsection{Resolution and uncertainty}{
#' The sine multitaper adaptive process 
#' introduces a variable resolution and error in the frequency domain. 
#' See documentation for \code{\link{spectral_properties}} details on
#' how these are computed.
#' }
#'
#' @docType package
#' @name psd-package
#' @aliases psd spec.psd
#' 
#' @author Andrew J. Barbour <andy.barbour@@gmail.com> and Robert L. Parker
#' 
#' 
#' @useDynLib psd
#' @exportPattern "^[[:alpha:]]+"
#' @import RColorBrewer signal zoo
#' @importFrom Rcpp evalCpp
#'
#'
#' @references Barbour, A. J. and R. L. Parker, (2014), 
#' psd: Adaptive, sine multitaper power spectral density estimation for R,
#' \emph{Computers and Geosciences}, \strong{63}, 1--8,
#' \url{http://dx.doi.org/10.1016/j.cageo.2013.09.015}
#'
#' @references Percival, D. B., and A.T. Walden (1993),
#' Spectral analysis for physical applications,
#' \emph{Cambridge University Press}
#'
#' @references Prieto, G. A., R. L. Parker, D. J. Thomson, F. L. Vernon, and R. L. Graham  (2007), 
#' Reducing the bias of multitaper spectrum estimates,
#' \emph{Geophysical Journal International}, \strong{171}, 1269--1281,
#' \url{http://gji.oxfordjournals.org/content/171/3/1269}
#' 
#' @references Riedel, K. S., & Sidorenko, A. (1995), 
#' Minimum bias multiple taper spectral estimation,
#' \emph{Signal Processing, IEEE Transactions on}, \strong{43}(1), 188--195.
#'
#' @seealso \code{\link{psdcore}} and \code{\link{riedsid}}
#'  
NULL
.psdEnvName = ".psdEnv"
.psdEnv = new.env()

##
## Datasets
##

#' A single line of Project MAGNET horizontal field intensity
#' 
#' The Project MAGNET mission 
#' provided a wealth of airborne-magnetometer data
#' spanning the globe (Coleman, 1992).  
#' This dataset represents a single track of horizontal field
#' intensities (a very small subset of the full collection!).
#'
#' \subsection{Raw and Clean Sets}{
#' There are non-real data points in raw MAGNET series; these are 
#' instrumental artefacts, and can severely affect
#' power spectral density (PSD) estimates.  
#' A clean series has been included
#' so that a comparison of PSDs may be made.
#'
#' Some command like \code{subset(magnet, abs(mdiff) > 0)}
#' can be used to identify the rows where edits have been made.
#' }
#' 
#' @name magnet
#' @docType data
#' @format A dataframe with 2048 observations on the following 4 variables.
#'
#' \describe{
#' \item{\code{km}}{Relative along-track distance, in kilometers. The first observation is at zero kilometers.}
#' \item{\code{raw}}{Raw intensities, in nanotesla.}
#' \item{\code{clean}}{Edited raw intensites, in nanotesla}
#' \item{\code{mdiff}}{The difference between \code{clean} and \code{raw} intensities, in nanotesla.}
#' }
#'
#' @seealso \code{\link{pspectrum}}, \code{\link{Tohoku}}, \code{\link{hfsnm}}
#'
#' @references Coleman, R. J. (1992),
#' Project Magnet high-level vector survey data reduction. 
#' In \emph{Types and Characteristics of Data for Geomagnetic Field Modeling},
#' \strong{3153}, pp. 215-248.
#' 
#' @source Project MAGNET page: \url{http://www.ngdc.noaa.gov/geomag/proj_mag.shtml}
#' @keywords datasets
#' @examples
#' data(magnet)
#' summary(magnet)
NULL

#' Observations of teleseismic strains from the 2011 Tohoku earthquake.
#'
#' The \eqn{M_w 9} Tohoku earthquake happend on March 11, 2011.  The seismic
#' waves were recorded at stations across the globe, including by strainmeters
#' in the Plate Boundary Observatory (PBO) borehole strainmeters.
#'
#' These data are for station B084, which is located approximately 8500 km away from
#' the epicenter. Because this distance is large, the seismic waves didn't arrive
#' at this station for more than 700 seconds after the origin time.  So there
#' is a record of pre-seismic noise included, the timeseries extends 6784 seconds
#' prior to the origin time, and 9215 seconds after.  
#'
#' The data are classified with the \code{"epoch"} variable, which separates
#' the series into pre-seismic and seismic data; this is defined relative
#' to the predicted P-wave arrival time from a traveltime model.
#'
#' The original dataset contained \code{NA} values, which were imputed
#' using \code{zoo::na.locf}, which fills \code{NA} with the last previous observation.
#'
#' @name Tohoku
#' @docType data
#' @format A dataframe with 16000 observations on the following 15 variables.
#'
#' \describe{
#' \item{\code{Dts}}{The original datetime string, in UTC.}
#' \item{\code{areal}}{Areal strains}
#' \item{\code{areal.tide}}{Tidal correction to the areal strains.}
#' \item{\code{areal.baro}}{Barometric correction to the areal strains.}
#' \item{\code{gamma1}}{Engineering differential extensional strain: \eqn{\gamma_1}}
#' \item{\code{gamma1.tide}}{Tidal correction for the \eqn{\gamma_1} strains.}
#' \item{\code{gamma1.baro}}{Barometric pressure correction to the \eqn{\gamma_1} strains.}
#' \item{\code{gamma2}}{Engineering shear strain: \eqn{\gamma_2}}.
#' \item{\code{gamma2.tide}}{Tidal correction for the \eqn{\gamma_2} strains.}
#' \item{\code{gamma2.baro}}{Barometric pressure correction to the \eqn{\gamma_2} strains.}
#' \item{\code{pressure.atm}}{Atmospheric pressure.}
#' \item{\code{pressure.pore}}{Pore-fluid pressure.}
#' \item{\code{Dt}}{The \code{Dts} information converted to POSIX datetime.}
#' \item{\code{Origin.secs}}{The number of seconds relative to the earthquake-origin time.}
#' \item{\code{epoch}}{Classification based on predicted P-wave arrival: preseismic or seismic.}
#' }
#'
#' and 2 attributes:
#' 
#' \describe{
#' \item{\code{units}}{A list of strings regarding the units of various physical quantities given here.}
#' \item{\code{iasp}}{A list of source and station characteristics, including the
#' the origin time, predicted
#' traveltimes for P and S waves, and the geodetic information used in the traveltime
#' calculation.}
#' }

#' @seealso \code{\link{pspectrum}}, \code{\link{hfsnm}}, \code{\link{magnet}}
#' @seealso \code{TauP.R} for an R-implementation of the traveltime calculations:
#' @seealso \url{http://cran.r-project.org/web/packages/TauP.R/}
#' @keywords datasets
#'
#' @references USGS summary page: 
#' @references \url{http://earthquake.usgs.gov/earthquakes/eqinthenews/2011/usc0001xgp/}
#' @source PBO High Frequency archive: 
#' @source \url{http://borehole.unavco.org/bsm/earthquakes/NeartheEastCoastofHonshuJapan_20110311}
#'
#' @examples
#' data(Tohoku)
#' str(Tohoku)
NULL

#' Noise levels found in PBO strainmeter data at seismic frequencies.
#'
#' These values represent noise levels in high frequency data (\eqn{10^{-3} - 10} Hz) from 2009, averaged over all
#' stations in the Anza cluster of the Plate Boundary Observatory (PBO) borehole
#' strainmeter network, and the UCSD-style longbase laser strainmeters.
#'
#' \code{NA} values in the series highlight frequency bands where the noise
#' levels are unreliable, due to a instrumental artifact.
#'
#' @name hfsnm
#' @docType data
#' @format A dataframe with 141 observations on the following 4 variables:
#'
#' \describe{
#' \item{\code{freq}}{Frequencies, in Hertz.}
#' \item{\code{P50}}{The 50th percentile (median) noise levels in decibels relative to \eqn{1 \epsilon^2 / } Hz.}
#' \item{\code{P10}}{The 10th percentile noise levels also in decibels.}
#' \item{\code{meter.type}}{The strainmeter design type.}
#' }
#'
#' and 2 attributes:
#' \describe{
#' \item{\code{source.doi}}{The DOI number of the source publication.}
#' \item{\code{generator}}{The structure of a function which will refresh the values from the supplemental files of the original publication.}
#' }
#'
#' @seealso \code{\link{pspectrum}}, \code{\link{Tohoku}}, \code{\link{magnet}}
#' @keywords datasets
#' @source Barbour, A. J., and Agnew, D. C. (2011), Noise Levels on Plate Boundary Observatory Borehole Strainmeters in Southern California,
#' \emph{Bulletin of the Seismological Society of America},
#' \strong{101}(5), 2453-2466, doi:10.1785/0120110062
#'
#' @examples
#' data(hfsnm)
#' str(hfsnm)
#' FUN <- attr(hfsnm, "generator")
#' try(dat <- FUN(molten=FALSE)) # may fail without library-access to BSSA
#' try(all.equal(dat[,1:4], hfsnm[,1:4]))
NULL
