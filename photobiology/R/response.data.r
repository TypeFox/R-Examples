#' @title Spectral response of a GaAsP photodiode
#'
#' @description A dataset containing wavelengths at a 1 nm interval and
#'   spectral response as $A W^{-1} nm^{-1}$ for GaAsP photodiode type
#'   G6262 from Hamamatsu. Data digitized from manufacturer's data sheet.
#'
#' @details \itemize{ \item w.length (nm). \item s.e.response
#' (A/W)  }
#'
#' @references
#' Hamamatsu (2011) Datasheet: GaAsP Photodiodes G5645 G5842 G6262. Hamamatsu
#' Photonics KK, Hamamatsu, City.
#' http://www.hamamatsu.com/resources/pdf/ssd/g5645_etc_kgpd1004e.pdf.
#' Visited 2016-02-01.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{response_spct} object with 311 rows and 2 variables
#' @name photodiode.spct
NULL

#' @title Spectral response of a back-thinned CCD image sensor.
#'
#' @description A dataset containing wavelengths at a 1 nm interval and
#'   spectral response as $A W^{-1} nm^{-1}$ for CCD sensor type
#'   S11071/S10420 from Hamamatsu (measured without a quartz window). These
#'   arrays are frequently used as sensors in high-UV-sensitivity array
#'   spectrometers. Data digitized from manufacturer's data sheet.
#'
#' @details \itemize{ \item w.length (nm). \item s.e.response
#' (A/W)  }
#'
#' @references
#' Hamamatsu (2014) Datasheet: CCD Image Sensors S11071/S10420-01 Series.
#' Hamamatsu Photonics KK, Hamamatsu, City.
#' http://www.hamamatsu.com/resources/pdf/ssd/s11071-1004_etc_kmpd1120e.pdf.
#' Visited 2016-02-01.
#'
#' @docType data
#' @keywords datasets
#' @format A \code{response_spct} object with 201 rows and 2 variables
#' @name ccd.spct
NULL
