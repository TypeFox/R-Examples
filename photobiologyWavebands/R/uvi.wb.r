#' Definition of CIE weighted waveband
#'
#' UVI (UV Index) based on Erythema BSWF
#'
#' @param std "NOAA" or "WMO"
#' @return a list defining the wavelength range, weighting and normalization.
#'   This is just the CIE98 BSWF but with the wavelength limits adjusted to
#'   those used for UVI.
#'
#' @description This is just a convenience function that returns the same
#'   weights as \code{CIE} as defined by \code{CIE_e_fun} but with no option to
#'   change normalization wavelength, and with the wavelength limits adjusted to
#'   those used for UVI. Using std="NOAA" follows the definition in
#'   \url{http://www.esrl.noaa.gov/gmd/grad/neubrew/docs/UVindex.pdf} but using
#'   CIE98 as SWF. NOAA definition discards wavelengths shorter than 286.5 nm as
#'   when calculated based on spectral data from Brewer instruments. "WMO" uses
#'   the internationally accepted lower limit at 250 nm. "NOAA" is the default
#'   as this is safer with noisy data for sunlight, and for solar radiation the
#'   value of UVI should be correct. When calculating UVI for radiation spectra
#'   from UV lamps, "WMO" should be used.
#'
#' @references WHO (2002) Global Solar UV Index: A Practical Guide. ISBN 92 4
#' 159007 6, WHO, Geneva. \url{http://uv.biospherical.com/Solar_Index_Guide.pdf}
#'
#' P. Kiedron, S. Stierle and K. Lantz (2007) Instantaneous UV Index and Daily
#' UV Dose Calculations. NOAA-EPA Brewer Network.
#' \url{http://www.esrl.noaa.gov/gmd/grad/neubrew/docs/UVindex.pdf}
#'
#' @seealso \code{\link{UVI}}
#'
#' @keywords internal

UVI_wb <- function(std="NOAA") {
  label <- "UVIndex"
  if (std=="NOAA") {w.low=286.5} else
  if (std=="WMO" | std=="WHO") {w.low=250.0} else
  {
    w.low=250.0
    warning("Unrecognized value for UVI std, using 'WMO'.")
  }
  new_waveband(w.low=w.low, w.high=400,
               weight="SWF", SWF.e.fun=CIE_e_fun, SWF.norm=298.0,
               norm=298.0, hinges=c(249.99, 250, 298, 328, 399.99, 400),
               wb.name=paste(label, std, sep="."), wb.label=label)
}
