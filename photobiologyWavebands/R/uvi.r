#' Calculate UVI from spectral irradiance
#'
#' @param spct a \code{source.spct} object
#' @param std "WWO" (250 nm to 400 nm), "NOAA" (286.5 nm to 400 nm)
#'
#' @return a numeric value for the unitless UVI (This is a value on a continuous scale,
#' rather than the discrete scale normally used.)
#'
#' @description UVI (UV Index) is a unitless quantity based on erythema BSWF,
#'   that gives an easy to interpret UV measure, mainly meant for informing
#'   general public about sunburn risk.
#'
#' @details Two different definitions of UV Index are implemented in this
#'   package. Setting std="NOAA" follows the definition in Kiedron et al. (2007)
#'   but using CIE98 as SWF. NOAA definition discards wavelengths shorter than
#'   286.5 nm as when calculated based on spectral data from Brewer instruments.
#'   "WMO" uses the internationally accepted lower limit at 250 nm (see WHO,
#'   2002). "NOAA" is the default as this is safer with noisy data for solar
#'   radiation measured at ground level, and in this case the value of UVI
#'   should be correct, and almost identical except for errors caused by noise
#'   at shorter wavelengths. However, when calculating UVI from radiation
#'   spectra from UV lamps, "WMO" should be used, as most UV lamps do emit some
#'   radiation between 250 nm and 286.5 nm.
#'
#' @references WHO (2002) Global Solar UV Index: A Practical Guide. ISBN 92 4
#' 159007 6, WHO, Geneva.
#' \url{http://www.who.int/uv/publications/globalindex/en/} or
#' \url{http://uv.biospherical.com/Solar_Index_Guide.pdf}.
#'
#' P. Kiedron, S. Stierle and K. Lantz (2007) Instantaneous UV Index and Daily
#' UV Dose Calculations. NOAA-EPA Brewer Network.
#' \url{http://www.esrl.noaa.gov/gmd/grad/neubrew/docs/UVindex.pdf}
#'
#' @export
#'
#' @examples
#' UVI(sun.spct)
#' UVI(sun.spct, "WMO")
#' round(UVI(sun.spct), 0)

UVI <- function(spct, std="NOAA") {
  e_irrad(spct, w.band=UVI_wb(std), wb.trim=TRUE) * 40.0
}
