#' Calculate heat index.
#'
#' \code{heat.index} creates a numeric vector of heat index values from
#'    numeric vectors of air temperature and either relative humidity or
#'    dew point temperature.
#'
#' @param t Numeric vector of air temperatures.
#' @param dp Numeric vector of dew point temperatures.
#' @param rh Numeric vector of relative humidity (in \%).
#' @param temperature.metric Character string indicating the temperature
#'    metric of air temperature and dew point temperature. Possible values
#'    are 'fahrenheit' or 'celsius'.
#' @param output.metric Character string indicating the metric into which
#'    heat index should be calculated. Possible values are 'fahrenheit' or
#'    'celcius'.
#' @param round Integer indicating the number of decimal places to
#'     round converted value.
#'
#' @details Include air temperature (\code{t}) and either dew point
#'    temperature (\code{dp}) or relative humdity (\code{rh}). You cannot
#'    specify both dew point temperature and relative humidity- this will
#'    return an error. Heat index is calculated as \code{NA} when impossible
#'    values of dew point temperature or humidity are input (e.g., humidity
#'    above 100\% or below 0\%, dew point temperature above air temperature).
#'
#' @return A numeric vector of heat index values in the metric specified
#'    by \code{output.metric}. (If \code{output.metric} is not specified,
#'    heat index will be returned in the same metric in which air
#'    temperature was input, specified by \code{temperature.metric}.)
#'
#' @note Equations are from the source code for the US National Weather
#'     Service's
#'     \href{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}{online heat index calculator}.
#'
#' @author
#' Brooke Anderson \email{brooke.anderson@@colostate.edu},
#' Roger Peng \email{rdpeng@@gmail.com}
#'
#' @references
#' Anderson GB, Bell ML, Peng RD. 2013. Methods to calculate the heat index
#'    as an exposure Metric in environmental health research.
#'    Environmental Health Perspectives 121(10):1111-1119.
#'
#' National Weather Service Hydrometeorological Prediction
#'    Center Web Team. Heat Index Calculator. 30 Jan 2015.
#'    \url{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}.
#'    Accessed 18 Dec 2015.
#'
#' Rothfusz L. 1990. The heat index (or, more than you ever wanted to know
#'    about heat index) (Technical Attachment SR 90-23). Fort Worth:
#'    Scientific Services Division, National Weather Service.
#'
#' R. Steadman, 1979. The assessment of sultriness. Part I: A
#'    temperature-humidity index based on human physiology and clothing
#'    science. Journal of Applied Meteorology, 18(7):861--873.
#'
#' @examples # Calculate heat index from temperature (in Fahrenheit)
#'# and relative humidity.
#'
#' data(suffolk)
#'suffolk$heat.index <- heat.index(t = suffolk$TemperatureF,
#'                                  rh = suffolk$Relative.Humidity)
#'suffolk
#'
#'# Calculate heat index (in Celsius) from temperature (in
#'# Celsius) and dew point temperature (in Celsius).
#'
#'data(lyon)
#'lyon$heat.index <- heat.index(t = lyon$TemperatureC,
#'                               dp = lyon$DewpointC,
#'                               temperature.metric = 'celsius',
#'                               output.metric = 'celsius')
#'lyon
#'
#' @export
heat.index <-
        function (t = NA, dp = c(), rh = c(), temperature.metric = "fahrenheit",
                  output.metric = NULL, round = 0)
        {
                if(length(output.metric) == 0){
                        output.metric <- temperature.metric
                }
                if (length(dp) == 0 & length(rh) == 0) {
                        stop("You must give values for either dew point temperature ('dp') or relative humidity ('rh').")
                }
                else if (length(dp) > 0 & length(rh) > 0) {
                        stop("You can give values for either dew point temperature ('dp') or relative humidity ('rh'), but you cannot specify both to this function.")
                }
                if (length(dp) != length(t) & length(rh) != length(t)) {
                        stop("The vectors for temperature ('t') and moisture (either relative humidity, 'rh', or dew point temperature, 'dp') must be the same length.")
                }
                if (length(dp) > length(rh)) {
                        rh <- dewpoint.to.humidity(t = t, dp = dp, temperature.metric = temperature.metric)
                }
                else if (length(rh[!is.na(rh) & (rh > 100 | rh < 0)]) > 0) {
                        rh[!is.na(rh) & (rh > 100 | rh < 0)] <- NA
                        warning("There were observations with an impossible values for relative humidity (below 0% or above 100%). For these observations, heat index was set to NA.")
                }
                if (temperature.metric == "celsius") {
                        t <- celsius.to.fahrenheit(t)
                }
                hi <- mapply(heat.index.algorithm, t = t, rh = rh)
                if (output.metric == "celsius") {
                        hi <- fahrenheit.to.celsius(hi)
                }
                hi <- round(hi, digits = round)
                return(hi)
        }

#' Algorithm for heat.index function.
#'
#' \code{heat.index.algorithm} converts a numeric scalar of temperature
#'    (in Fahrenheit) and a numeric scalar of relative humidity (in \%)
#'    to heat index (in Fahrenheit). This function is not meant to be used
#'    outside of the \code{\link{heat.index}} function.
#'
#' @param t Numeric scalar of air temperature, in Fahrenheit.
#' @param rh Numeric scalar of relative humidity, in \%.
#'
#' @details If an impossible value of relative humidity is given
#'    (below 0\% or above 100\%), heat index is returned as \code{NA}.
#'
#' @return A numeric scalar of heat index, in Fahrenheit.
#'
#' @note Equations are from the source code for the US National Weather
#'     Service's
#'     \href{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}{online heat index calculator}.
#'
#' @author
#' Brooke Anderson \email{brooke.anderson@@colostate.edu},
#' Roger Peng \email{rdpeng@@gmail.com}
#'
#' @references
#' Anderson GB, Bell ML, Peng RD. 2013. Methods to calculate the heat index
#'    as an exposure Metric in environmental health research.
#'    Environmental Health Perspectives 121(10):1111-1119.
#'
#' National Weather Service Hydrometeorological Prediction
#'    Center Web Team. Heat Index Calculator. 30 Jan 2015.
#'    \url{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}.
#'    Accessed 18 Dec 2015.
#'
#' Rothfusz L. 1990. The heat index (or, more than you ever wanted to know
#'    about heat index) (Technical Attachment SR 90-23). Fort Worth:
#'    Scientific Services Division, National Weather Service.
#'
#' R. Steadman, 1979. The assessment of sultriness. Part I: A
#'    temperature-humidity index based on human physiology and clothing
#'    science. Journal of Applied Meteorology, 18(7):861--873.
#'
#' @seealso \code{\link{heat.index}}
#'
#' @export
heat.index.algorithm <-
        function (t = NA, rh = NA)
        {
                if (is.na(rh) | is.na(t)) {
                        hi <- NA
                } else if (t <= 40) {
                        hi <- t
                } else {
                        alpha <- 61 + ((t - 68) * 1.2) + (rh * 0.094)
                        hi <- 0.5*(alpha + t)
                        if (hi > 79) {
                                hi <- -42.379 + 2.04901523 * t + 10.14333127 * rh -
                                        0.22475541 * t * rh - 6.83783 * 10^-3 * t^2 -
                                        5.481717 * 10^-2 * rh^2 + 1.22874 * 10^-3 * t^2 *
                                        rh + 8.5282 * 10^-4 * t * rh^2 - 1.99 * 10^-6 *
                                        t^2 * rh^2
                                if (rh <= 13 & t >= 80 & t <= 112) {
                                        adjustment1 <- (13 - rh)/4
                                        adjustment2 <- sqrt((17 - abs(t - 95))/17)
                                        total.adjustment <- adjustment1 * adjustment2
                                        hi <- hi - total.adjustment
                                } else if (rh > 85 & t >= 80 & t <= 87) {
                                        adjustment1 <- (rh - 85)/10
                                        adjustment2 <- (87 - t)/5
                                        total.adjustment <- adjustment1 * adjustment2
                                        hi <- hi + total.adjustment
                                }
                        }
                }
                return(hi)
        }
