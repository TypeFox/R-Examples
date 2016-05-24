#' Calculate relative humidity.
#'
#' \code{dewpoint.to.humidity} creates a numeric vector of relative humidity
#'   from numerical vectors of air temperature and dew point temperature.
#'
#' @param dp Numeric vector of dew point temperatures.
#' @param t Numeric vector of air temperatures.
#' @param temperature.metric Character string indicating the temperature
#'    metric of air temperature and dew point temperature. Possible values
#'    are \code{fahrenheit} or \code{celsius}.
#'
#' @details Dew point temperature and temperature must be in the same
#'    metric (i.e., either both in Celsius or both in Fahrenheit). If
#'    necessary, use \code{\link{fahrenheit.to.celsius}} or
#'    \code{\link{celsius.to.fahrenheit}} to convert before using this
#'    function.
#'
#' @return A numeric vector of relative humidity (in \%)
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
#' National Weather Service Hydrometeorological Prediction
#'    Center Web Team. Heat Index Calculator. 30 Jan 2015.
#'    \url{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}.
#'    Accessed 18 Dec 2015.
#'
#' @seealso \code{\link{humidity.to.dewpoint},
#'    \link{fahrenheit.to.celsius},
#'    \link{celsius.to.fahrenheit}}
#'
#' @examples # Calculate relative humidity from air temperature and
#' # dew point temperature.
#'
#' data(lyon)
#' lyon$RH <- dewpoint.to.humidity(t = lyon$TemperatureC,
#'                                 dp = lyon$DewpointC,
#'                                 temperature.metric = 'celsius')
#' lyon
#'
#' @export
dewpoint.to.humidity <-
        function (dp = NA, t = NA,
                  temperature.metric = "fahrenheit")
        {
                if (!(temperature.metric %in% c("celsius", "fahrenheit"))) {
                        stop("The 'temperature.metric' option can onnly by 'celsius' or 'fahrenheit'.")
                }
                if (length(dp) != length(t)) {
                        stop("The vectors for temperature('t') and dewpoint temperature ('dp') must have the same length.")
                }
                if (length(dp[dp > t & !is.na(dp) & !is.na(t)]) > 0) {
                        dp[dp > t] <- NA
                        warning("For some observations, dew point temperature was higher than temperature. Since dew point temperature cannot be higher than air temperature, relative humidty for these observations was set to 'NA'.")
                }
                if (temperature.metric == "fahrenheit") {
                        t <- fahrenheit.to.celsius(t)
                        dp <- fahrenheit.to.celsius(dp)
                }
                beta <- (112 - (0.1 * t) + dp)/(112 + (0.9 * t))
                relative.humidity <- 100 * beta^8
                return(relative.humidity)
        }

#' Calculate dew point temperature.
#'
#' \code{humidity.to.dewpoint} creates a numeric vector of dew point
#'    temperature from numeric vectors of air temperature and relative
#'    humidity.
#'
#' @param rh Numeric vector of relative humidity (in \%).
#' @param t Numeric vector of air temperatures.
#' @param temperature.metric Character string indicating the temperature
#'    metric of air temperature. Possible values are \code{fahrenheit} or
#'    \code{celsius}.
#'
#' @details Dew point temperature will be calculated in the same metric as
#'    the temperature vector (as specified by the \code{temperature.metric}
#'    option). If you'd like dew point temperature in a different metric,
#'     use the function \code{\link{celsius.to.fahrenheit}} or
#'     \code{\link{fahrenheit.to.celsius}} on the output from this function.
#'
#' @return A numeric vector of dew point temperature, in the same metric
#'    as the temperature vector (as specified by the \code{temperature.metric}
#'    option).
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
#' National Weather Service Hydrometeorological Prediction
#'    Center Web Team. Heat Index Calculator. 30 Jan 2015.
#'    \url{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}.
#'    Accessed 18 Dec 2015.
#'
#' @seealso \code{\link{dewpoint.to.humidity},
#'    \link{fahrenheit.to.celsius},
#'    \link{celsius.to.fahrenheit}}
#'
#' @examples # Calculate dew point temperature from relative humidity and
#' # air temperature.
#'
#' data(newhaven)
#' newhaven$DP <- humidity.to.dewpoint(t = newhaven$TemperatureF,
#'                                     rh = newhaven$Relative.Humidity,
#'                                     temperature.metric = 'fahrenheit')
#' newhaven
#'
#' @export
humidity.to.dewpoint <-
        function (rh = NA, t = NA, temperature.metric = "fahrenheit")
        {
                if (!(temperature.metric %in% c("celsius", "fahrenheit"))) {
                        stop("The 'temperature.metric' option can onnly by 'celsius' or 'fahrenheit'.")
                }
                if (length(rh) != length(t)) {
                        stop("The vectors for temperature('t') and relative humidity ('rh') must have the same length.")
                }
                if (length(rh[!is.na(rh) & (rh < 0 | rh > 100)]) > 0) {
                        rh[!is.na(rh) & (rh < 0 | rh > 100)] <- NA
                        warning("For some observations, relative humidity was below 0% or above 100%. Since these values are impossible for relative humidity, dew point temperature for these observations was set to 'NA'.")
                }
                if (temperature.metric == "fahrenheit") {
                        t <- fahrenheit.to.celsius(t)
                }
                dewpoint <- (rh/100)^(1/8) * (112 + (0.9 * t)) - 112 + (0.1 * t)
                if (temperature.metric == "fahrenheit") {
                        dewpoint <- celsius.to.fahrenheit(dewpoint)
                }
                return(dewpoint)
        }
