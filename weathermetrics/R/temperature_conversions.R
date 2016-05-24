#' Convert from Celsius to Fahrenheit.
#'
#' \code{celsius.to.fahrenheit} creates a numeric vector of temperatures in
#'    Fahrenheit from a numeric vector of temperatures in Celsius.
#'
#' @param T.celsius Numeric vector of temperatures in Celsius.
#' @param round Integer indicating the number of decimal places to
#'     round converted value.
#'
#' @return A numeric vector of temperature values in Fahrenheit.
#'
#' @note Equations are from the source code for the US National Weather
#'     Service's
#'     \href{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}{online heat index calculator}.
#'
#' @author
#' Brooke Anderson \email{brooke.anderson@@colostate.edu},
#' Roger Peng \email{rdpeng@@gmail.com}
#'
#' @seealso \code{\link{fahrenheit.to.celsius}}
#'
#' @examples # Convert from Celsius to Fahrenheit.
#' data(lyon)
#' lyon$TemperatureF <- celsius.to.fahrenheit(lyon$TemperatureC)
#' lyon
#'
#' @export
celsius.to.fahrenheit <-
        function (T.celsius, round = 2)
        {
                T.fahrenheit <- (9/5) * T.celsius + 32
                T.fahrenheit <- round(T.fahrenheit, digits = round)
                return(T.fahrenheit)
        }

#' Convert from Fahrenheit to Celsius.
#'
#' \code{fahrenheit.to.celsius} creates a numeric vector of temperatures in
#'    Celsius from a numeric vector of temperatures in Fahrenheit.
#'
#' @param T.fahrenheit Numeric vector of temperatures in Fahrenheit.
#' @param round Integer indicating the number of decimal places to
#'     round converted value.
#'
#' @return A numeric vector of temperature values in Celsius.
#'
#' @note Equations are from the source code for the US National Weather
#'     Service's
#'     \href{http://www.wpc.ncep.noaa.gov/html/heatindex.shtml}{online heat index calculator}.
#'
#' @author
#' Brooke Anderson \email{brooke.anderson@@colostate.edu},
#' Roger Peng \email{rdpeng@@gmail.com}
#'
#' @seealso \code{\link{celsius.to.fahrenheit}}
#'
#' @examples # Convert from Fahrenheit to Celsius.
#' data(norfolk)
#' norfolk$TempC <- fahrenheit.to.celsius(norfolk$TemperatureF)
#' norfolk
#'
#' @export
fahrenheit.to.celsius <-
        function (T.fahrenheit, round = 2)
        {
                T.celsius <- (5/9) * (T.fahrenheit - 32)
                T.celsius <- round(T.celsius, digits = round)
                return(T.celsius)
        }
