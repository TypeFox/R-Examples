
#' Continuum Interpolation of Discrete Temperatures Values
#'
#' Builds an interpolator of the temperature at each time and its first
#' derivative.
#' First derivatives are approximated using forward finite
#' differences (\code{\link{approxfun}}). It is assumed that temperature is 0
#' and constant outside the time interval provided.
#'
#' @param temperature_data data frame with the values of the temperatures at
#'        each value of time. It need to have 2 columns, named time and
#'        temperature.
#'
#' @return a list with with two elements: \itemize{
#'      \item temp, the interpolator of the temperature and
#'      \item dtemp, the interpolator of its first derivative
#' }
#'
#' @seealso \code{\link{approxfun}}
#'
#' @importFrom stats approxfun
#'
build_temperature_interpolator <- function(temperature_data){

    temp_profile <- approxfun(temperature_data$time, temperature_data$temperature, rule=2)

    slopes <- with(temperature_data, {

        diff(temperature) / diff(time)

    })

    dtemp_times <- with(temperature_data, {

        c(time[[1]] - 1e-6, time)  # Points are added so dtemp equals zero out of the interval.

    })

    dtemp_profile <- approxfun(dtemp_times, c(0, slopes,0), method="constant", rule=2)

    return(list(temp = temp_profile, dtemp = dtemp_profile))

}
