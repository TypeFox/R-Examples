######################################################################
## CODE AND DOCUMENTATION FOR PACKAGE lunar
######################################################################

#' @name lunar-package
#' @docType package
#' @description
#'   Provides functions to calculate the phase of the moon,
#'   its distance from the earth, the season and possibly
#'   other environmental factors, based on date and location.
#' @details
#'   This package is used by the author to calculate covariates
#'   for studies of lunar effect on health and healthcare.
#'
#'   References forthcoming.
NULL

#' @title
#'   Lunar Phase Categories (8)
#' @description
#'   Return 8 category labels for lunar phases.
#' @details
#'   These are category names corresponding to phases of
#'   the moon.
#'   Moon phase category names may be returned in
#'   the output of the
#'   \code{\link{lunar.phase}} function if its
#'   \code{name} option is set to \code{TRUE}.
#' @keywords lunar moon phase
#' @examples
#' print(lunar.8phases)
#' @seealso \code{\link{lunar.phase}}
#' @export
lunar.8phases <- c("New",
                   "Waxing crescent",
                   "First quarter",
                   "Waxing gibbous",
                   "Full",
                   "Waning gibbous",
                   "Last quarter",
                   "Waning crescent")

#' @title
#'   Lunar Phase Categories (4)
#' @description
#'   Return 4 category labels for lunar phases.
#' @details
#'   These are category names corresponding to phases of
#'   the moon.
#'   Moon phase category names may be returned in
#'   the output of the
#'   \code{\link{lunar.phase}} function if its
#'   \code{name} option is set to \code{TRUE}.
#' @keywords lunar moon phase
#' @examples
#' print(lunar.4phases)
#' @seealso \code{\link{lunar.phase}}
#' @export
lunar.4phases <- c("New",
                   "Waxing",
                   "Full",
                   "Waning")

#' @title
#'   Lunar Distance Categories
#' @description
#'   Return category labels for lunar distances.
#' @details
#'   These are category names corresponding to distances
#'   from the center of the Earth to the center of its moon.
#'   Distance category names are used in the output of the
#'   \code{\link{lunar.distance}} function if its
#'   \code{name} option is set to \code{TRUE}.
#'
#'   The perigee occurs at roughly \eqn{56} Earth radii, and
#'   the apogee at about \eqn{62.8} Earth radii.
#'   These categories are not determined according to any
#'   common standard.
#'   They may have different precise definitions
#'   for their bounds in different analyses.
#' @examples
#' print(lunar.distances)
#' @keywords lunar moon distance
#' @seealso \code{\link{lunar.distance}}
#' @export
lunar.distances <- c("Apogee",
                     "Far",
                     "Average",
                     "Near",
                     "Perigee")

#' @title
#'   Terrestrial Season Categories
#' @description
#'   Return category labels for the seasons on Earth.
#' @details
#'   These are category names corresponding to the
#'   seasons of the planet Earth.
#' @examples
#' print(terrestrial.seasons)
#' @keywords earth season
#' @export
terrestrial.seasons <- c("Winter",
                         "Spring",
                         "Summer",
                         "Autumn")
                     
#' @title
#'   Lunar Phase
#' @description
#'   Returns the lunar phase on specified dates.
#' @details
#'   Adapted from function 
#'   \code{moon.illumination} in from the
#'   R4MFCL project (not an R package), which was developed
#'   by the Secretariat of the Pacific Community (SPC).
#'   The R4MFCL project was led by Simon Hoyle, and also includes code
#'   by Shelton Harley, Nick Davies, and Adam Langley of the SPC,
#'   and Pierre Kleiber of the US National Marine Fisheries Service.
#'   Pierre Kleiber is the author of the
#'   \code{moon.illumination} function.
#'
#'   Code from project R4MFCL is distributed under the MIT License:
#'
#'     \url{http://opensource.org/licenses/mit-license.php}
#'
#'   Here is a link to the R4MFCL project:
#'
#'     \url{https://code.google.com/p/r4mfcl/}
#'
#'   The R4MFCL code was modified as follows:
#'   \itemize{
#'     \item Changed function name from \code{moonphase} to \code{lunar.phase}.
#'     \item Changed input date to be in \code{\link[base]{Date}} format
#'           (as opposed to \code{\link[base]{POSIXct}}).
#'     \item Removed reliance on other R4MFCL functions.
#'     \item Changed name of primary input from \code{ptime} to \code{x}.
#'     \item Added optional \code{shift} term (in hours) relative to 12h UT.
#'     \item Added optional \code{name} term to control whether phase names
#'           as opposed to radians should be returned.
#'     \item Changed the documentation.
#'   }
#'
#'   Where radians are returned:
#'   \itemize{
#'     \item 0 refers to the new moon
#'     \item \eqn{\pi/2} refers to the first quarter
#'     \item \eqn{\pi} refers to the full moon
#'     \item \eqn{3\pi/2} refers to the last quarter
#'   }
#'
#'   Adapted originally from Stephen R. Schmitt: Sky & Telescope,
#'   Astronomical Computing, April 1994 and
#'   \url{http://mysite.verizon.net/res148h4j/zenosamples/zs_lunarphasecalc.html},
#'   which references
#'   Jean Meeus, Astronomical Algorithms. Willmann-Bell, Inc. (1991) 429p.
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @param shift
#'   The number of hours by which to shift the calculation of lunar phase.
#'   By default lunar phase is calculated at 12 noon UT.
#' @param name
#'   Optional parameter indicating whether the return is a factor variable
#    consisting of a lunar phase label, or the lunar phase in radians.
#'   By default lunar phase is returned in radians.
#'   If assigned the value 8, it returns a factor variable with
#'   8 phase levels.
#'   If \code{TRUE} or any value other than 0 or 8,
#'   it returns a factor variable with 4 phase labels.
#' @param ...
#'   Other optional arguments are ignored.
#' @examples
#' lunar.phase(as.Date("2013-05-06"))
#' @keywords lunar moon phase
#' @seealso \code{\link{lunar.4phases}}
#' @seealso \code{\link{lunar.8phases}}
#' @export
lunar.phase <- function(
  x,
  shift = 0,
  ...,
  name = FALSE
) {
  
  ##========================================================================
  ## Returns lunar phase in radians on a given date.
  ## 0 => new; pi/2 => first quarter; pi => full; 3*pi/2 => last quarter.
  ## Adapted from Stephen R. Schmitt: Sky & Telescope, Astronomical
  ## Computing, April 1994 and
  ## http://mysite.verizon.net/res148h4j/zenosamples/zs_lunarphasecalc.html
  ## which refs.
  ## Jean Meeus. 1991. Astronomical Algorithms. Willmann-Bell, Inc. 429p.
  ##========================================================================
  ## Adapted from function moonphase() in R project R4MFCL, which
  ## was developed by the Secretariat of the Pacific Community (SPC).
  ## The R4MFCL project was led by Simon Hoyle, and also includes code
  ## by Shelton Harley, Nick Davies, and Adam Langley of the SPC,
  ## and Pierre Kleiber of the US National Marine Fisheries Service.
  ##
  ## Code from project R4MFCL is distributed under the MIT License:
  ##   http://opensource.org/licenses/mit-license.php
  ## Here is a link to the R4MFCL project:
  ##   https://code.google.com/p/r4mfcl/
  ##========================================================================
  ## Modified as follows:
  ##   * changed function name from moonphase() to lunar.phase()
  ##   * changed input date to be in Date format (as opposed to POSIXct)
  ##   * removed reliance on other functions in R4MFCL
  ##   * changed name of primary input from ptime to x
  ##   * added optional shift term (in hours) relative to 12h UT
  ##   * added optional name term to control whether phase names
  ##     as opposed to radians should be returned
  ##   * changed documentation
  ##========================================================================

  Y <- as.numeric(format(x, "%Y"))
  M <- as.numeric(format(x, "%m"))
  D <- as.numeric(format(x, "%d"))

  ## calculate the Julian date at 12h UT + shift
  YY <- Y - floor( ( 12 - M ) / 10 )
  MM <- M + 9
  MM <- ifelse(MM >= 12, MM-12, MM)

  K1 <- floor( 365.25 * ( YY + 4712 ) )
  K2 <- floor( 30.6001 * MM + 0.5 )
  K3 <- floor(floor((YY/100)+49)*0.75) - 38

  ## for dates in Julian calendar
  JD <- K1 + K2 + D + 59 + (shift / 24)
  ## adjust for Gregorian calendar
  JD <- ifelse(JD > 2299160, JD - K3, JD)

  IP <- (JD - 2451550.1)/29.530588853
  IP <- IP-floor(IP)
  IP <- ifelse(IP<0, IP+1, IP)

  if(name) {
    ## get moon age in days
    MC <- 29.53
    AG <- IP * MC

    if(name == 8) {
      MC8 <- MC/16
      retval <- factor(
              ifelse(AG < MC8, lunar.8phases[1],
              ifelse(AG < 3*MC8, lunar.8phases[2],
              ifelse(AG < 5*MC8, lunar.8phases[3],
              ifelse(AG < 7*MC8, lunar.8phases[4],
              ifelse(AG < 9*MC8, lunar.8phases[5],
              ifelse(AG < 11*MC8, lunar.8phases[6],
              ifelse(AG < 13*MC8, lunar.8phases[7],
              ifelse(AG < 15*MC8, lunar.8phases[8],
                  lunar.8phases[1])))))))),
          levels = lunar.8phases)
    } else {
      MC4 <- MC/8
      retval <- factor(
              ifelse(AG < MC4, lunar.4phases[1],
              ifelse(AG < 3*MC4, lunar.4phases[2],
              ifelse(AG < 5*MC4, lunar.4phases[3],
              ifelse(AG < 7*MC4, lunar.4phases[4],
                  lunar.4phases[1])))),
          levels = lunar.4phases)
    }
  } else {
    ## Convert phase to radians
    retval <- 2 * pi * IP
  }
  
  return(retval)
}

#' @title
#'   Lunar Illumination
#' @description
#'   Returns the proportion of lunar illumination on specified dates.
#' @details
#'   Adapted from function 
#'   \code{moon.illumination} in from the
#'   R4MFCL project (not an R package), which was developed
#'   by the Secretariat of the Pacific Community (SPC).
#'   The R4MFCL project was led by Simon Hoyle, and also includes code
#'   by Shelton Harley, Nick Davies, and Adam Langley of the SPC,
#'   and Pierre Kleiber of the US National Marine Fisheries Service.
#'   Pierre Kleiber is the author of the
#'   \code{moon.illumination} function.
#'
#'   Code from project R4MFCL is distributed under the MIT License:
#'
#'     \url{http://opensource.org/licenses/mit-license.php}
#'
#'   Here is a link to the R4MFCL project:
#'
#'     \url{https://code.google.com/p/r4mfcl/}
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @param shift
#'   The number of hours by which to shift the calculation of lunar phase.
#'   By default lunar phase is calculated at 12 noon UT.
#' @examples
#' lunar.illumination(as.Date("2004-03-24"))
#' @keywords lunar moon illumination light
#' @seealso \code{\link{lunar.illumination.mean}}
#' @export
lunar.illumination <- function(
  x,
  shift = 0
) {
  
  ##========================================================================
  ## Returns the proportion of the moon illuminated on a given date.
  ##========================================================================
  ## Adapted from function moon.illumination() in R project R4MFCL, which
  ## was developed by the Secretariat of the Pacific Community (SPC).
  ## The R4MFCL project was led by Simon Hoyle, and also includes code
  ## by Shelton Harley, Nick Davies, and Adam Langley of the SPC,
  ## and Pierre Kleiber of the US National Marine Fisheries Service.
  ##
  ## Code from project R4MFCL is distributed under the MIT License:
  ##   http://opensource.org/licenses/mit-license.php
  ## Here is a link to the R4MFCL project:
  ##   https://code.google.com/p/r4mfcl/
  ##========================================================================
  
  return((1-sin(pi/2-lunar.phase(x, shift = shift)))/2)
}

#' @title
#'   Lunar Distance
#' @description
#'   Returns the distance of the moon from the earth on specified dates.
#' @details
#'   Distance to the moon is returned in units of earth radii, or
#'   as a 5-level factor variable referring to the moon's
#'   perigee (at about \eqn{56} earth radii) and
#'   apogee (at about \eqn{63.8} earth radii).
#'
#'   Adapted from Stephen R. Schmitt: Lunar Phase Computation:
#'   \url{http://mysite.verizon.net/res148h4j/zenosamples/zs_lunarphasecalc.html}.
#'   Last accessed: 1 September 2014.
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @param shift
#'   The number of hours by which to shift the distance calculation.
#'   By default distance is calculated at 12 noon UT.
#' @param name
#'   Optional parameter indicating whether the return is a factor variable
#'   consisting of a lunar distance label, or the lunar distance in earth
#'   radii.  By default lunar phase is returned in earth radii.
#' @param strict
#'   Optional parameter indicating whether the return should employ strict
#'   definitions for distance labels, that is, with apogee and perigee
#'   within 5% of their extreme values.  By default the alternative
#'   definition breaks the distance categories evenly into 20%-iles.
#'   The 'average' category is the same in both definitions.
#' @param ...
#'   Other optional arguments are ignored.
#' @examples
#' lunar.distance(as.Date("2004-03-24"))
#' @keywords lunar moon
#' @seealso \code{\link{lunar.distances}}
#' @export
lunar.distance <- function(
  x,
  shift = 0,
  ...,
  name = FALSE,
  strict = FALSE
) {
  
  ##========================================================================
  ## Returns the distance of the moon from the earth on a given date.
  ##========================================================================
  ## Adapted from Stephen R. Schmitt: Lunar Phase Computation:
  ## http://mysite.verizon.net/res148h4j/zenosamples/zs_lunarphasecalc.html
  ##========================================================================

  Y <- as.numeric(format(x, "%Y"))
  M <- as.numeric(format(x, "%m"))
  D <- as.numeric(format(x, "%d"))

  ## calculate the Julian date at 12h UT + shift
  YY <- Y - floor( ( 12 - M ) / 10 )
  MM <- M + 9
  MM <- ifelse(MM >= 12, MM-12, MM)

  K1 <- floor( 365.25 * ( YY + 4712 ) )
  K2 <- floor( 30.6001 * MM + 0.5 + (shift / 24) )
  K3 <- floor(floor((YY/100)+49)*0.75) - 38

  JD <- K1 + K2 + D + 59 ## for dates in Julian calendar
  JD <- ifelse(JD > 2299160, JD - K3, JD) ## for Gregorian calendar

  IP <- (JD - 2451550.1)/29.530588853
  IP <- IP-floor(IP)
  IP <- ifelse(IP<0, IP+1, IP)
  IP <- 2 * pi * IP
  
  DP <- (JD - 2451562.2 ) / 27.55454988
  DP <- DP-floor(DP)
  DP <- ifelse(DP<0, DP+1, DP)
  DP <- 2 * pi * DP
  
  DI <- 60.4 - 3.3*cos(DP) - 0.6*cos(2*IP - DP) - 0.5*cos(2*IP)

  if(name) {
    if(strict) retval <- factor(
                   ifelse(DI <  56 + (63.8-56)/20, lunar.distances[5],
                   ifelse(DI <  56 + 2 * (63.8-56)/5, lunar.distances[4],
                   ifelse(DI <  56 + 3 * (63.8-56)/5, lunar.distances[3],
                   ifelse(DI <  56 + 19 * (63.8-56)/20, lunar.distances[2],
                       lunar.distances[1])))),
                   levels = lunar.distances)
    else retval <- factor(
             ifelse(DI <  56 + (63.8-56)/5, lunar.distances[5],
             ifelse(DI <  56 + 2 * (63.8-56)/5, lunar.distances[4],
             ifelse(DI <  56 + 3 * (63.8-56)/5, lunar.distances[3],
             ifelse(DI <  56 + 4 * (63.8-56)/5, lunar.distances[2],
                 lunar.distances[1])))),
             levels = lunar.distances)
  } else {
    retval <- DI
  }

  return(retval)
}

#' @title
#'   Average Lunar Metrics
#' @description
#'   Returns an average measurement around specified dates.
#' @details
#'   This in an internal support function that integrates a lunar
#'   measurement over time using step sizes of days or hours.
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @param towards
#'   The directed window size from \code{x} in days.
#'   By default the window looks back 7 days including \code{x}.
#' @param by
#'   The exposure interval and integration basis.
#'   The default is to represent a day's illumination by
#'   the illumination at 12 noon UT.  The other options integrate
#'   midrange illuminations over hours.  Options 'day' and 'night'
#'   are not currently implemented, but will be used to limit
#'   exposure intervals.
#'   The use of an unimplemented option
#'   in a function call will result in a \code{NULL} value
#'   being returned.
#' @param type
#'   Whether illumination or distance metrics are to be returned.
#'   The use of an unimplemented option
#'   in a function call will result in a \code{NULL} value
#'   being returned.
#' @param ...
#'   Other optional arguments are ignored.
#' @examples
#' \dontrun{
#' lunar.metric.mean(as.Date("2004-03-24"), type="illumination")
#' }
#' @keywords lunar moon metric covariate
#' @seealso \code{\link{lunar.illumination}}
#' @seealso \code{\link{lunar.distance}}
lunar.metric.mean <- function(
  x,
  towards = -6,
  ...,
  by = c("date", "hours", "day", "night"),
  type = c("illumination", "distance")
) {

  retval <- NULL

  if("date" %in% by) {
    retval <- unlist(lapply(x, function(X) {
        my.metric.dates <-
            seq(from = X + ifelse(towards <= -1, towards, 0),
               to = X + ifelse(towards >= 1, towards, 0), by = "days")
        if("illumination" %in% type)
          my.metric.values <- lunar.illumination(my.metric.dates)
        else if("distance" %in% type)
          my.metric.values <- lunar.distance(my.metric.dates)
        else return(NULL)
        mean(my.metric.values) }))
  } else if("hours" %in% by) {
    retval <- unlist(lapply(x, function(X) {
        my.metric.dates <-
            seq(from = X + ifelse(towards <= -1, towards, 0),
                to = X + ifelse(towards >= 1, towards, 0), by = "days")
        if("illumination" %in% type)
          my.metric.values <-
            unlist(lapply(my.metric.dates, function(Y) {
                mean(c(lunar.illumination(Y, shift = (-12):11),
                       lunar.illumination(Y, shift = (-11):12))) }))
        else if("distance" %in% type)
          my.metric.values <-
            unlist(lapply(my.metric.dates, function(Y) {
                mean(c(lunar.distance(Y, shift = (-12):11),
                       lunar.distance(Y, shift = (-11):12))) }))
        else return(NULL)
        mean(my.metric.values) }))
  }
  
  return(retval)
}

#' @title
#'   Average Lunar Illumination
#' @description
#'   Returns the average lunar illumination around specified dates.
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @param towards
#'   The directed window size from \code{x} in days.
#'   By default the window looks back 7 days including \code{x}.
#' @param by
#'   The exposure interval and integration basis.
#'   The default is to represent a day's illumination by
#'   the illumination at 12 noon UT.  The other options integrate
#'   midrange illuminations over hours.  Options 'day' and 'night'
#'   are not currently implemented, but will be used to limit
#'   exposure intervals.
#'   The use of an unimplemented option
#'   in a function call will result in a \code{NULL} value
#'   being returned.
#' @param ...
#'   Other optional arguments are ignored.
#' @examples
#' lunar.illumination.mean(as.Date("2004-03-24"))
#' @keywords lunar moon illumination light
#' @seealso \code{\link{lunar.illumination}}
#' @export
lunar.illumination.mean <- function(
  x,
  towards = -6,
  ...,
  by = c("date", "hours", "day", "night")
) {

  lunar.metric.mean(x = x,
                    towards = towards,
                    by = by,
                    type = "illumination")
  
}

#' @title
#'   Average Lunar Distance
#' @description
#'   Returns the average lunar distance around specified dates.
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @param towards
#'   The directed window size from \code{x} in days.
#'   By default the window looks back 7 days including \code{x}.
#' @param by
#'   The exposure interval and integration basis.
#'   The default is to represent a day's distance by
#'   the distance at 12 noon UT.  The other options integrate
#'   midrange distances over hours.
#' @param ...
#'   Other optional arguments are ignored.
#' @examples
#' lunar.distance.mean(as.Date("2004-03-24"))
#' @keywords lunar moon distance
#' @seealso \code{\link{lunar.distance}}
#' @export
lunar.distance.mean <- function(
  x,
  towards = -6,
  ...,
  by = c("date", "hours", "day", "night")
) {

  lunar.metric.mean(x = x,
                    towards = towards,
                    by = by,
                    type = "distance")
  
}

#' @title
#'   Terrestrial Season
#' @description
#'   Returns the season on specified dates.
#' @details
#'   The definitions for non-leap years are as follows
#'   (dates are inclusive):
#'   \describe{
#'     \item{Winter:}{21 December through 21 March}
#'     \item{Spring:}{22 March through 21 June}
#'     \item{Summer:}{22 June through 20 September}
#'     \item{Autumn:}{21 September through 20 December}
#'   }
#'   In leap years spring comes a day early!
#' @param x
#'   A vector of \code{\link[base]{Date}} values.
#' @examples
#' terrestrial.season(as.Date("2004-03-24"))
#' @keywords earth season
#' @seealso \code{\link{terrestrial.seasons}}
#' @export
terrestrial.season <- function(
  x
) {

  # Check that we have a proper Date
  if(!is(x, "Date")) stop("Must provide a Date.")
  
  my.day <- as.POSIXlt(x)$yday
 
  return(factor(
      ifelse(my.day < 81, terrestrial.seasons[1],
          ifelse(my.day < 173, terrestrial.seasons[2],
              ifelse(my.day < 284, terrestrial.seasons[3],
                  ifelse(my.day < 355, terrestrial.seasons[4],
                         terrestrial.seasons[1])))),
      levels = terrestrial.seasons))
}
