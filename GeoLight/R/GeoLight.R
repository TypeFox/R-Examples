#' @name GeoLight-package
#' @aliases GeoLight
#' @title The GeoLight Package
#' @description This is a summary of all features of \bold{\code{GeoLight}}, a \code{R}-package for 
#' analyzing light based geolocator data
#' @details \bold{\code{GeoLight}} is a package to derive geographical positions from daily light intensity pattern. 
#' Positioning and calibration methods are based on the threshold-method (Ekstrom 2004, Lisovski \emph{et al.} 2012). 
#' A changepoint model from the \code{R} package \code{changepoint} is implemented to distinguish between periods of 
#' residency and movement based on the sunrise and sunset times. Mapping functions are implemented 
#' using the \code{R} package \code{maps}.
#' @section Getting Started: 
#' We refrain from giving detailed background on the (several steps of) 
#' analysis of light-based geolocator data here but strongly recommend the key-publications below. 
#' @section Updates:
#' We advise all users to update their installation of \bold{\code{GeoLight}} regularly.
#' Type \code{news(package="GeoLight")} to read news documentation about changes to the recent and all previous version of the package
#' @section Important notes:
#' Most functions in \bold{\code{GeoLight}} require the same initial units and mostly the format and object type is mandatory:
#' \tabular{rll}{
#' \tab \code{tFirst} \tab yyyy-mm-dd hh:mm "UTC" (see: \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones})\cr
#' \tab \code{tSecond} \tab as \emph{tFirst} (e.g. 2008-12-01 17:30) \cr
#' \tab \code{type} \tab either 1 or 2 depending on wheter \emph{tFirst} is sunrise (1) or sunset (2)\cr
#' \tab \code{coord} \tab \code{SpatialPoints} or a \code{matrix}, containing x and y coordinates (in that order) \cr
#' \tab \code{degElevation} \tab a \code{vector} or a single \code{value} of sun elevation angle(s) in degrees (e.g. -6)
#' }
#' @section FUNCTIONS AND DATASETS:
#' In the following, we give a summary of the main functions and sample datasets in the \bold{\code{GeoLight}} package. 
#' Alternatively a list of all functions and datasets in alphabetical order is available by typing \code{library(help=GeoLight)}. 
#' For further information on any of these functions, type \code{help(function name)}.
#'
#'  @section CONTENTS:
#' \tabular{rll}{
#' \tab {I.} \tab Determination of sunset and sunrise \cr
#' \tab {II.} \tab Residency analysis \cr
#' \tab {III.} \tab Calibration \cr
#' \tab {IV.} \tab {Positioning}\cr
#' \tab {V.} \tab {Data visualisation}\cr
#' \tab {VI.} \tab {Examples}
#' }
#' @section I. Determination of sunset and sunrise:
#' \tabular{rll}{
#' \tab \code{\link{gleTrans}} \tab transformation of already defined twilight events* \cr
#' \tab \code{\link{glfTrans}} \tab transformation of light intensity measurements over time* \cr
#' \tab \code{\link{luxTrans}} \tab transformation of light intensity measurements over time** \cr
#' \tab \code{\link{lightFilter}} \tab filter to remove noise in light intensity measurements during the night \cr
#' \tab \code{\link{twilightCalc}} \tab definition of twilight events (\emph{sunrise, sunset}) from light intensity measurements \cr
#' }
#' 
#' * written for data recorded by geolocator devices from the \bold{Swiss Ornithological Institute} \cr
#' ** written for data recorded by geolocator devices from \bold{Migrate Technology Ltd}
#' 
#' @section II. Residency Analysis:
#' \tabular{rll}{
#' \tab \code{\link{changeLight}} \tab function to distinguish between residency and movement periods \cr
#' \tab \code{\link{schedule}} \tab function to produce a data frame summerizing the residency and movement pattern \cr
#' }
#' 
#' @section III. Calibration:
#' See Lisovski \emph{et al.} 2012 for all implemented calibration methods.
#' \tabular{rll}{
#' \tab \code{\link{getElevation}} \tab function to calculate the sun elevation angle for data with known position \cr
#' \tab \code{\link{HillEkstromCalib}} \tab \emph{Hill-Ekstrom calibration} for one or more defined stationary periods \cr
#' }
#' 
#' @section IV. Positioning:
#' \tabular{rll}{
#' \tab \code{\link{coord}} \tab main function to derive a \code{matrix} of spatial coordinates \cr
#' \tab \code{\link{distanceFilter}} \tab filter function to reduce unrealistic positions (not recommended, since the filtering ignore positioning error) \cr
#' \tab \code{\link{loessFilter}} \tab filter function to define outliers in sunrise and sunset times (defined twilight events) \cr
#' }
#' 
#' @section V. Data visualisation:
#' \tabular{rll}{
#' \tab \code{\link{tripMap}} \tab function to map the derived positions and combine the coordinates in time order\cr
#' \tab \code{\link{siteMap}} \tab function to show the results of the residency analysis on a map \cr
#' }
#' 
#' @section IV. Examples:
#' \tabular{rll}{
#' \tab \code{\link{calib1}} \tab data for calibration: light intensities \cr
#' \tab \code{\link{calib2}} \tab  data for calibration: Calculated twilight events (from \code{\link{calib1}} by \code{\link{twilightCalc}}) \cr
#' \tab \code{\link{hoopoe1}} \tab light intensity measurements over time recorded on a migratory bird \cr
#' \tab \code{\link{hoopoe2}} \tab sunrise and sunset times: From light intensity measurement (from \code{\link{hoopoe1}}) \cr
#' }
#' 
#' @section \bold{\code{R}} Packages for Further Spatial Ananlyses:
#' \code{spatstat} \cr
#' \code{adehabitat} \cr
#' \code{gstat} \cr
#' \code{trip} \cr
#' \code{tripEstimation} \cr
#' \code{move} \cr
#' ...
#' 
#' @section Acknowledgements:
#' Steffen Hahn, Felix Liechti, Fraenzi Korner-Nievergelt, Andrea Koelzsch, Eldar Rakhimberdiev, Erich Baechler, Eli Bridge,  Andrew Parnell, Richard Inger
#'
#' @section Authors:
#' Simeon Lisovski, Simon Wotherspoon, Michael Sumner, Silke Bauer, Tamara Emmenegger \cr
#' Maintainer: Simeon Lisovski <simeon.lisovski(at)gmail.com>
#' 
#' @section References:
#' Ekstrom, P.A. (2004) An advance in geolocation by light. \emph{Memoirs of the National Institute of Polar Research}, Special Issue, \bold{58}, 210-226.
#' 
#' Fudickar, A.M., Wikelski, M., Partecke, J. (2011) Tracking migratory songbirds: accuracy of light-level loggers (geolocators) in forest habitats. \emph{Methods in Ecology and Evolution}, DOI: 10.1111/j.2041-210X.2011.00136.x.
#'
#' Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step: Latitude. \emph{Electronic Tagging and Tracking in Marine Fisheries} (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The Netherlands.
#'
#' Hill, R.D. (1994) Theory of geolocation by light levels. \emph{Elephant Seals: Population Ecology, Behavior, and Physiology} (eds L. Boeuf, J. Burney & R.M. Laws), pp. 228-237. University of California Press, Berkeley.
#'  
#' Lisovski, S. and Hahn, S. (2012) GeoLight - processing and analysing light-based geolocator data in R. \emph{Methods in Ecology and Evolution}, doi: 10.1111/j.2041-210X.2012.00248.x
#' 
#' Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and precision affected by environmental factors. \emph{Methods in Ecology and Evolution}, doi: 10.1111/j.2041-210X.2012.00185.x
#' 
#' Wilson, R.P., Ducamp, J.J., Rees, G., Culik, B.M. & Niekamp, K. (1992) Estimation of location: global coverage using light intensity. \emph{Wildlife telemetry: remote monitoring and tracking of animals} (eds I.M. Priede & S.M. Swift), pp. 131-134. Ellis Horward, Chichester.
NULL


i.argCheck <- function(y) {
  
  if(any(sapply(y, function(x) class(x))=="data.frame")) {
    ind01 <- which(sapply(y, function(x) class(x))=="data.frame")
    if(!all(ind02 <- c("tFirst", "tSecond", "type")%in%names(y[[ind01]]))) {
      whc <- paste("The following columns in data frame twl are missing with no default: ", paste(c("tFirst", "tSecond", "type")[!ind02], collapse = " and "), sep = "")
      stop(whc , call. = F)
    } 
    out <- y[[ind01]]
  } else {
    if(!all(c("tFirst", "tSecond", "type")%in%names(y))) {
      ind03 <- c("tFirst", "tSecond", "type")%in%names(y)
      stop(sprintf(paste(paste(c("tFirst", "tSecond", "type")[!ind03], collapse = " and "), "is missing with no default.")))
    } else {
      out <- data.frame(tFirst = y$tFirst, tSecond = y$tSecond, type = y$type)
    }
  }
  if(any(c(class(out[,1])[1], class(out[,2])[1])!="POSIXct")) {
    stop(sprintf("Date and time inforamtion (e.g. tFirst and tSecond) need to be provided as POSIXct class objects."), call. = F)
  }
  out  
  
}


##' Estimate location from consecutive twilights
##'
##' This function estimates the location given the times at which 
##' the observer sees two successive twilights.
##' 
##' Longitude is estimated by computing apparent time of local noon
##' from sunrise and sunset, and determining the longitude for which
##' this is noon. Latitude is estimated from the required zenith and
##' the sun's hour angle for both sunrise and sunset, and averaged.
##'
##' When the solar declination is near zero (at the equinoxes)
##' latitude estimates are extremely sensitive to errors.  Where the
##' sine of the solar declination is less than \code{tol}, the
##' latitude estimates are returned as \code{NA}.
##' 
##' The format (date and time) of \emph{tFirst} and \emph{tSecond} has to be
##' "yyyy-mm-dd hh:mm" corresponding to Universal Time Zone UTC (see:
##' \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones})
##' 
##' 
##' @title Simple Threshold Geolocation Estimates
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param degElevation the sun elevation angle (in degrees) that defines twilight (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst} or \code{nrow(x)}.
##' @param tol tolerance on the sine of the solar declination (only implemented in method 'NOAA').
##' @param note \code{logical}, if TRUE a notation of how many positions could
##' be calculated in proportion to the number of failures will be printed at the
##' end.
##' @param method Defines the method for the location estimates. 'NOAA' is based on
##' code and the excel spreadsheet from the NOAA site (http://www.esrl.noaa.gov/gmd/grad/solcalc/),
##' 'Montenbruck' is based on Montenbruck, O. & Pfleger, T. (2000) Astronomy on the Personal
##' Computer. \emph{Springer}, Berlin.
##' @return A matrix of coordinates in decimal degrees. First column are
##' longitudes, expressed in degrees east of Greenwich. Second column contains
##' the latitudes in degrees north the equator.
##' @author Simeon Lisovski, Simon Wotherspoon, Michael Sumner
##' @examples
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' crds <- coord(hoopoe2, degElevation=-6, tol = 0.2)
##' ## tripMap(crds, xlim=c(-20,20), ylim=c(5,50), main="hoopoe2")
##' @export   
coord  <- function(tFirst, tSecond, type, twl, degElevation = -6, tol = 0, method = "NOAA",  note = TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  rise <- ifelse(tab$type==1, tab$tFirst, tab$tSecond)
  set <- ifelse(tab$type==1, tab$tSecond, tab$tFirst)
  
  if(method == "NOAA") {
    rad <- pi/180
    sr <- solar(rise)
    ss <- solar(set)
    cosz <- cos(rad*(90-degElevation))
    lon <- -(sr$solarTime+ss$solarTime+ifelse(sr$solarTime<ss$solarTime,360,0))/2
    lon <- (lon+180)%%360-180
    
    ## Compute latitude from sunrise
    hourAngle <- sr$solarTime+lon-180
    a <- sr$sinSolarDec
    b <- sr$cosSolarDec*cos(rad*hourAngle)
    x <- (a*cosz-sign(a)*b*suppressWarnings(sqrt(a^2+b^2-cosz^2)))/(a^2+b^2)
    lat1 <- ifelse(abs(a)>tol,asin(x)/rad,NA)
    
    ## Compute latitude from sunset
    hourAngle <- ss$solarTime+lon-180
    a <- ss$sinSolarDec
    b <- ss$cosSolarDec*cos(rad*hourAngle)
    x <- (a*cosz-sign(a)*b*suppressWarnings(sqrt(a^2+b^2-cosz^2)))/(a^2+b^2)
    lat2 <- ifelse(abs(a)>tol,asin(x)/rad,NA)
    
    ## Average latitudes
    out <- cbind(lon=lon,lat=rowMeans(cbind(lat1,lat2),na.rm=TRUE))
  } 
  
  if(method == "Montenbruck") {
    
    out <- coord2(tab$tFirst, tab$tSecond, tab$type, degElevation)
    
  }
  
  
  if(note) cat(paste("Note: Out of ", nrow(out)," twilight pairs, the calculation of ", sum(is.na(out[,2]))," latitudes failed ","(",
                     floor(sum(is.na(out[,2])*100)/nrow(out))," %)",sep=""))
  out
}

coord2 <- function(tFirst, tSecond, type, degElevation=-6) {
  
  # if noon, RadHourAngle = 0, if midnight RadHourAngle = pi
  # --------------------------------------------------------
  RadHourAngle <- numeric(length(type))
  index1 <- type==1
  if (sum(index1)>0) RadHourAngle[index1] <- 0
  RadHourAngle[!index1] <- pi
  # --------------------------------------------------------
  
  tSunTransit <- as.character(as.POSIXct(tFirst, tz = "GMT") + as.numeric(difftime(as.POSIXct(tSecond, tz = "GMT"), 
                                                                                   as.POSIXct(tFirst, tz = "GMT"), 
                                                                                   units="secs")/2))
  
  index0 <- (nchar(tSunTransit) <= 10)
  if(sum(index0)>0) tSunTransit[index0] <- as.character(paste(tSunTransit[index0]," ","00:00",sep=""))
  
  # Longitude
  jD  <- i.julianDate(as.numeric(substring(tSunTransit,1,4)),as.numeric(substring(tSunTransit,6,7)),
                      as.numeric(substring(tSunTransit,9,10)),as.numeric(substring(tSunTransit,12,13)),as.numeric(substring(tSunTransit,15,16)))
  jD0 <- i.julianDate(as.numeric(substring(tSunTransit,1,4)),as.numeric(substring(tSunTransit,6,7)),
                      as.numeric(substring(tSunTransit,9,10)),rep(0,length(tFirst)),rep(0,length(tFirst)))
  jC  <- i.JC2000(jD)
  jC0 <- i.JC2000(jD0)
  
  radObliquity         <- i.radObliquity(jC)
  radEclipticLongitude <- i.radEclipticLongitude(jC)
  radRightAscension    <- i.radRightAscension(radEclipticLongitude,radObliquity)
  radGMST              <- i.radGMST(jD,jD0,jC,jC0)
  
  degLongitude <- i.setToRange(-180,180,i.deg(RadHourAngle + radRightAscension - radGMST))
  
  
  # Latitude
  jD  <- i.julianDate(as.numeric(substring(tFirst,1,4)),as.numeric(substring(tFirst,6,7)),
                      as.numeric(substring(tFirst,9,10)),as.numeric(substring(tFirst,12,13)),
                      as.numeric(substring(tFirst,15,16)))
  jD0 <- i.julianDate(as.numeric(substring(tFirst,1,4)),as.numeric(substring(tFirst,6,7)),
                      as.numeric(substring(tFirst,9,10)),rep(0,length(tFirst)),rep(0,length(tSecond)))
  jC  <- i.JC2000(jD)
  jC0 <- i.JC2000(jD0)
  
  
  radElevation         <- if(length(degElevation)==1) rep(i.rad(degElevation),length(jD)) else i.rad(degElevation)
  
  sinElevation         <- sin(radElevation)
  radObliquity         <- i.radObliquity(jC)
  radEclipticLongitude <- i.radEclipticLongitude(jC)
  radDeclination       <- i.radDeclination(radEclipticLongitude,radObliquity)
  sinDeclination       <- sin(radDeclination)
  cosDeclination       <- cos(radDeclination)
  
  radHourAngle         <- i.radGMST(jD,jD0,jC,jC0) + i.rad(degLongitude) - i.radRightAscension(radEclipticLongitude,radObliquity)
  cosHourAngle         <- cos(radHourAngle)
  
  term1                <- sinElevation/(sqrt(sinDeclination^2 + (cosDeclination*cosHourAngle)^2))
  term2                <- (cosDeclination*cosHourAngle)/sinDeclination
  
  degLatitude <- numeric(length(radElevation))
  
  index1  <- (abs(radElevation) > abs(radDeclination))
  if (sum(index1)>0) degLatitude[index1] <- NA
  
  index2  <- (radDeclination > 0 & !index1)
  if (sum(index2)>0)  degLatitude[index2] <- i.setToRange(-90,90,i.deg(asin(term1[index2])-atan(term2[index2])))
  
  index3  <- (radDeclination < 0 & !index1)
  if (sum(index3)>0)  degLatitude[index3] <- i.setToRange(-90,90,i.deg(acos(term1[index3])+atan(1/term2[index3])))
  
  degLatitude[!index1&!index2&!index3] <- NA
  
  cbind(degLongitude, degLatitude)
}



##' Function to calculate the median sun elevation angle for light measurements at a
##' known location and the choosen light threshold.
##' 
##' Optionally, shape and scale paramters of the twiligth error (in minutes) can be estimated. The error is assumed
##' to follow a log-normal distribution and 0 (elev0) is set 0.1 below the minimum sun elevation angle of estimated twilight times.
##' Those parameters might be of interest for sensitivity analysis or further processing using the R Package SGAT (https://github.com/SWotherspoon/SGAT).
##'
##' @title Calculate the appropriate sun elevation angle for known location
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param known.coord a \code{SpatialPoint} or \code{matrix} object, containing
##' known x and y coordinates (in that order) for the selected measurement
##' period.
##' @param plot \code{logical}, if TRUE a plot will be produced.
##' @param lnorm.pars \code{logical}, if TRUE shape and scale parameters of the twilight error (log-normal distribution) 
##' will be estimated and included in the output (see Details).
##' @author Simeon Lisovski
##' @references Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt,
##' F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
##' precision affected by environmental factors. \emph{Methods in Ecology and
##' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
##' @examples
##' data(calib2)
##'   calib2$tFirst  <- as.POSIXct(calib2$tFirst, tz = "GMT")
##'   calib2$tSecond <- as.POSIXct(calib2$tSecond, tz = "GMT")
##' getElevation(calib2, known.coord = c(7.1,46.3), lnorm.pars = TRUE)
##' @importFrom MASS fitdistr
##' @importFrom graphics arrows par hist plot lines mtext
##' @importFrom stats dlnorm median
##' @export getElevation
getElevation <- function(tFirst, tSecond, type, twl, known.coord, plot=TRUE, lnorm.pars = FALSE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  tab <- geolight.convert(tab[,1], tab[,2], tab[,3])  
  
  sun <- solar(as.POSIXct(tab[,1], "GMT"))
  z   <- 90-refracted(zenith(sun, known.coord[1], known.coord[2]))
  
  tab$z.tm <- as.POSIXct("1900-01-01 00:00:01", "GMT")
  tab$z.tm[tab[,2]] <- twilight(tab[tab[,2], 1], known.coord[1], known.coord[2], rise = TRUE, zenith = (min(z)-0.1)-90 ,iters = 3) 
  tab$z.tm[!tab[,2]] <- twilight(tab[!tab[,2],1], known.coord[1], known.coord[2], rise = FALSE, zenith = (min(z)-0.1)-90 ,iters = 3)
  
  tab$diff <- NA 
  
  tab$diff[tab[,2]] <- as.numeric(difftime(tab[tab[,2],1], tab[tab[,2],3], units = "mins"))
  tab$diff[!tab[,2]] <- as.numeric(difftime(tab[!tab[,2],3], tab[!tab[,2],1], units = "mins"))
  
  
  if(plot) {
    opar <- par(mfrow = c(1, 2), mar = c(7, 7, 5, 1), oma = c(0, 0, 0, 2), 
                cex.lab = 1.5, cex.axis = 1.5, las = 1, mgp = c(4.8, 2, 1))
    
    hist(z, breaks =  seq(min(z)-0.5, max(z)+0.5, length = 18), main = "", 
         col = "grey60", xlab = "Sun elevaion angle")
    arrows(median(z), -0.75, median(z), -0.1,  lwd = 3, col = "cornflowerblue", xpd = T)
    mtext(paste("median elevation", round(median(z),2)), font = 3, col = "cornflowerblue", cex = 1.2, line = 1.6)
    
    hist(tab$diff, freq = F, breaks = seq(min(tab$diff)-2, max(tab$diff)+2, length = 18), 
         main = "", xlab = "Twilight error (minutes)", col = "grey95", ylab = "Probability density")
    
    if(lnorm.pars) {
      seq1 <- seq(0, max(tab$diff), length = 100)
      fit <- fitdistr(tab$diff, "log-Normal")
      lines(seq1, dlnorm(seq1, fit$estimate[1], fit$estimate[2]), col = "firebrick", lwd = 3, lty = 2)
      
      mtext(paste("elev0 =", round((min(z)-0.1), 2), "\nshape =", round(fit$estimate[1],2), 
                  "\nscale =", round(fit$estimate[2],2)), 
            font = 3, col = "firebrick", cex = 1.2, line = 0.5)
    }
    par(opar)
  }
  
  if(lnorm.pars) c(med.elev=median(z), shape = as.numeric(fit$estimate[1]), 
                   scale = as.numeric(fit$estimate[2])) else c(med.elev = median(z))
}


##' Residency analysis using a changepoint model
##'
##' Function to discriminate between periods of residency and movement based on
##' consecutive sunrise and sunset data. The calculation is based on a
##' changepoint model (\bold{\pkg{R}} Package \code{\link{changepoint}}:
##' \code{\link{cpt.mean}}) to find multiple changepoints within the
##' data.
##'
##' The \code{cpt.mean} from the \code{R} Package \code{changepoint} is a
##' function to find a multiple changes in mean for data where no assumption is
##' made on their distribution. The value returned is the result of finding the
##' optimal location of up to Q changepoints (in this case as many as possible)
##' using the cumulative sums test statistic.
##'
##'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param quantile probability threshold for stationary site selection. Higher
##' values (above the defined quantile of all probabilities) will be considered
##' as changes in the behavior. Argmuent will only be considered if either \code{rise.prob} and/or
##' \code{set.prob} remain unspecified.
##' @param rise.prob the probability threshold for \bold{sunrise}: greater or
##' equal values indicates changes in the behaviour of the individual.
##' @param set.prob the probability threshold for \bold{sunset}: higher and
##' equal values indicates changes in the behaviour of the individual.
##' @param days a threshold for the length of stationary period. Periods smaller
##' than "days" will not be considered as a residency period
##' @param plot logical, if \code{TRUE} a plot will be produced
##' @param summary logical, if \code{TRUE} a summary of the results will be
##' printed
##' @return A \code{list} with probabilities for \emph{sunrise} and
##' \emph{sunset} the user settings of the probabilities and the resulting
##' stationary periods given as a \code{vector}, with the residency sites as
##' positiv numbers in ascending order (0 indicate movement/migration).
##' @note The sunrise and/or sunset times shown in the graph (if
##' \code{plot=TRUE}) represent hours of the day. However if one or both of the
##' twilight events cross midnight during the recording period the values will
##' be formed to avoid discontinuity.
##' @author Simeon Lisovski & Tamara Emmenegger
##' @seealso \code{\link{changepoint}}, \code{\link{cpt.mean}}
##' @references Taylor, Wayne A. (2000) Change-Point Analysis: A Powerful New
##' Tool For Detecting Changes.
##'
##' M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis.
##' \emph{Wiley}.
##'
##' Chen, J. and Gupta, A. K. (2000) Parametric statistical change point
##' analysis. \emph{Birkhauser}.
##' @examples
##'
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' residency <- changeLight(hoopoe2, quantile=0.9)
##'
##' @importFrom changepoint cpt.mean cpts.full pen.value.full
##' @importFrom stats na.omit quantile
##' @importFrom graphics abline axis layout mtext par plot rect
##' @export changeLight
changeLight <- function(tFirst, tSecond, type, twl, quantile=0.9, rise.prob=NA, set.prob=NA, days=5, plot=TRUE, summary=TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type = c(tab$type, ifelse(tab$type == 1, 2, 1)))
  tw <- tw[!duplicated(tw$datetime),]
  tw <- tw[order(tw[,1]),]
  
  hours <- as.numeric(format(tw[,1],"%H"))+as.numeric(format(tw[,1],"%M"))/60
  
  for(t in 1:2){
    cor <- rep(NA, 24)
    for(i in 0:23){
      cor[i+1] <- max(abs((c(hours[tw$type==t][1],hours[tw$type==t])+i)%%24 -
                            (c(hours[tw$type==t],hours[tw$type==t][length(hours)])+i)%%24),na.rm=T)
    }
    hours[tw$type==t] <- (hours[tw$type==t] + (which.min(round(cor,2)))-1)%%24
  }
  
  sr <- tw[tw[,2]==1,1]
  ss <- tw[tw[,2]==2,1]
  
  rise <- hours[tw[,2]==1]
  set <- hours[tw[,2]==2]
  
  # start: Change Point Model
  # max. possible Change Points (length(sunrise)/2)
  CPs1 <- suppressWarnings(cpt.mean(rise, method= "BinSeg", Q=length(rise)/2, penalty="Manual",pen.value=0.001, 
                                    test.stat = "CUSUM",param.estimates=FALSE))
  CPs2 <- suppressWarnings(cpt.mean(set, method="BinSeg", Q=length(set)/2, penalty="Manual",pen.value=0.001, 
                                    test.stat = "CUSUM",param.estimates=FALSE))
  
  N1 <- seq(1,length(rise))
  N2 <- seq(1,length(set))
  
  tab1 <- merge(data.frame(N=N1,prob=NA),data.frame(N=cpts.full(CPs1)[nrow(cpts.full(CPs1)),],
                                                    prob=pen.value.full(CPs1)/2),by.x="N",by.y="N",all.x=T)[,-2]
  tab1[is.na(tab1[,2]),2] <- 0
  tab2 <- merge(data.frame(N=N2,prob=NA),data.frame(N=cpts.full(CPs2)[nrow(cpts.full(CPs2)),],
                                                    prob=pen.value.full(CPs2)/2),by.x="N",by.y="N",all.x=T)[,-2]
  tab2[is.na(tab2[,2]),2] <- 0
  # end: Change Point Model
  
  # quantile calculation
  if(is.na(rise.prob) & is.na(set.prob)) {
    rise.prob <- as.numeric(round(as.numeric(quantile(tab1[tab1[,2]!=0,2],probs=quantile,na.rm=TRUE)), digits=5))
    set.prob  <- as.numeric(round(as.numeric(quantile(tab2[tab2[,2]!=0,2],probs=quantile,na.rm=TRUE)), digits=5))
  }
  
  
  riseProb <- ifelse(tab1[,2]>=rise.prob, NA, TRUE)
  setProb  <- ifelse(tab2[,2]>=set.prob, NA, TRUE)
  
  tmp02 <- rbind(data.frame(time = tw[tw[,2]==1,1], prob = tab1[,2], cut = riseProb), 
                 data.frame(time = tw[tw[,2]==2,1], prob = tab2[,2], cut = setProb))[order(c(sr, ss)),]
  tmp02 <- cbind(tmp02, NA)
  
  s <- 1
  for(i in 2:nrow(tmp02)) {
    if(is.na(tmp02[i-1, 3]) & !is.na(tmp02[i, 3])) {
      s <- s+1
      tmp02[i, 4] <- s
    }
    if(!is.na(tmp02[i-1, 3]) & !is.na(tmp02[i, 3])) tmp02[i, 4] <- s
  }
  
  ind01 <- tapply(as.numeric(tmp02[,1]), tmp02[,4], FUN = function(x) ((x[length(x)]-x[1])/60/60/24)>days)
  ind02 <- as.numeric(names(ind01)[ind01])
  
  tmp02[,4] <- ifelse(tmp02[,4]%in%ind02, tmp02[,4], NA)
  
  s <- 1
  for(i in ind02) {
    tmp02[!is.na(tmp02[,4]) & tmp02[,4]==i, 4] <- s
    s <- s+1
  }
  
  
  t02 <- schedule(tab$tFirst, tab$tSecond, tmp02[1:nrow(tab),4])
  
  arr <- tmp02[!is.na(tmp02[,4]) & !duplicated(tmp02[,4]),]
  dep <- tmp02[!is.na(tmp02[,4]) & !duplicated(tmp02[,4], fromLast = T),]
  
  t02$P.start <- arr$prob
  t02$P.end   <- dep$prob
  
  t02$Days <- apply(t02, 1, function(x) round(as.numeric(difftime(x[3], x[2], units = "days")), 1))
  ds <- t02[,c(1:3, 6, 4:6,4)]
  
  out <- list(riseProb = tab1[,2], setProb = tab2[,2], rise.prob = rise.prob, set.prob = set.prob, 
              site = ifelse(!is.na(tmp02[,4]), tmp02[,4], 0)[1:nrow(tab)],
              migTable = ds)
  
  
  if(plot){
    def.par <- par(no.readonly = TRUE)
    nf <- layout(matrix(c(4,1,2,3),nrow=4,byrow=T),heights=c(0.5,1,0.5,0.5))
    
    par(mar=c(2,4.5,2,5),cex.lab=1.5,cex.axis=1.5,bty="o")
    plot(sr, rise, type="o",cex=0.2,col="firebrick",ylab="Sunrise (red)", 
         xlim = range(sr),xaxt="n")
    par(new=T)
    plot(ss, set, type="o",cex=0.2,col="cornflowerblue",xaxt="n",yaxt="n",xlab="",
         ylab="",xlim=range(ss))
    axis(4)
    mtext("Sunset (blue)",4,line=2.7,cex=1)
    axis(1,at=seq(min(ss),max(ss), by=(10*24*60*60)),labels=F)
    axis(1,at=seq(min(ss),max(ss), by=(30*24*60*60)),lwd.ticks=2,
         labels=trunc(seq(min(ss),max(ss), by=(30*24*60*60)), "days"),cex.axis=1)
    
    par(mar=c(1.5,4.5,0.8,5),bty="n")
    plot(sr, tab1[,2], type = "h", lwd = 4, col = "firebrick", ylab = "", xaxt = "n",
         xlim= range(sr) ,ylim = c(0, max(na.omit(c(tab1[,2],tab2[,2])))))
    if(is.numeric(rise.prob)) abline(h = rise.prob, lty=2, lwd = 1.5)
    
    opar <- par(mar=c(1.5,4.5,0.8,5),bty="n")
    plot(ss, tab2[,2], type = "h", lwd = 4, col = "cornflowerblue", ylab = "", xaxt = "n",
         xlim= range(ss) ,ylim = c(0, max(na.omit(c(tab1[,2],tab2[,2])))))
    if(is.numeric(set.prob)) abline(h = set.prob, lty=2, lwd = 1.5)
    
    mtext("Probability of change", side=2, at = max(na.omit(c(tab1[,2],tab2[,2]))), line=3)
    
    par(mar=c(1,4.5,1,5),bty="o")
    mig <- out$site
    mig[mig>0] <- 1
    plot(tab[,1] + (tab[,2] - tab[,1])/2, 
         ifelse(out$site>0, 1, 0), type = "l", yaxt = "n", ylab = NA, ylim=c(0,1.5))
    rect(tw[out$site > 0 & !duplicated(out$site), 1], 1.1, 
         tw[out$site > 0 & !duplicated(out$site, fromLast = T), 1], 1.4, lwd = 0, 
         col = "grey")
    par(def.par)
  }
  
  
  if(summary){i.sum.Cl(out)}
  
  return(out)
}



#' Function to merge sites
#'
#' The \code{\link{changeLight}} functions provides a vector grouping the twilight times
#' into stationary (>0) and movement (0) periods. This function was written to enable the user
#' to merge sites based on the distance between consequtive sites. NOTE: The function requires
#' position estimate and desicison on whether sites should be merged will be made based on
#' the defined \code{distance}, the \code{cutoff} values and the provided positions. The analysis
#' is this dependent on the accuracy of the position estiamtes and should be applied to positons that
#' were estimated using a sensible sun elevation angle.
#'
#' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
#' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
#' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
#' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type}
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0. This \code{vector} will be used as the initial state.
#' @param degElevation the sun elevation angle (in degrees) that defines twilight (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst} or \code{nrow(x)}.
#' @param distThreshold a \code{numerical} value defining the threshold of the distance under 
#' which consequtive sites should be merged (in km).
#' @param alpha mean and standard variation for position optimization process.
#' @param plot \code{logical}, if TRUE a plot comparing the inital and the final site selection.
#' @return A \code{vector} with the merged site numbers
#' @author Simeon Lisovski
#'
#' @export mergeSites
#' @importFrom fields rdist.earth
#' @importFrom graphics abline axis lines mtext par plot points rect 
#' @importFrom stats optim dnorm
mergeSites <- function(tFirst, tSecond, type, twl, site, degElevation, distThreshold = 250, 
                       alpha = c(0, 15), plot = TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), 
                                                  FUN = function(x) any(class(x) != "name"))])
  site0 <- site
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type = c(tab$type, ifelse(tab$type == 1, 2, 1)))
  tw <- tw[!duplicated(tw$datetime), ]
  tw <- tw[order(tw[, 1]), ]
  tw <- tw[1:nrow(tab),]
  
  crds0 <- coord(tab, degElevation = degElevation, note = F)
  tw$lon <- crds0[,1]
  tw$lat <- crds0[,2]
  
  lonlim <- range(crds0[, 1], na.rm = T)
  lon.seq <- seq(lonlim[1] - 1, lonlim[2] + 1, by = 1)
  latlim <- range(crds0[, 2], na.rm = T)
  lat.seq <- seq(latlim[1] - 1, latlim[2] + 1, by = 1)
  mod <- function(x) {
    loglik <- function(crds) {
      t.tw <- twilight(x[, 1], lon = crds[1], lat = crds[2], 
                       rise = ifelse(x[, 2] == 1, TRUE, FALSE), zenith = 90 - 
                         degElevation)
      diff <- as.numeric(difftime(x[, 1], t.tw, units = "mins"))
      -sum(dnorm(diff, alpha[1], alpha[2], log = T), na.rm = T)
    }
    fit0 <- optim(cbind(median(x[,3], na.rm = T), median(x[,4], na.rm = T)), loglik, lower = cbind(lonlim[1], 
                                                                                                   latlim[1]), upper = cbind(lonlim[2], latlim[2]), 
                  method = "L-BFGS-B", hessian = T)
    fit  <- optim(cbind(fit0$par[1], fit0$par[2]), loglik, lower = cbind(lonlim[1], 
                                                                         latlim[1]), upper = cbind(lonlim[2], latlim[2]), 
                  method = "L-BFGS-B", hessian = T)
    fisher_info <- solve(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    lon.lower <- c(fit$par[1] - 1.96 * prop_sigma)[1]
    lat.lower <- c(fit$par[2] - 1.96 * prop_sigma)[4]
    lon.upper <- c(fit$par[1] + 1.96 * prop_sigma)[1]
    lat.upper <- c(fit$par[2] + 1.96 * prop_sigma)[4]
    matrix(c(fit$par[1], fit$par[2], lon.lower, lat.lower, 
             lon.upper, lat.upper), ncol = 6)
  }
  out <- data.frame(site = unique(site[site != 0]), t(sapply(split(tw[site != 0, ], f = site[site != 0]), mod)))
  rep = TRUE
  ite = 1
  repeat {
    for (i in 1:(length(site[site != 0 & !duplicated(site)]) - 1)) {
      dist0 <- rdist.earth(out[i:(i + 1), 2:3])[2, 1]
      if (dist0 <= distThreshold) 
        break
    }
    if (i < (length(site[site != 0 & !duplicated(site)]) - 1)) {
      site[(which(site == i)[1]):(which(site == (i + 1))[sum(site == (i + 1))])] <- i
      site[which(site > i)] <- site[which(site > i)] - 1
    } else rep = FALSE
    out <- data.frame(site = unique(site[site != 0]), t(sapply(split(tw[site !=  0, ], f = site[site != 0]), mod)))
    if (!rep) 
      break
    else ite <- ite + 1
  }
  
  
  if (plot) {
    hours0 <- as.numeric(format(tw[, 1], "%H")) + as.numeric(format(tw[, 1], "%M"))/60
    crd0 <- out[match(site, out$site), 2:3]
    crd0[!is.na(crd0[,1]),] <- crds0[!is.na(crd0[,1]),]
    hours1 <- twilight(tw[, 1], rise = ifelse(tw[, 2] == 1, TRUE, FALSE), zenith = 90 - degElevation, 
                       lon = out[match(site, out$site), 2], lat = out[match(site, out$site), 3])
    hours1 <- as.numeric(format(hours1, "%H")) + as.numeric(format(hours1, "%M"))/60
    hours2 <- twilight(tw[, 1], rise = ifelse(tw[, 2] == 1, TRUE, FALSE), zenith = 90 - degElevation, 
                       lon = out[match(site, out$site), 4], lat = out[match(site, out$site), 3])
    
    hours2 <- as.numeric(format(hours2, "%H")) + as.numeric(format(hours2,"%M"))/60
    hours3 <- twilight(tw[, 1], rise = ifelse(tw[, 2] == 1, TRUE, FALSE), zenith = 90 - degElevation, 
                       lon = out[match(site, out$site), 6], lat = out[match(site, out$site), 3])
    hours3 <- as.numeric(format(hours3, "%H")) + as.numeric(format(hours3,"%M"))/60
    
    
    for (t in 1:2) {
      cor <- rep(NA, 24)
      for (i in 0:23) {
        cor[i + 1] <- max(abs((c(hours0[tw$type == t][1], 
                                 hours0[tw$type == t]) + i)%%24 - 
                                (c(hours0[tw$type == t], 
                                   hours0[tw$type == t][length(hours0)]) +i)%%24), 
                          na.rm = T)
      }
      hours0[tw$type == t] <- (hours0[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
      hours1[tw$type == t] <- (hours1[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
      hours2[tw$type == t] <- (hours2[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
      hours3[tw$type == t] <- (hours3[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
    }
    opar <- par(mfrow = c(5, 1), oma = c(5, 0, 0, 0), mar = c(1.5,5, 1, 1))
    mig1 <- site0
    mig1[mig1 > 0] <- 1
    mig2 <- site
    mig2[mig2 > 0] <- 1
    plot(tw[, 1], ifelse(mig2 > 0, 1, 0), type = "l", yaxt = "n", 
         ylab = NA, ylim = c(0, 1.5), col = "firebrick", lwd = 2, 
         xaxt = "n")
    lines(tw[, 1], ifelse(mig1 > 0, 1, 0), type = "l", lty = 2)
    rect(tw[site > 0 & !duplicated(site), 1], 1.1, 
         tw[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, lwd = 0, 
         col = "grey")
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[tw[, 2] == 1, 1], hours1[tw[, 2] == 1], type = "l", 
         lwd = 2, col = "firebrick", ylab = "Sunrise (red)", 
         xlim = range(tw[, 1]), ylim = range(hours0[tw[, 2] == 1]), xaxt = "n")
    lines(tw[tw[, 2] == 1, 1], hours2[tw[, 2] == 1], type = "l", 
          lwd = 1, lty = 2)
    lines(tw[tw[, 2] == 1, 1], hours3[tw[, 2] == 1], type = "l", 
          lwd = 1, lty = 2)
    points(tw[tw[, 2] == 1, 1], hours0[tw[, 2] == 1], cex = 0.5, 
           pch = 21, col = "black", bg = "firebrick", lwd = 0.5)
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[tw[, 2] == 2, 1], hours1[tw[, 2] == 2], type = "l", 
         lwd = 2, col = "cornflowerblue", ylab = "Sunset (blue)", 
         xlim = range(tw[, 1]), ylim = range(hours0[tw[, 2] == 
                                                      2]), xaxt = "n")
    lines(tw[tw[, 2] == 2, 1], hours2[tw[, 2] == 2], type = "l", 
          lwd = 1, lty = 2)
    lines(tw[tw[, 2] == 2, 1], hours3[tw[, 2] == 2], type = "l", 
          lwd = 1, lty = 2)
    points(tw[tw[, 2] == 2, 1], hours0[tw[, 2] == 2], cex = 0.5, 
           pch = 21, col = "black", bg = "cornflowerblue", lwd = 0.5)
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[, 1], crds0[, 1], type = "o", pch = 16, cex = 0.5, 
         xaxt = "n", ylab = "Longitude", cex.lab = 1.7, xlab = "")
    abline(v = c(tw[site0 > 0 & !duplicated(site0), 1], 
                 tw[site0 > 0 & !duplicated(site0, fromLast = T), 1]), lty = 2)
    abline(v = c(tw[site > 0 & !duplicated(site), 1], 
                 tw[site > 0 & !duplicated(site, fromLast = T), 1]), lwd = 1.5, 
           col = "firebrick")
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[, 1], crds0[, 2], type = "o", pch = 16, cex = 0.5, 
         xaxt = "n", ylab = "Latitude", cex.lab = 1.7, xlab = "")
    abline(v = c(tw[site0 > 0 & !duplicated(site0), 1], 
                 tw[site0 > 0 & !duplicated(site0, fromLast = T), 1]), lty = 2)
    abline(v = c(tw[site > 0 & !duplicated(site), 1], 
                 tw[site > 0 & !duplicated(site, fromLast = T), 1]), lwd = 1.5, 
           col = "firebrick")
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = format(seq(tw[1, 1], tw[nrow(tw), 1], length = 10), "%d-%b"))
    mtext("Date", 1, outer = T, line = 1.6, cex = 1.2)
    par(opar)
  }
  names(out) <- c("site", "Lon", "Lat", "Lon.upper", "Lat.upper", 
                  "Lon.lower", "Lat.lower")
  list(site = site, summary = out)
}






##' Filter for unrealistic positions within a track based on distance
##'
##' The filter identifies unrealistic positions. The maximal distance per
##' hour/day can be set corresponding to the particular species.
##'
##' Note that this type of filter significantly depends on the calibration
##' (\code{degElevation}). Especially during equinox periods. In contrast, the
##' (\code{loessFilter}) is independent from positions (uses twilight times) 
##' and therefore superior.
##'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
##' twilight")
##' @param distance the maximal distance in km per \code{units}. Distances above
##' will be considered as unrealistic.
##' @param units the time unite corresponding to the distance. Default is
##' "hour", alternative option is "day".
##' @return Logical \code{vector}. TRUE means the particular position passed the filter.
##' @author Simeon Lisovski, Fraenzi Korner-Nievergelt
##' @importFrom stats loess
##' @export distanceFilter
distanceFilter <- function(tFirst, tSecond, type, twl, degElevation = -6, distance, units = "hour") {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  if(units=="days") units <- distance/24
  
  tFirst <- as.POSIXct(tab$tFirst, tz = "GMT")
  tSecond<- as.POSIXct(tab$tSecond, tz = "GMT")
  type <- tab$type
  tSunTransit <- tFirst + (tSecond-tFirst)/2
  crds  <- coord(tab, degElevation, note=FALSE)
  
  crds[is.na(crds[,2]),2] <- 999
  
  difft <- as.numeric(difftime(tSunTransit[-length(tSunTransit)],tSunTransit[-1],units="hours"))
  diffs <- abs(as.numeric(i.loxodrom.dist(crds[-nrow(crds), 1], crds[-nrow(crds), 2], crds[-1, 1], crds[-1, 2]))/difft)
  
  index <- rep(TRUE,length(tFirst))
  index[diffs>distance] <- FALSE
  index[crds[,2]==999] <- TRUE
  
  
  cat(paste("Note: ",length(index[!index])," of ",length(index[crds[,2]!=999])," positions were filtered (",floor((length(index[!index])*100)/length(index[crds[,2]!=999]))," %)",sep=""))
  index
}


#' Transformation of *.gle files
#'
#' Function to transform *.gle files derived by the software GeoLocator for
#' further analyses in \bold{\code{GeoLight}}.
#'
#' The *.gle file derived by the software "GeoLocator" (Swiss Ornithological
#' Institute) is a table with interpolated light intensities over time gathered
#' from the *.glf file. Furthermore a column defines wether the light intensity
#' passes the defined light intensity threshold in the morning (sunrise) or in
#' the evening (sunset). This information is used in \code{gleTrans()}, to
#' create a table with two subsequent twilight events (\emph{tFirst, tSecond})
#' and \emph{type} defining wether \emph{tFirst} refers to sunrise (1) or
#' sunset (2). Date and time information will be transferred into a
#' \code{GeoLight} appropriate format (see: \code{\link{as.POSIXct}}).
#'
#' @param file the full patch and filename with suffix of the *.gle file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{glfTrans}}
#' @importFrom utils read.table
#' @export gleTrans
gleTrans <- function(file) {
  
  
  gle1 <- read.table(file,sep="\t",skip=16,col.names=c("date","light","1","2","3","4","5","6","7","8","type")) # read file
  gle1 <- subset(gle1,gle1$type>0,select=c("date","type"))
  
  # Date transformation
  year   <- as.numeric(substring(gle1$date,7,10))
  month  <- as.numeric(substring(gle1$date,4,5))
  day    <- as.numeric(substring(gle1$date,1,2))
  hour   <- as.numeric(substring(gle1$date,12,13))
  min    <- as.numeric(substring(gle1$date,15,16))
  gmt.date <- paste(year,"-", month,"-",day," ",hour,":",min,":",0,sep="")
  gmt.date <- as.POSIXct(strptime(gmt.date, "%Y-%m-%d %H:%M:%S"), "UTC")
  
  gle <- data.frame(date=gmt.date,type=gle1$type)
  
  
  opt <- data.frame(tFirst=as.POSIXct("1900-01-01 01:01","UTC"),tSecond=as.POSIXct("1900-01-01 01:01","UTC"),type=0)
  
  row <- 1
  for (i in 1:(length(gmt.date)-1))
  {
    if (abs(as.numeric(difftime(gle$date[i],gle$date[i+1],units="hours")))< 18 & gle$date[i] != gle$date[i+1])
    {
      opt[row,1] <- gle$date[i]
      opt[row,2] <- gle$date[i+1]
      opt[row,3] <- gle$type[i]
      
      row <- row+1
    }
  }
  
  return(opt)
}

#' Transformation of *.glf files
#'
#' Transform *.glf files derived by the software GeoLocator for further
#' analyses in \bold{\code{GeoLight}}.
#'
#' The *.glf files produced by the software GeoLocator (Swiss Ornithological
#' Institute) is a table with light intensity measurements over time.
#' \code{glfTrans} produces a table with these measurements and transfer the
#' data and time information into the format required by \bold{\code{GeoLight}}
#' format (see: \code{\link{as.POSIXct}}).
#'
#' @param file the full patch and filename with suffix of the *.glf file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{gleTrans}}; \code{\link{luxTrans}} for transforming
#' *.lux files produced by \emph{Migrate Technology Ltd}
#' @importFrom utils read.table
#' @export glfTrans
glfTrans <- function(file="/path/file.glf") {
  
  
  glf1 <- read.table(file,sep="\t",skip=6,col.names=c("datetime","light","1","2","3")) # read file
  
  # Date transformation
  year   <- as.numeric(substring(glf1$datetime,7,10))
  month  <- as.numeric(substring(glf1$datetime,4,5))
  day    <- as.numeric(substring(glf1$datetime,1,2))
  hour   <- as.numeric(substring(glf1$datetime,12,13))
  min    <- as.numeric(substring(glf1$datetime,15,16))
  gmt.date <- paste(year,"-", month,"-",day," ",hour,":",min,":",0,sep="")
  gmt.date <- as.POSIXct(strptime(gmt.date, "%Y-%m-%d %H:%M:%S"), "UTC")
  
  glf <- data.frame(datetime=gmt.date,light=glf1$light)
  
  return(glf)
}

#' Hill-Ekstrom calibration
#'
#' Hill-Ekstrom calibration for one or multiple stationary periods.
#'
#' The \emph{Hill-Ekstrom calibration} has been suggested by Hill & Braun
#' (2001) and Ekstrom (2004), and allows for calibrating data during stationary
#' periods at unknown latitudinal positions. The Hill-Ekstrom calibration bases
#' on an increasing error range in latitudes with an increasing mismatch
#' between light level threshold and the used sun angle. This error is strongly
#' amplified with proximity to the equinox times due to decreasing slope of day
#' length variation with latitude. Furthermore, the sign of the error switches
#' at the equinox, i.e. latitude is overestimated before the equinox and
#' underestimated after the equinox (or vice versa depending on autumnal/vernal
#' equinox, hemisphere, and sign of the mismatch between light level threshold
#' and sun angle). When calculating the positions of a stationary period, the
#' variance in latitude is minimal if the sun elevation angle fits to the
#' defined light level threshold. Moreover, the accuracy of positions increases
#' with decreasing variance in latitudes. \bold{However, the method is only
#' applicable for stationary periods and under stable shading intensities}. The
#' plot produced by the function may help to judge visually if the calculated
#' sun elevation angles are realistic (e.g. site 2 in the example below) or not
#' (e.g. site 3 in the example below.
#'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param site a \code{numerical vector} assigning each row to a particular
##' period. Stationary periods in numerical order and values >0,
##' migration/movement periods 0
##' @param start.angle a single sun elevation angle. The combined process of
##' checking for minimal variance in resulting latitude, which is the initial
##' value for the sun elevation angle in the iterative process of identifying
##' the latitudes with the least variance
##' @param distanceFilter logical, if TRUE the \code{\link{distanceFilter}} will
##' be used to filter unrealistic positions
##' @param distance if \code{distanceFilter} is set \code{TRUE} a threshold
##' distance in km has to be set (see: \code{\link{distanceFilter}})
##' @param plot logical, if TRUE the function will give a plot with all relevant
##' information
##' @return A \code{vector} of sun elevation angles corresponding to the
##' Hill-Ekstrom calibration for each defined period.
##' @author Simeon Lisovski
##' @references Ekstrom, P.A. (2004) An advance in geolocation by light.
##' \emph{Memoirs of the National Institute of Polar Research}, Special Issue,
##' \bold{58}, 210-226.
##'
##' Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step:
##' Latitude. In: \emph{Electronic Tagging and Tracking in Marine Fisheries}
##' (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The
##' Netherlands.
##'
##' Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F.,
##' Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
##' precision affected by environmental factors. \emph{Methods in Ecology and
##' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
##' @examples
##'
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' residency <- with(hoopoe2, changeLight(tFirst,tSecond,type, rise.prob=0.1, 
##'                   set.prob=0.1, plot=FALSE, summary=FALSE))
##' HillEkstromCalib(hoopoe2,site = residency$site)
##'
##' @importFrom grDevices grey.colors
##' @importFrom graphics abline axis mtext legend lines par plot points text
##' @importFrom stats var na.omit
##' @export HillEkstromCalib
HillEkstromCalib <- function(tFirst, tSecond, type, twl, site, start.angle=-6, distanceFilter=FALSE, distance, plot=TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  tFirst <- tab$tFirst
  tSecond <- tab$tSecond
  type <- tab$type
  
  sites <- as.numeric(length(levels(as.factor(site[as.numeric(site)!=0]))))
  
  HECalib <- rep(NA,sites)
  
  
  for(j in 1:sites) {
    b <- 0
    start <- start.angle
    
    repeat{
      # backwards
      if(start-((b*0.1)-0.1) < -9) {
        HECalib[j] <- NA
        break
      }
      t0 <- var(na.omit(coord(tab[site==j,], degElevation = start-((b*0.1)-0.1), note = F)[,2]))
      t1 <- var(na.omit(coord(tab[site==j,], degElevation = start-(b*0.1), note = F)[,2]))
      t2 <- var(na.omit(coord(tab[site==j,], degElevation = start-((b*0.1)+0.1), note = F)[,2]))
      if(sum(is.na(c(t0,t1,t2)))>0) {
        HECalib[j] <- NA
        break
      }
      if(t0>t1 & t1<t2) {
        HECalib[j] <- start-(b*0.1)
        break}
      
      # forward
      if(start-((b*0.1)-0.1) > 9) {
        HECalib[j] <- NA
        break
      }
      f0 <- var(na.omit(coord(tab[site==j,], degElevation = start+((b*0.1)-0.1), note = F)[,2]))
      f1 <- var(na.omit(coord(tab[site==j,], degElevation = start+(b*0.1), note = F)[,2]))
      f2 <- var(na.omit(coord(tab[site==j,], degElevation = start+((b*0.1)+0.1), note = F)[,2]))
      if(sum(is.na(c(f0,f1,f2)))>0) {
        HECalib[j] <- NA
        break
      }
      if(f0>f1 & f1<f2) {
        HECalib[j] <- start+(b*0.1)
        break}
      
      b <- b+1
    }
  }
  
  if(plot) {

      opar <- par(mfrow = c(sites, 2), oma = c(3, 5, 6, 2), mar = c(2, 2, 2, 2))
      for(j in 1:sites){
        
        if(is.na(HECalib[j])){
          
          plot(0,0,cex=0,pch=20,col="white",ylab="",xlab="",xaxt="n",yaxt="n", bty = "n")
          text(0,0,"NA",cex=2)
          plot(0,0,cex=0,pch=20,col="white",ylab="",xlab="",xaxt="n",yaxt="n", bty="n")
          
        } else {
          
          angles <- c(seq(HECalib[j]-2,HECalib[j]+2,0.2))
          
          latM <- matrix(ncol=length(angles),nrow=length(tFirst[site==j]))
          
          for(i in 1:ncol(latM)){
            latM[,i] <- coord(tab[site==j,], degElevation = c(angles[i]),note=F)[,2]
          }
          
          latT <- latM
          var1 <- rep(NA,ncol(latT))
          n1   <- rep(NA,ncol(latT))
          min  <- rep(NA,ncol(latT))
          max  <- rep(NA,ncol(latT))
          for(t in 1:length(var1)){
            var1[t] <- var(na.omit(latT[,t]))
            n1[t]   <- length(na.omit(latT[,t]))
            min[t]  <- if(length(na.omit(latT[,t]))<=1) NA else min(na.omit(latT[,t]))
            max[t]  <- if(length(na.omit(latT[,t]))<=1) NA else max(na.omit(latT[,t]))
          }
          
          
          colors <- grey.colors(length(angles))
          plot(tFirst[site==j],latT[,1],ylim=c(min(na.omit(min)),max(na.omit(max))),type="o",cex=0.7,pch=20,col=colors[1],ylab="",xlab="")
          for(p in 2:ncol(latT)){
            lines(tFirst[site==j],latM[,p],type="o",cex=0.7,pch=20,col=colors[p])
          }
          lines(tFirst[site==j],coord(tab[site==j,], degElevation = HECalib[j], note=F)[,2],col="tomato2",type="o",lwd=2,cex=1,pch=19)
          if(j==sites) mtext("Latitude",side=2,line=3)
          if(j==sites) mtext("Date",side=1,line=2.8)
          
          plot(angles,var1,type="o",cex=1,pch=20,ylab="")
          lines(angles,var1,type="p",cex=0.5,pch=7)
          abline(v=HECalib[j],lty=2,col="red",lwd=1.5)
          par(new=T)
          plot(angles,n1,type="o",xaxt="n",yaxt="n",col="blue",pch=20,ylab="")
          points(angles,n1,type="p",cex=0.5,pch=8,col="blue")
          axis(4)
          if(j==sites) mtext("Sun elevation angles",side=1,line=2.8)
          legend("topright",c(paste(HECalib[j]," degrees",sep=""),"sample size","variance"),pch=c(-1,20,20),lty=c(2,2,2),lwd=c(1.5,0,0),col=c("red","blue","black"),bg="White")
        }
      }
      
      mtext("Hill-Ekstrom Calibration",cex=1.5,line=0.8,outer=T)
      par(opar)

  }
  
  
  return(HECalib)
}



i.JC2000 <- function(jD) {
  
  #--------------------------------------------------------------------------------------------------------
  # jD: julian Date
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  
  jC<- (jD - 2451545)/36525
  
  return(jC)
  
}

i.deg <- function(Rad) {
  
  #------------------------------------------------------------
  # Deg: 	The input angle in radian
  #------------------------------------------------------------
  
  options(digits=10)
  
  
  return(Rad * (180/pi))
  
}

i.frac <- function(In) {
  
  #------------------------------------------------------------
  # In: 	numerical Number
  #------------------------------------------------------------
  
  options(digits=10)
  
  
  return(In - floor(In))
  
}

i.get.outliers<-function(residuals, k=3) {
  x <- residuals
  # x is a vector of residuals
  # k is a measure of how many interquartile ranges to take before saying that point is an outlier
  # it looks like 3 is a good preset for k
  QR<-quantile(x, probs = c(0.25, 0.75))
  IQR<-QR[2]-QR[1]
  Lower.band<-QR[1]-(k*IQR)
  Upper.Band<-QR[2]+(k*IQR)
  delete<-which(x<Lower.band |  x>Upper.Band)
  return(as.vector(delete))
}

i.julianDate <- function(year,month,day,hour,min) {
  
  #--------------------------------------------------------------------------------------------------------
  # Year:		Year as numeric e.g. 2010
  # Month:	Month as numeric e.g.   1
  # Day:		Day as numeric e.g.		1
  # Hour:		Hour as numeric e.g.   12
  # Min:		Minunte as numeric e.g. 0
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=15)
  
  fracOfDay	<- hour/24 + min/1440
  
  # Julian date (JD)
  # ------------------------------------
  
  index1 <- month <= 2
  if(sum(index1) > 0)
  {
    year[index1]  <- year[index1] -1
    month[index1] <- month[index1] +12
  }
  
  index2 <- (year*10000)+(month*100)+day <= 15821004
  
  JD <- numeric(length(year))
  if(sum(index2)>0)
  {
    JD[index2] <- floor(365.25*(year[index2]+4716)) + floor(30.6001*(month[index2]+1)) + day[index2] + fracOfDay[index2] - 1524.5
  }
  index3 <- year*10000+month*100+day >= 15821015
  if (sum(index3)>0)
  {
    a <- floor(year/100)
    b <- 2 - a + floor(a/4)
    
    JD[index3] <- floor(365.25*(year[index3]+4716)) + floor(30.6001*(month[index3]+1)) + day[index3] + fracOfDay[index3] + b[index3] - 1524.5
  }
  JD[!index2&!index3] <- 1
  
  return(JD)
  
}

i.loxodrom.dist <- function(x1, y1, x2, y2, epsilon=0.0001){
  dis<-numeric(length(x1))
  rerde<-6368
  deltax<-abs(x2*pi/180-x1*pi/180)
  deltay<-abs(y2*pi/180-y1*pi/180)
  tga<-deltax/(log(tan(pi/4+y2*pi/360))-log(tan(pi/4+y1*pi/360)))
  
  dis[abs(x1-x2)<epsilon&abs(y1-y2)<epsilon]<-0
  dis[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]<-abs(cos(y1[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]*pi/180)*deltax[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)])
  dis[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]<-abs(deltay[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]/cos((pi-atan(tga[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]))))
  dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]<--deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]))
  dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]<-abs(deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)])))
  dis[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)])))
  dis[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)])))
  dis*rerde
}

i.preSelection <- function(datetime, light, LightThreshold){
  
  dt <- cut(datetime,"1 hour")
  st <- as.POSIXct(levels(dt),"UTC")
  
  raw <- data.frame(datetime=dt,light=light)
  
  h  <- tapply(light,dt,max)
  df1 <- data.frame(datetime=st+(30*60),light=as.numeric(h))
  
  smooth <- i.twilightEvents(df1[,1], df1[,2], LightThreshold)
  smooth <- data.frame(id=1:nrow(smooth),smooth)
  raw    <- i.twilightEvents(datetime, light, LightThreshold)
  raw <- data.frame(id=1:nrow(raw),raw)
  
  ind2 <- rep(NA,nrow(smooth))
  for(i in 1:nrow(smooth)){
    tmp <- subset(raw,datetime>=(smooth[i,2]-(90*60)) & datetime<=(smooth[i,2]+(90*60)))
    
    if(smooth[i,3]==1) ind3 <- tmp$id[which.min(tmp[,2])]
    if(smooth[i,3]==2) ind3 <- tmp$id[which.max(tmp[,2])]
    ind2[i] <- ind3
  }
  
  
  res <- data.frame(raw,mod=1)
  res$mod[ind2] <- 0
  
  return(res)
}

i.rad <- function(Deg) {
  
  #------------------------------------------------------------
  # Deg: 	The input angle in degrees
  #------------------------------------------------------------
  
  options(digits=10)
  
  
  return(Deg * (pi/180))
  
}

i.radDeclination <- function(radEclipticalLength,radObliquity) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane.
  # EadObliquity:
  #-------------------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  dec <- asin(sin(radObliquity)*sin(radEclipticalLength))
  
  return(dec)
  
}

i.radEclipticLongitude <- function(jC) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00
  #-------------------------------------------------------------------------------------------------------------------
  
  
  
  options(digits=10)
  
  radMeanAnomaly <- 2*pi*i.frac(0.993133 + 99.997361*jC)
  EclipticLon    <- 2*pi*i.frac(0.7859452 + radMeanAnomaly/(2*pi) + (6893*sin(radMeanAnomaly) + 72*sin(2*radMeanAnomaly) + 6191.2*jC) / 1296000)
  
  return(EclipticLon)
  
}

i.radGMST <- function(jD,jD0,jC,jC0) {
  
  #--------------------------------------------------------------------------------------------------------
  # jD:  Julian Date with Hour and Minute
  # jD0: Julan Date at t 0 UT
  # jC:  Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
  # jC0: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00) at t 0 UT
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  
  UT  <- 86400 * (jD-jD0)
  st0 <- 24110.54841 + 8640184.812866*jC0 + 1.0027379093*UT + (0.093104 - 0.0000062*jC0)*jC0*jC0
  gmst<- (((2*pi)/86400)*(st0%%86400))
  
  return(gmst)
  
}

i.radObliquity <- function(jC) {
  
  #--------------------------------------------------------------------------------------------------------
  # jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  
  degObliquity <- 23.43929111 - (46.8150 + (0.00059 - 0.001813 *jC)*jC) *jC/3600
  radObliquity <- i.rad(degObliquity)
  
  return(radObliquity)
  
}

i.radRightAscension <- function(RadEclipticalLength,RadObliquity) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane.
  # EadObliquity:
  #-------------------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  index1 <- (cos(RadEclipticalLength) < 0)
  
  res <- numeric(length(RadEclipticalLength))
  if (sum(index1)>0)
  {
    res[index1] <- (atan((cos(RadObliquity[index1])*sin(RadEclipticalLength[index1]))/cos(RadEclipticalLength[index1])) + pi)
  }
  
  index2 <- (cos(RadEclipticalLength) >= 0)
  if (sum(index2)>0)
  {
    res[index2] <- (atan((cos(RadObliquity[index2])*sin(RadEclipticalLength[index2]))/cos(RadEclipticalLength[index2])))
  }
  
  
  return(res)
  
}

i.setToRange <- function(Start,Stop,Angle) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # Start:	 Minimal value of the range in degrees
  # Stop: 	 Maximal value of the range in degrees
  # Angle:	 The angle that should be fit to the range
  #-------------------------------------------------------------------------------------------------------------------
  
  options(digits=15)
  
  angle <- Angle
  range <- Stop - Start
  
  
  
  index1 <- angle >= Stop
  if (sum(index1, na.rm = T)>0) angle[index1] <- angle[index1] - (floor((angle[index1]-Stop)/range)+1)*range
  
  index2 <- angle < Start
  if(sum(index2, na.rm = T)>0) angle[index2] <- angle[index2]  + ceiling(abs(angle[index2] -Start)/range)*range
  
  return(angle)
  
}

i.sum.Cl <- function(object) {
  
  if(all(names(object)%in%c("riseProb","setProb","rise.prob","set.prob","site","migTable"))){
    cat("\n")
    cat("Probability threshold(s):")
    cat(rep("\n",2))
    if(!is.na(object$rise.prob)) cat(paste("	Sunrise: ",object$rise.prob))
    if(!is.na(object$set.prob)) cat(paste("	Sunset: ",object$set.prob))
    cat(rep("\n",3))
    cat("Migration schedule table:")
    cat(rep("\n",2))
    
    print(object$migTable,quote=FALSE)
  } else {
    cat("Error: List must be the output list of the changeLight function.")
  }
}


#' Filter to remove noise in light intensity measurements during the night
#'
#' The filter identifies and removes light intensities oczillating around the
#' baseline or few light intensities resulting in a short light peak during the
#' night. Such noise during the night will increase the calculated twilight
#' events using the function \code{\link{twilightCalc}} and therewith the
#' manual work to remove these false twilight events.
#'
#' The filter searches for light levels above the baseline and compares the
#' prior and posterior levels. If these values are below the threshold the
#' particular light level will be reduced to the baseline. A few (usually two)
#' iterations might be enough to remove most noise during the night (however,
#' not if such noise occurs at the begining or at the end were not enough prior
#' or posterior values are available).
#'
#' @param light \code{numerical} value of the light intensity (usually
#' arbitrary units).
#' @param baseline the light intensity baseline (no light). If \code{Default},
#' it will be calculated as the most frequent value below the mean light
#' intensities.
#' @param iter a \code{numerical} value, specifying how many iterations should
#' be computed (see details).
#' @return numerical \code{vector} with the new light levels. Same length as
#' the initial light vector.
#' @author Simeon Lisovski
#' @examples
#'
#' night <- rep(0,50); night[runif(4,0,50)] <- 10; night[runif(4,0,50)] <- -5
#' nightday <- c(night,rep(30,50))
#' plot(nightday,type="l",ylim=c(-5,30),ylab="light level",xlab="time (time)")
#' light2 <- lightFilter(nightday, baseline=0, iter=4)
#' lines(light2,col="red")
#' legend("bottomright",c("before","after"),lty=c(1,1),col=c("black","red"),bty="n")
#'
#' @export lightFilter
lightFilter <- function(light, baseline=NULL, iter=2){
  
  r <- as.data.frame(table(light))
  r[,1] <- as.character(r[,1])
  nr <- as.numeric(which.max(r$Freq[as.numeric(r[,1])<mean(light)]))
  LightThreshold <- ifelse(is.null(baseline),as.numeric(r[nr,1]),baseline)
  
  light[light<LightThreshold] <- LightThreshold
  
  
  index <- which(light<mean(light) & light!= LightThreshold)
  
  rep   <- rep(FALSE,length(light))
  
  for(i in 1:iter){
    
    for(i in index[index>5 & index<(length(light)-5)]){
      
      back=FALSE
      if(any(light[seq(i-5,i)]==LightThreshold)) back <- TRUE
      forw=FALSE
      if(any(light[seq(i,i+5)]==LightThreshold)) forw <- TRUE
      
      if(back & forw) rep[i] <- TRUE
      
    }
    
    light[rep] <- LightThreshold
    
  }
  
  light
}


##' Filter to remove outliers in defined twilight times based on smoother function
##'
##' This filter defines outliers based on residuals from a local polynomial
##' regression fitting provcess (\code{\link{loess}}).
##'
##'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param k a measure of how many interquartile ranges to take before saying
##' that a particular twilight event is an outlier
##' @param plot codelogical, if TRUE a plot indicating the filtered times will
##' be produced.
##' @return Logical \code{vector} matching positions that pass the filter.
##' @author Simeon Lisovski & Eldar Rakhimberdiev
##' @importFrom graphics axis mtext legend lines par plot points
##' @importFrom stats loess predict residuals
##' @export loessFilter
loessFilter <- function(tFirst, tSecond, type, twl, k = 3, plot = TRUE){
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type = c(tab$type, ifelse(tab$type == 1, 2, 1)))
  tw <- tw[!duplicated(tw$datetime),]
  tw <- tw[order(tw[,1]),]
  
  hours <- as.numeric(format(tw[,1],"%H"))+as.numeric(format(tw[,1],"%M"))/60
  
  for(t in 1:2){
    cor <- rep(NA, 24)
    for(i in 0:23){
      cor[i+1] <- max(abs((c(hours[tw$type==t][1],hours[tw$type==t])+i)%%24 -
                            (c(hours[tw$type==t],hours[tw$type==t][length(hours)])+i)%%24),na.rm=T)
    }
    hours[tw$type==t] <- (hours[tw$type==t] + (which.min(round(cor,2)))-1)%%24
  }
  
  dawn <- data.frame(id=1:sum(tw$type==1),
                     datetime=tw$datetime[tw$type==1],
                     type=tw$type[tw$type==1],
                     hours = hours[tw$type==1], filter=FALSE)
  
  dusk <- data.frame(id=1:sum(tw$type==2),
                     datetime=tw$datetime[tw$type==2],
                     type=tw$type[tw$type==2],
                     hours = hours[tw$type==2], filter=FALSE)
  
  
  for(d in seq(30,k,length=5)){
    
    predict.dawn <- predict(loess(dawn$hours[!dawn$filter]~as.numeric(dawn$datetime[!dawn$filter]),span=0.1))
    predict.dusk <- predict(loess(dusk$hours[!dusk$filter]~as.numeric(dusk$datetime[!dusk$filter]),span=0.1))
    
    del.dawn <-	i.get.outliers(as.vector(residuals(loess(dawn$hours[!dawn$filter]~
                                                           as.numeric(dawn$datetime[!dawn$filter]),span=0.1))),k=d)
    del.dusk <-	i.get.outliers(as.vector(residuals(loess(dusk$hours[!dusk$filter]~
                                                           as.numeric(dusk$datetime[!dusk$filter]),span=0.1))),k=d)
    
    if(length(del.dawn)>0) dawn$filter[!dawn$filter][del.dawn] <- TRUE
    if(length(del.dusk)>0) dusk$filter[!dusk$filter][del.dusk] <- TRUE
  }
  
  if(plot){
    opar <- par(mfrow=c(2,1),mar=c(3,3,0.5,3),oma=c(2,2,0,0))
    plot(dawn$datetime[dawn$type==1],dawn$hours[dawn$type==1],pch="+",cex=0.6,xlab="",ylab="",yaxt="n")
    lines(dawn$datetime[!dawn$filter], predict(loess(dawn$hours[!dawn$filter]~as.numeric(dawn$datetime[!dawn$filter]),span=0.1)) , type="l")
    points(dawn$datetime[dawn$filter],dawn$hours[dawn$filter],col="red",pch="+",cex=1)
    axis(2,labels=F)
    mtext("Sunrise",4,line=1.2)
    
    plot(dusk$datetime[dusk$type==2],dusk$hours[dusk$type==2],pch="+",cex=0.6,xlab="",ylab="",yaxt="n")
    lines(dusk$datetime[!dusk$filter], predict(loess(dusk$hours[!dusk$filter]~as.numeric(dusk$datetime[!dusk$filter]),span=0.1)), type="l")
    points(dusk$datetime[dusk$filter],dusk$hours[dusk$filter],col="red",pch="+",cex=1)
    axis(2,labels=F)
    legend("bottomleft",c("OK","Filtered"),pch=c("+","+"),col=c("black","red"),
           bty="n",cex=0.8)
    mtext("Sunset",4,line=1.2)
    mtext("Time",1,outer=T)
    mtext("Sunrise/Sunset hours (rescaled)",2,outer=T)
    par(opar)
  }
  all <- rbind(subset(dusk,filter),subset(dawn,filter))
  
  filter <- rep(FALSE,length(tab$tFirst))
  filter[as.POSIXct(tab$tFirst, "GMT")%in%all$datetime | as.POSIXct(tab$tSecond, "GMT")%in%all$datetime] <- TRUE
  
  return(!filter)
}


#' Transformation of *.lux files
#'
#' Transform *.lux files derived from \emph{Migrate Technology Ltd} geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' The *.lux files produced by \emph{Migrate Technology Ltd} are table with
#' light intensity measurements over time. \code{luxTrans} produces a table
#' with these measurements and transfer the data and time information into the
#' format required by \bold{\code{GeoLight}} format (see:
#' \code{\link{as.POSIXct}}).
#'
#' @param file the full patch and filename with suffix of the *.lux file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{gleTrans}} for transforming *.glf files produced by the
#' software GeoLocator (\emph{Swiss Ornithological Institute})
#' @importFrom utils read.table
#' @export luxTrans
luxTrans <- function(file) {
  
  lux1 <- read.table(file,sep="\t",skip=21,col.names=c("datetime","time")) # read file
  lux <- data.frame(datetime=as.POSIXct(strptime(lux1[,1],format="%d/%m/%Y %H:%M:%S",tz="UTC")),light=lux1[,2])
  
  return(lux)
}


#' Transformation of *.lig files
#'
#' Transform *.lig files derived from geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' @param file the full patch and filename with suffix of the *.lux file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Tamara Emmenegger
#' @seealso \code{\link{gleTrans}} for transforming *.glf files produced by the
#' software GeoLocator (\emph{Swiss Ornithological Institute})
#' @importFrom utils read.csv
#' @export ligTrans
ligTrans <- function(file) {
  df <- cbind(read.csv(file,row.names=NULL)[,c(2,4)])
  colnames(df) <- c("datetime", "light")
  df$datetime <- as.POSIXct(strptime(df$datetime,format="%d/%m/%y %H:%M:%S",tz="UTC"))
  return(df)
}

#' Transformation of staroddi files
#'
#' Transform staroddi files derived from geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' @param file the full patch and filename of the staroddi file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Tamara Emmenegger
#' @importFrom utils read.table
#' @export staroddiTrans
staroddiTrans <- function(file){
  df <- cbind(read.table(file,sep="\t")[,c(2,4)])
  colnames(df) <- c("datetime", "light")
  df$datetime <- as.POSIXct(strptime(df$datetime,format="%d/%m/%y %H:%M:%S",tz="UTC"))
  return(df)
}


#' Transformation of *.trn files
#'
#' Transform *.trn files derived from geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' @param file the full patch and filename of the *.trn file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Tamara Emmenegger
#' @importFrom utils read.table
#' @export trnTrans
trnTrans<-function(file){
  data<-read.table(file,sep=",")
  tFirst <- vector("numeric",length=(length(data[,1])-1))
  tSecond <- vector("numeric",length=(length(data[,1])-1))
  type <- vector("numeric",length=(length(data[,1])-1))
  for (i in 1:(length(data[,1])-1)){
    date1<-as.Date(substr(as.character(data$V1[i]),1,8),format="%d/%m/%y")
    tFirst[i] <- as.POSIXct(paste(as.character(date1),prefix=substr(as.character(data$V1[i]),10,17)),tz="UTC")
    date2<-as.Date(substr(as.character(data$V1[i+1]),1,8),format="%d/%m/%y")
    tSecond[i] <- as.POSIXct(paste(as.character(date1),prefix=substr(as.character(data$V1[i+1]),10,17)),tz="UTC")
    if(as.character(data$V2[i])=="Sunrise") type[i] <- 1
    if(as.character(data$V2[i])=="Sunset") type[i] <- 2
  }
  output <- data.frame(tFirst=as.POSIXlt(tFirst,origin="1970-01-01",tz="UTC"),tSecond=as.POSIXlt(tSecond,origin="1970-01-01",tz="UTC"),type=type)
  return(output)
}



## ' Summary of migration/movement pattern
##'
##' Function for making a data frame summarising residency and movement pattern.
##'
##' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
##' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
##' @param site a \code{vector}, indicating the residency period of a particular
##' day (see output: \code{\link{changeLight}})
##' @return A \code{data.frame} with end and start date (yyyy-mm-dd hh:mm, UTC)
##' for each stationary period.
##' @author Simeon Lisovski
##' @examples
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' residency <- changeLight(hoopoe2, rise.prob=0.1, set.prob=0.1, plot=FALSE, summary=FALSE)
##' schedule(hoopoe2[,1], hoopoe2[,2], site = residency$site)
##' @export schedule
schedule <- function(tFirst, tSecond, site) {    
  tm <- tFirst + (tSecond - tFirst)/2
  
  arr <- tm[!is.na(site) & !duplicated(site)]
  dep <- tm[!is.na(site) & !duplicated(site, fromLast = T)]
  
  out <- data.frame(Site =  letters[1:length(arr)], Arrival = arr, Departure = dep)
  if(!is.na(site[1])) out$Arrival[1] <- NA
  if(!is.na(site[length(tFirst)])) out$Departure[nrow(out)] <- NA
  
  out
}



i.twilightEvents <- function (datetime, light, LightThreshold) 
{
  df <- data.frame(datetime, light)
  ind1 <- which((df$light[-nrow(df)] < LightThreshold & df$light[-1] > 
                   LightThreshold) | (df$light[-nrow(df)] > LightThreshold & 
                                        df$light[-1] < LightThreshold) | df$light[-nrow(df)] == 
                  LightThreshold)
  bas1 <- cbind(df[ind1, ], df[ind1 + 1, ])
  bas1 <- bas1[bas1[, 2] != bas1[, 4], ]
  x1 <- as.numeric(unclass(bas1[, 1]))
  x2 <- as.numeric(unclass(bas1[, 3]))
  y1 <- bas1[, 2]
  y2 <- bas1[, 4]
  m <- (y2 - y1)/(x2 - x1)
  b <- y2 - (m * x2)
  xnew <- (LightThreshold - b)/m
  type <- ifelse(bas1[, 2] < bas1[, 4], 1, 2)
  res <- data.frame(datetime = as.POSIXct(xnew, origin = "1970-01-01", 
                                          tz = "UTC"), type)
  return(res)
}

#' Draws sites of residency and adds a convex hull
#'
#' Draw a map (from the \code{R} Package \code{maps}) showing the defined
#' stationary sites
#'
#'
#' @param crds a \code{SpatialPoints} or \code{matrix} object, containing x
#' and y coordinates (in that order).
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0.
#' @param type either \emph{points}, or \emph{cross} to show all points for each site or only show the mean position of the site with standard deviation.
#' @param quantiles the quantile of the error bars (\emph{cross}) around the median.
#' @param hull \code{logical}, if TRUE a convex hull will be plotted around the points of each site.
#' @param map.range some possibilities to choose defined areas ("World
#' (default)", "EuroAfrica","America","AustralAsia").
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @details Standard graphical paramters like \code{pch}, \code{cex}, \code{lwd}, \code{lty} and \code{col} are implemented. 
#' The color can be specified as either a vector of colors (e.g. c("blue", "red", ...)) or as a character string indicating a color ramp (at the moment only "random" and "rainbow" is available )
#' @author Simeon Lisovski & Tamara Emmenegger
#' @examples
#' data(hoopoe2)
#'  hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
#'  hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
#' crds <- coord(hoopoe2, degElevation = -6)
#' filter <- distanceFilter(hoopoe2, distance = 30)
#' site <- changeLight(hoopoe2, rise.prob = 0.1, set.prob = 0.1, plot = FALSE, 
#'  summary = FALSE)$site
#' siteMap(crds[filter,], site[filter], xlim=c(-20,20), ylim=c(0,60), 
#'  lwd=2, pch=20, cex=0.5, main="hoopoe2")
#'
#' @importFrom maps map map.axes
#' @importFrom grDevices chull col2rgb rainbow rgb
#' @importFrom graphics legend mtext par plot points segments
#' @export siteMap
siteMap <- function(crds, site, type = "points", quantiles = c(0.25, 0.75), hull = T, 
                    map.range = c("EuroAfrica", "AustralAsia", "America", "World"), ...) {  
  args <- list(...)
  
  if(all(map.range==c("EuroAfrica", "AustralAsia", "America", "World")) & sum(names(args)%in%c("xlim", "ylim"))!=2) {
    range <- c(-180, 180, -80, 90)
  } 
  if(all(map.range=="EuroAfrica")) range <- c(-24, 55, -55, 70)
  if(all(map.range=="AustralAsia")) range <- c(60, 190, -55, 78)
  if(all(map.range=="America")) range <- c(-170, -20, -65, 78)
  if(all(map.range=="World")) range <- c(-180, 180, -75, 90)
  
  if(sum(names(args)%in%c("xlim", "ylim"))==2) range <- c(args$xlim, args$ylim)
  
  # colors for sites
  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  
  
  if(any(names(args)%in%"col")) {
    allcolors <- rainbow(60,v=0.85)[1:50]
    if(any(args$col=="rainbow")) {
      col <- allcolors[round(seq(1,50,length.out=length(unique(site))))]
    }
    if(any(args$col=="random")) {
      col <- sample(allcolors[round(seq(1,50,length.out=length(unique(site))))],length(unique(site)),replace = FALSE)
    }
    if(!any(args$col%in%c("random", "rainbow"))) {
      if(!any(areColors(args$col))) stop("invalid colors", call. = FALSE)
      if(length(args$col)!=length(unique(site))) {
        col = rep(args$col, length(unique(site)))[1:length(unique(site))]
        warning("Length of color vector is not equal to number of sites!", call. = FALSE)
      }
    }
  } else {
    allcolors <- rainbow(60,v=0.85)[1:50]
    col <- sample(allcolors[round(seq(1,50,length.out=length(unique(site))))],length(unique(site)),replace = FALSE)
  }
  
  
  if(sum(names(args)%in%"add")==1) add <- args$add else add = FALSE
  
  
  if(!add) {
    opar <- par(mar = c(6,5,1,1))
    plot(NA, xlim=c(range[1],range[2]), ylim=c(range[3],range[4]), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    map(xlim=c(range[1],range[2]), ylim=c(range[3],range[4]), fill=T, lwd=0.01, col=c("grey90"), add=TRUE)
    mtext(ifelse(sum(names(args)%in%"xlab")==1, args$xlab, "Longitude"), side=1, line=2.2, font=3)
    mtext(ifelse(sum(names(args)%in%"ylab")==1, args$ylab, "Latitude"), side=2, line=2.5, font=3)
    map.axes()
    mtext(ifelse(sum(names(args)%in%"main")==1, args$main, ""), line=0.6, cex=1.2)
  }
  
  
  if(type=="points") {
    points(crds[site>0, ], 
           cex = ifelse(any(names(args)=="cex"), args$cex, 0.5),
           pch = ifelse(any(names(args)=="pch"), args$pch, 16),
           col = col[as.numeric(site)])
  }
  if(type=="cross") {
    for (i in 1:max(unique(site))){
      if(!all(is.na(crds[site==i, 2]))){
        tmp.lon <- quantile(crds[site == i, 1], probs = c(quantiles, 0.5), na.rm = T)
        tmp.lat <- quantile(crds[site == i, 2], probs = c(quantiles, 0.5), na.rm = T)
        points(tmp.lon[3], tmp.lat[3],
               col = col[i],
               cex = ifelse(any(names(args)=="cex"), args$cex, 1),
               pch = ifelse(any(names(args)=="pch"), args$pch, 16))
        segments(tmp.lon[1], tmp.lat[3], tmp.lon[2], tmp.lat[3],
               col = col[i],
               lwd = ifelse(any(names(args)=="lwd"), args$lwd, 2))
        segments(tmp.lon[3], tmp.lat[1], tmp.lon[3], tmp.lat[2],
               col = col[i],
               lwd = ifelse(any(names(args)=="lwd"), args$lwd, 2))
      }
    }
  }    
  
  
  if(hull) {
    for(j in unique(site)){
          if(!all(is.na(crds[site==j, 2]))){
          if(j>0){
            X <- na.omit(crds[site==j,])
            hpts <- chull(X)
            hpts <- c(hpts,hpts[1])
            lines(X[hpts,], 
                 lty = ifelse(sum(names(args)%in%"lty")==1, args$lty, 1),
                       lwd = ifelse(sum(names(args)%in%"lwd")==1, args$lwd, 1),
                       col = col[j])
        }
      }
    }
  }
  
  legend("bottomright", letters[1:max(site)], 
         pch = ifelse(sum(names(args)%in%"pch")==1, args$pch, 16),
         col=col[1:max(as.numeric(site))])
  
  if(!add) par(opar)
}



##' Write a file which plots a trip in Google Earth
##'
##' This function creates a .kml file from light intensity measurements over
##' time that can ve viewed as a trip in Google Earth.
##'
##'
##' @param file A character expression giving the whole path and the name of the
##' resulting output file including the .kml extension.
##' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
##' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
##' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
##' respectively
##' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst}.
##' @param col.scheme the color scheme used for the points. Possible color
##' schemes are: \code{\link{rainbow}}, \code{\link{heat.colors}},
##' \code{\link{topo.colors}}, \code{\link{terrain.colors}}.
##' @param point.alpha a \code{numerical value} indicating the transparency of
##' the point colors on a scale from 0 (transparent) to 1 (opaque).
##' @param cex \code{numerical value} for the size of the points.
##' @param line.col An character expression (any of \code{\link{colors}} or
##' hexadecimal notation), or numeric indicating the color of the line
##' connecting the point locations.
##' @return This function returns no data. It creates a .kml file in the in the
##' defined path.
##' @author Simeon Lisovski and Michael U. Kemp
##' @examples
##' \donttest{
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' filter <- distanceFilter(hoopoe2,distance=30)
##' ## takes time
##' ## trip2kml("trip.kml", hoopoe2$tFirst[filter], hoopoe2$tSecond[filter], hoopoe2$type[filter],
##' ##		degElevation=-6, col.scheme="heat.colors", cex=0.7,
##' ##		line.col="goldenrod")
##'}
##' @importFrom grDevices rgb
##' @export trip2kml
trip2kml <- function(file, tFirst, tSecond, type, degElevation, col.scheme="heat.colors", point.alpha=0.7, cex=1, line.col="goldenrod")
{
  if((length(tFirst)+length(type))!=(length(tSecond)+length(type))) stop("tFirst, tSecond and type must have the same length.")
  
  coord   <- coord(tFirst,tSecond,type,degElevation,note=F)
  index <- !is.na(coord[,2])
  datetime <- as.POSIXct(strptime(paste(ifelse(type==1,substring(tFirst,1,10),substring(tSecond,1,10)),
                                        " ",ifelse(type==1,"12:00:00","00:00:00"),sep=""),format="%Y-%m-%d %H:%M:%S"),"UTC")
  
  coord   <- coord[index,]
  longitude <- coord[,1]
  latitude <- coord[,2]
  
  date <- unlist(strsplit(as.character(datetime[index]), split = " "))[seq(1,
                                                                           ((length(datetime[index]) * 2) - 1), by = 2)]
  time <- unlist(strsplit(as.character(datetime[index]), split = " "))[seq(2,
                                                                           ((length(datetime[index]) * 2)), by = 2)]
  
  if(length(!is.na(coord[,2]))<1) stop("Calculation of coordinates results in zero spatial information.")
  
  if ((col.scheme%in% c("rainbow", "heat.colors", "terrain.colors", "topo.colors",
                        "cm.colors"))==F) stop("The col.scheme has been misspecified.")
  
  seq   <- seq(as.POSIXct(datetime[1]),as.POSIXct(datetime[length(datetime)]),by=12*60*60)
  index2<- ifelse(!is.na(merge(data.frame(d=datetime[index],t=TRUE),data.frame(d=seq,t=FALSE),by="d",all.y=T)[,2]),TRUE,FALSE)
  
  usable.colors <- strsplit(eval(parse(text = paste(col.scheme,
                                                    "(ifelse(length(index2) < 1000, length(index2), 1000), alpha=point.alpha)",
                                                    sep = ""))), split = "")[index2]
  
  usable.line.color <- strsplit(rgb(col2rgb(line.col)[1,
                                                      1], col2rgb(line.col)[2, 1], col2rgb(line.col)[3,
                                                                                                     1], col2rgb(line.col, alpha = 1)[4, 1], maxColorValue = 255),
                                split = "")
  
  date <- unlist(strsplit(as.character(datetime), split = " "))[seq(1,
                                                                    ((length(datetime) * 2) - 1), by = 2)]
  time <- unlist(strsplit(as.character(datetime), split = " "))[seq(2,
                                                                    ((length(datetime) * 2)), by = 2)]
  scaling.parameter <- rep(cex, length(latitude))
  
  data.variables <- NULL
  filename <- file
  write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", filename)
  write("<kml xmlns=\"http://www.opengis.net/kml/2.2\">", filename,
        append = TRUE)
  write("<Document>", filename, append = TRUE)
  write(paste("<name>", filename, "</name>", sep = " "),
        filename, append = TRUE)
  
  write("  <open>1</open>", filename, append = TRUE)
  write("\t<description>", filename, append = TRUE)
  write("\t  <![CDATA[Generated using <a href=\"http://simeonlisovski.wordpress.com/geolight\">GeoLight</a>]]>",
        filename, append = TRUE)
  write("\t</description>", filename, append = TRUE)
  write("<Folder>", filename, append = TRUE)
  write("  <name>Points</name>", filename, append = TRUE)
  write("<open>0</open>", filename, append = TRUE)
  for (i in 1:length(latitude)) {
    write("<Placemark id='point'>", filename, append = TRUE)
    write(paste("<name>", as.character(as.Date(datetime[i])), "</name>", sep = ""),
          filename, append = TRUE)
    write("  <TimeSpan>", filename, append = TRUE)
    write(paste("    <begin>", date[i], "T", time[i], "Z</begin>",
                sep = ""), filename, append = TRUE)
    write(paste("    <end>", date[ifelse(i == length(latitude),
                                         i, i + 1)], "T", time[ifelse(i == length(latitude),
                                                                      i, i + 1)], "Z</end>", sep = ""), filename, append = TRUE)
    write("  </TimeSpan>", filename, append = TRUE)
    write("<visibility>1</visibility>", filename, append = TRUE)
    write("<description>", filename, append = TRUE)
    write(paste("<![CDATA[<TABLE border='1'><TR><TD><B>Variable</B></TD><TD><B>Value</B></TD></TR><TR><TD>Date/Time</TD><TD>",
                datetime[i], "</TD></TR><TR><TD>lat long</TD><TD>",
                paste(latitude[i], longitude[i], sep = " "),
                "</TABLE>]]>", sep = "", collapse = ""), filename,
          append = TRUE)
    write("</description>", filename, append = TRUE)
    write("\t<Style>", filename, append = TRUE)
    write("\t<IconStyle>", filename, append = TRUE)
    write(paste("\t\t<color>", paste(noquote(usable.colors[[i]][c(8,
                                                                  9, 6, 7, 4, 5, 2, 3)]), collapse = ""), "</color>",
                sep = ""), filename, append = TRUE)
    write(paste("  <scale>", scaling.parameter[i], "</scale>",
                sep = ""), filename, append = TRUE)
    write("\t<Icon>", filename, append = TRUE)
    write("\t\t<href>http://maps.google.com/mapfiles/kml/pal2/icon26.png</href>",
          filename, append = TRUE)
    write("\t</Icon>", filename, append = TRUE)
    write("\t</IconStyle>", filename, append = TRUE)
    write("\t</Style>", filename, append = TRUE)
    write("\t<Point>", filename, append = TRUE)
    write(paste("\t<altitudeMode>", "relativeToGround", "</altitudeMode>",
                sep = ""), filename, append = TRUE)
    write("<tesselate>1</tesselate>", filename, append = TRUE)
    write("<extrude>1</extrude>", filename, append = TRUE)
    write(paste("\t  <coordinates>", longitude[i], ",", latitude[i],
                ",", 1,
                "</coordinates>", sep = ""), filename, append = TRUE)
    write("\t</Point>", filename, append = TRUE)
    write(" </Placemark>", filename, append = TRUE)
  }
  
  write("</Folder>", filename, append = TRUE)
  write("<Placemark>", filename, append = TRUE)
  write("  <name>Line Path</name>", filename, append = TRUE)
  write("  <Style>", filename, append = TRUE)
  write("    <LineStyle>", filename, append = TRUE)
  write(paste("\t<color>", paste(noquote(usable.line.color[[1]][c(8,
                                                                  9, 6, 7, 4, 5, 2, 3)]), collapse = ""), "</color>", sep = ""),
        filename, append = TRUE)
  
  
  write(paste("      <width>1</width>", sep = ""), filename,
        append = TRUE)
  write("    </LineStyle>", filename, append = TRUE)
  write("  </Style>", filename, append = TRUE)
  write("  <LineString>", filename, append = TRUE)
  write("    <extrude>0</extrude>", filename, append = TRUE)
  write("    <tessellate>1</tessellate>", filename, append = TRUE)
  write(paste("\t<altitudeMode>clampToGround</altitudeMode>",
              sep = ""), filename, append = TRUE)
  write(paste("     <coordinates>", noquote(paste(longitude,
                                                  ",", latitude, sep = "", collapse = " ")), "</coordinates>",
              sep = ""), filename, append = TRUE)
  write("    </LineString>", filename, append = TRUE)
  write("</Placemark>", filename, append = TRUE)
  write("</Document>", filename, append = TRUE)
  write("</kml>", filename, append = TRUE)
}

#' Draw the positions and the trip on a map
#'
#' Draw a map (from the \code{R} Package \code{maps}) with calculated positions
#' connected by a line
#'
#'
#' @param crds a \code{SpatialPoints} or \code{matrix} object, containing x
#' and y coordinates (in that order).
#' @param equinox logical; if \code{TRUE}, the equinox period(s) is shown as a
#' broken blue line.
#' @param map.range some possibilities to choose defined areas (default:
#' "World").
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @param legend \code{logical}; if \code{TRUE}, a legend will be added to the plot.
#' @author Simeon Lisovski
#' @examples
#'
#' data(hoopoe2)
#'  hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
#'  hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
#' crds <- coord(hoopoe2, degElevation = -6)
#' tripMap(crds, xlim = c(-20,20), ylim = c(0,60), main="hoopoe2")
#'
#' @importFrom maps map map.axes
#' @importFrom graphics lines legend mtext par plot points
#' @export tripMap
tripMap <- function(crds, equinox=TRUE, map.range=c("EuroAfrica","AustralAsia","America","World"), legend = TRUE, ...) {
  
  args <- list(...)
  if(all(map.range==c("EuroAfrica", "AustralAsia", "America", "World")) & sum(names(args)%in%c("xlim", "ylim"))!=2) {
    range <- c(-180, 180, -80, 90)
  } 
  if(all(map.range=="EuroAfrica")) range <- c(-24, 55, -55, 70)
  if(all(map.range=="AustralAsia")) range <- c(60, 190, -55, 78)
  if(all(map.range=="America")) range <- c(-170, -20, -65, 78)
  if(all(map.range=="World")) range <- c(-180, 180, -75, 90)
  
  if(sum(names(args)%in%c("xlim", "ylim"))==2) range <- c(args$xlim, args$ylim)
  
  if(sum(names(args)%in%"add")==1) add <- args$add else add = FALSE
  
  if(!add) {
    opar <- par(mar = c(6,5,1,1))
    plot(NA,xlim=c(range[1],range[2]),ylim=c(range[3],range[4]), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]), fill=T,lwd=0.01,col=c("grey90"),add=TRUE)
    map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),interior=TRUE,col=c("darkgrey"),add=TRUE)
    mtext(ifelse(sum(names(args)%in%"xlab")==1, args$xlab, "Longitude"), side=1, line=2.2, font=3)
    mtext(ifelse(sum(names(args)%in%"ylab")==1, args$ylab, "Latitude"), side=2, line=2.5, font=3)
    map.axes()
    
    mtext(ifelse(sum(names(args)%in%"main")==1, args$main, ""), line=0.6, cex=1.2)
  }
  
  points(crds, 
         pch = ifelse(sum(names(args)%in%"pch")==1, args$pch, 3),
         cex = ifelse(sum(names(args)%in%"cex")==1, args$cex, 0.7))
  lines(crds,
        lwd = ifelse(sum(names(args)%in%"lwd")==1, args$lwd, 0.5),
        col = ifelse(sum(names(args)%in%"col")==1, args$col, "grey10"))
  
  if(equinox){
    nrow <- 1
    repeat{
      while(is.na(crds[nrow,2])==FALSE) {
        nrow <- nrow + 1
        if(nrow==nrow(crds)) break
      }
      if(nrow==nrow(crds)) break
      start   <- nrow-1
      while(is.na(crds[nrow,2])) {
        nrow <- nrow + 1
        if(nrow==nrow(crds)) break
      }
      if(nrow==nrow(crds)) break
      end    <- nrow
      
      lines(c(crds[start,1], crds[end,1]),c(crds[start,2], crds[end,2]), col="blue", lwd=3, lty=1)
    }
    
    if(legend) legend("bottomright", lty=c(0,1,1), pch=c(3,-1,-1), lwd=c(1,0.5,3), col=c("black",ifelse(sum(names(args)%in%"col")==1, args$col, "grey10"),"blue"),c("Positions","Trip","Equinox"),bty="n",bg="grey90",border="grey90",cex=0.8)
  } else {
    if(legend) 	  legend("bottomright",lty=c(0,1),pch=c(3,-1),lwd=c(1,0.5),col=c("black",ifelse(sum(names(args)%in%"col")==1, args$col, "grey10"),"blue"),c("Positions","Trip"),bty="n",bg="grey90",border="grey90",cex=0.8)
  }
  if(!add) par(opar)
  
}

#' Calculate twilight events (sunrise/sunset) from light intensity measurements
#' over time
#'
#' Defines twilight events (sunrise/sunset) at times when the light intensity
#' measurements (\emph{light}) pass the defined light intensity threshold. An
#' interactive plot can be drawn to assess the calculations and improve e.g.
#' select only the realistic events.
#'
#'
#' @param datetime date and time of light intensity measurements e.g.
#' 2008-12-01 08:30 "UTC" (see:
#' \code{\link{as.POSIXct}},\link[=Sys.timezone]{time zones}).
#' @param light \code{numerical} value of the light intensity (usually
#' arbitrary units).
#' @param preSelection codelogical, if TRUE a pre selection of all calculated
#' twiligth events will be offered within the interactive process (ask=TRUE).
#' @param LightThreshold the light intensity threshold for the twilight event
#' calibration. If \code{Default}, it will be set slightly above (3 units) the
#' baseline level (measurement during the night).
#' @param maxLight if the geolocator record the maximum light value of a
#' certain time span, give the interval of maximum recordings in minutes (e.g.
#' 5).
#' @param ask \code{logical}, if TRUE the interactive plot will start after the
#' calculation.
#' @param nsee number of points to plot per screen.
#' @param allTwilights \code{logical}, if TRUE the function returns a list with
#' two tables
#' @return if allTwilights=FALSE, a \code{data frame}. Each row contains two
#' subsequent twilight events (\emph{tFirst, tSecond}) and \emph{type} defining
#' wether \emph{tFirst} refers to sunrise (1) or sunset (2). If
#' allTwilights=TRUE, a \code{list} with the data frame described in the
#' previous sentence and a data frame with all light intensities and a column
#' describing whether each row refers to sunrise (1), sunset (2) or to none of
#' these categories (0).
#' @note Depending on shading during light intensity measurements (e.g. due to
#' vegetation, weather, etc., see Lisovski et \emph{al.} 2012) the light
#' intensities may pass the light intensity threshold several times during the
#' day, resulting false sunrises and sunsets. It is highly recommended to check
#' the derived events visually (\code{ask=TRUE}).Twilight events can be deleted
#' and undeleted by clicking the (first) mouse button at the particular
#' position in the graph. The second mouse buttom (or esc) moves the time
#' series forward. Note, that a backward option is not included.
#' @author Simeon Lisovski
#' @export twilightCalc
#' @importFrom grDevices graphics.off
#' @importFrom graphics abline axis identify legend plot points 
twilightCalc <- function(datetime, light, LightThreshold=TRUE, preSelection=TRUE, maxLight=NULL, ask=TRUE, nsee=500, allTwilights=FALSE)
{
  if(class(datetime)[1]!="POSIXct") {
    stop(sprintf("datetime need to be provided as POSIXct class object."), call. = F)
    
  } else {
    bas <- data.frame(datetime=as.POSIXct(as.character(datetime),"UTC"),light)
  }
  
  
  if (is.numeric(LightThreshold))
  {
    LightThreshold <- as.numeric(LightThreshold)
    min <- min(bas$light)
  } else {
    # Basic level
    r <- as.data.frame(table(bas$light))
    nr <- as.numeric(which.max(r$Freq[as.numeric(r[,1])<mean(bas$light)]))
    LightThreshold <- (as.numeric(as.character(r[nr,1])))+3
  }
  
  out <- i.preSelection(bas$datetime,bas$light, LightThreshold)[,-1]
  
  if(!preSelection) out$mod <- 0
  
  if(ask)
  {
    n   <- nrow(bas)
    nn  <- n%/%nsee + 1
    cutsub <- cut(1:n, nn)
    picks <- NULL
    
    for(i in 1:nn){
      sub <- cutsub == levels(cutsub)[i]
      
      repeat{
        plot(bas[sub,1],bas[sub,2],type="o",cex=0.6,pch=20,ylab="Light intensity",xaxs="i",xaxt="n",xlab="",
             main=paste(as.Date(min(bas[sub,1]))," to ", as.Date(max(bas[sub,1]))," (end: ",as.Date(max(bas[,1])),")",sep=""))
        abline(h=LightThreshold,col="blue",lty=2)
        abline(v=out[out$mod==0,1],col="orange",lty=ifelse(out[out$mod==0,2]==1,1,2))
        points(out[,1],rep(LightThreshold,nrow(out[,])),col=ifelse(out$mod==0,"orange","grey"),pch=20,cex=0.8)
        
        axis(1,at=out[seq(from=1,to=nrow(out),length.out=(nrow(out)%/%2)),1],
             labels=substring(as.character(out[seq(from=1,to=nrow(out),length.out=(nrow(out)%/%2)),1]),6,16),cex=0.7)
        
        legend("topright",lty=c(3,1,2),lwd=c(1.3,2,2),col=c("blue",rep("orange",2)),c("Light\nThreshold","sunrise","sunset"),cex=1,bg="white")
        
        nr <- identify(out[,1],rep(LightThreshold,nrow(out)),n=1,plot=F)
        ifelse(length(nr)>0,ifelse(out$mod[nr]==0,out$mod[nr]<-1,out$mod[nr]<-0),break)
      }
    }
    
    cat("Thank you!\n\n")
    graphics.off()
  }
  
  results <- list()
  
  
  out <- subset(out,out$mod==0)[,-3]
  
  raw <- data.frame(datetime=c(as.POSIXct(datetime,"UTC"),as.POSIXct(out$datetime,"UTC")),
                    light=c(light,rep(LightThreshold,nrow(out))),type=c(rep(0,length(datetime)),out$type))
  
  raw <- raw[order(raw$type),]
  raw <- raw[-which(duplicated(as.character(raw$datetime),fromLast=T)),]
  raw <- raw[order(raw$datetime),]
  
  results$allTwilights <- raw
  
  opt <- data.frame(tFirst=as.POSIXct("1900-01-01 01:01","UTC"),tSecond=as.POSIXct("1900-01-01 01:01","UTC"),type=0)
  row <- 1
  for (k in 1:(nrow(out)-1))
  {
    if (as.numeric(difftime(out[k,1],out[k+1,1]))< 24 & out[k,1] != out[k+1,1])
    {
      opt[row,1] <- out[k,1]
      opt[row,2] <- out[k+1,1]
      opt[row,3] <- out$type[k]
      
      row <- row+1
    }
  }
  
  if(is.numeric(maxLight))
  {
    opt$tFirst[opt$type==2] <- opt$tFirst[opt$type==2] - (maxLight*60)
    opt$tSecond[opt$type==1] <- opt$tSecond[opt$type==1] - (maxLight*60)
  }
  
  if(allTwilights) {
    results$consecTwilights <- opt
    return(results)
  } else {
    return (opt)}
}

#' Example data for calibration: Light intensities and twilight events
#'
#' Light intensity measurements over time (calib1) recorded at the rooftop of
#' the Swiss Ornithological Institute (Lon: 8.0, Lat: 47.01). Defined twilight
#' events from calib1 (calib2). These data serve as an example for calculating
#' the sun elevation angle of an additional data set, which is subsequently
#' used to calibrate the focal dataset.
#'
#' @name calib1
#' @docType data
#' @aliases calib1 calib2 calib
#' @references Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt,
#' F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
#' precision affected by environmental factors. \emph{Methods in Ecology and
#' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
#' @examples
#'
#' data(calib2)
#'  calib2$tFirst <- as.POSIXct(calib2$tFirst, tz = "GMT")
#'  calib2$tSecond <- as.POSIXct(calib2$tSecond, tz = "GMT")
#' getElevation(calib2, known.coord = c(8,47.01))
#'
NULL


#' Light intensity measurements over time recorded on a migratory bird
#'
#' Sunlight intensity measurements over time recorded during the first part of
#' the annual migration of a European Hoopoe (\cite{Upupa epops}). All
#' dates/times are measured in Universal Time Zone (UTC).
#'
#'
#' @name hoopoe1
#' @docType data
#' @format A table with 24474 rows and 2 columns, rows corresponding to light
#' measurements recorded in ten-minute intervals (datetime).
#' @source Baechler, E., Hahn, S., Schaub, M., Arlettaz, R., Jenni, L., Fox,
#' J.W., Afanasyev, V. & Liechti, F. (2010) Year-Round Tracking of Small
#' Trans-Saharan Migrants Using Light-Level Geolocators. \emph{Plos One},
#' \bold{5}.
NULL


#' Sunrise and sunset times: From light intensity measurement (hoopoe1)
#'
#' Sunrise and sunset times derived from light intensity measurements over time
#' (\code{\link{hoopoe1}}). The light measurements corresponding to the first
#' part of the annual migration of a European Hoopoe (\emph{Upupa epops}).
#'
#' @name hoopoe2
#' @docType data
#' @format A table with 340 rows and 3 columns. Each row corresponds to
#' subsequent twilight events ("tFirst" and "tSecond"). The third column
#' ("type") indicates weather the first event is sunrise (1) or sunset (2). All
#' dates/times are measured in Universal Time Zone (UTC).
#' @source Baechler, E., Hahn, S., Schaub, M., Arlettaz, R., Jenni, L., Fox,
#' J.W., Afanasyev, V. & Liechti, F. (2010) Year-Round Tracking of Small
#' Trans-Saharan Migrants Using Light-Level Geolocators. \emph{Plos One},
#' \bold{5}.
#' @examples
#'
#' data(hoopoe2)
#'  hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
#'  hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
#' coord <- coord(hoopoe2, degElevation=-6)
#' ## plot in a map using package maps
#' # par(oma=c(5,0,0,0))
#' # map(xlim=c(-20,40),ylim=c(-10,60),interior=F,col="darkgrey")
#' # map(xlim=c(-20,40),ylim=c(-10,60),boundary=F,lty=2,col="darkgrey",add=T)
#' # mtext(c("Longitude (degrees)","Latitude (degrees)"),side=c(1,2),line=c(2.2,2.5),font=3)
#' # map.axes()
#' # points(coord,col="brown",cex=0.5,pch=20)
#'
NULL