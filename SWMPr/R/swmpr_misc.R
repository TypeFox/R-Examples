######
#' Create a swmpr object
#' 
#' Wrapper for creating a swmpr object
#' 
#' @param  stat_in \code{data.frame} of swmp data
#' @param  meta_in chr string for station code (7 or 8 characters), can be multiple stations if data are combined
#' 
#' @export swmpr
#' 
#' @return Returns a swmpr object to be used with S3 methods
#' 
#' @details 
#' This function is a simple wrapper to \code{\link[base]{structure}} that is used internally within other functions to create a swmpr object.  The function does not have to be used explicitly.  Attributes of a swmpr object include \code{names}, \code{row.names}, \code{class}, \code{station}, \code{parameters}, \code{qaqc_cols}, \code{date_rng}, \code{timezone}, \code{stamp_class}, \code{metabolism} (if present), and \code{metab_units} (if present). 
#' 
swmpr <- function(stat_in, meta_in){
    
  if(!is.data.frame(stat_in)) 
    stop('stat_in must be data.frame')
  
  # qaqc attribute
  qaqc_cols <- FALSE
  if(any(grepl('^f_', names(stat_in)))) qaqc_cols <- TRUE
  
  # parameters attribute
  parameters <- grep('datetimestamp|^f_', names(stat_in), invert = TRUE, value = TRUE)
  
  # get stations, param_types attribtues
  param_types <- param_names()
  param_types <- unlist(lapply(param_types, function(x) any(x %in% parameters)))
  param_types <- names(param_names())[param_types]
  station <- grep(paste0(param_types, collapse = '|'), meta_in, value = TRUE)

  # timezone using time_vec function
  timezone <- time_vec(station_code = station, tz_only = TRUE)

  # create class, with multiple attributes
  structure(
    .Data = stat_in, 
    class = c('swmpr', 'data.frame'), 
    station = station,
    parameters = parameters, 
    qaqc_cols = qaqc_cols,
    date_rng = range(stat_in$datetimestamp),
    timezone = timezone, 
    stamp_class = class(stat_in$datetimestamp),
    metabolism = NULL, 
    metab_units = NULL
    )
  
}

######
#' Parse web results for swmpr
#' 
#' Parsing function for objects returned from CDMO web services
#' 
#' @param  resp_in web object returned from CDMO server, response class from httr package
#' @param  parent_in chr string of parent nodes to parse
#' 
#' @import XML
#' 
#' @export
#' 
#' @details 
#' This function parses XML objects returned from the CDMO server, which are further passed to \code{\link{swmpr}}.  It is used internally by the data retrieval functions, excluding \code{\link{import_local}}.  The function does not need to be called explicitly.
#' 
#' @return Returns a \code{data.frame} of parsed XML nodes
parser <- function(resp_in, parent_in = 'data'){

  # convert to XMLDocumentContent for parsing
  raw <- xmlTreeParse(resp_in, useInternalNodes = TRUE)

  # get parent data nodes
  parents <- xpathSApply(
    raw,
    paste0('//', parent_in)
  )
  
  # get children nodes from data parents
  out <- lapply(parents, 
    function(x) getChildrenStrings(x)
    )
  out <- do.call('rbind', out)
  out <- data.frame(out)
  names(out) <- tolower(names(out))
  
  # return output
  return(out)
  
}

######
#' Format SWMP datetimestamp
#'
#' Format the datetimestamp column of SWMP data
#' 
#' @param  chr_in chr string of datetimestamp vector
#' @param  station_code is chr string for station (three or more characters)
#' @param  tz_only logical that returns only the timezone, default \code{FALSE}
#' 
#' @export
#' 
#' @return  Returns a POSIX vector if \code{tz_only} is true, otherwise the timezone for a station is returned as a chr string
#' 
#' @details 
#' This function formats the datetimestamp column of SWMP data to the \code{\link[base]{POSIXct}} format and the correct timezone for a station.  Note that SWMP data do not include daylight savings and the appropriate location based on GMT offsets is used for formatting.  This function is used internally within data retrieval functions and does not need to be called explicitly.
time_vec <- function(chr_in = NULL, station_code, tz_only = FALSE){
  
  # lookup table for time zones based on gmt offset - no DST!
  gmt_tab <- data.frame(
    gmt_off=c(-4,-5,-6,-8,-9),
    tz = c('America/Virgin', 'America/Jamaica', 'America/Regina',
      'Pacific/Pitcairn', 'Pacific/Gambier'),
    stringsAsFactors = FALSE
    )
  
  # hard-coded gmt offset for each site, from metadata direct from CDMO
  sites <- c('ace', 'apa', 'cbm', 'cbv', 'del', 'elk', 
    'gnd', 'grb', 'gtm', 'hud', 'jac', 'job', 'kac', 
    'lks', 'mar', 'nar', 'niw', 'noc', 'owc', 'pdb', 
    'rkb', 'sap', 'sfb', 'sos', 'tjr', 'wel', 'wkb',
    'wqb')
  gmt_offsets <- c(-5, -5, -5, -5, -5, -8, -6, -5, -5, -5, -5, -4, 
    -9, -6, -6, -5, -5, -5, -5, -8, -5, -5, -8, -8, -8,
    -5, -6, -5)  
  
  # get timezone from above information
  gmt_offset <- gmt_offsets[which(sites %in% substr(station_code, 1, 3))]
  tzone <- gmt_tab[gmt_tab$gmt_off %in% gmt_offset, 'tz']

  # timezone only if T
  if(tz_only) return(tzone)
  
  # format datetimestamp
  out <- as.POSIXct(chr_in, tz = tzone, format = '%m/%d/%Y %H:%M')
  
  # return output
  return(out)
  
}

######
#' Get parameters of a given type
#'
#' Get parameter column names for each parameter type
#' 
#' @param  param_type chr string specifying \code{'nut'}, \code{'wq'}, or \code{'met'}.  Input can be one to three types.
#' 
#' @export
#' 
#' @return Returns a named list of parameters for the \code{param_type}.  The parameter names are lower-case strings of SWMP parameters and corresponding qaqc names (\code{'f_'} prefix)
#' 
#' @details
#' This function is used internally within several functions to return a list of the expected parameters for a given parameter type: nutrients, water quality, or meteorological.  It does not need to be called explicitly. 
param_names <- function(param_type = c('nut', 'wq', 'met')){
  
  # sanity check
  if(any(!param_type %in% c('nut', 'wq', 'met')))
    stop('param_type must chr string of nut, wq, or met')
  
  nut_nms <- c('po4f', 'chla_n', 'no3f', 'no2f', 'nh4f', 'no23f', 'ke_n',
    'urea')
  nut_nms <- paste0(c('', 'f_'), rep(nut_nms, each = 2))
  
  wq_nms <- c('temp', 'spcond', 'sal', 'do_pct', 'do_mgl', 'depth', 
    'cdepth', 'level', 'clevel', 'ph', 'turb', 'chlfluor')
  wq_nms <- paste0(c('', 'f_'), rep(wq_nms, each = 2))
  
  met_nms <- c('atemp', 'rh', 'bp', 'wspd', 'maxwspd', 'wdir', 'sdwdir',
    'totpar', 'totprcp', 'cumprcp', 'totsorad')
  met_nms <- paste0(c('', 'f_'), rep(met_nms, each = 2))
  
  # get names for a given parameter type
  out <- sapply(param_type, function(x) get(paste0(x, '_nms')), simplify = FALSE)
  
  return(out)
  
}

######
#' Locations of NERRS sites
#'
#' Location of NERRS sites in decimal degress and time offset from Greenwich mean time.  Only active sites as of January 2015 are included.  Sites are identified by five letters indexing the reserve and site names.  The dataset is used to plot locations with the \code{\link{map_reserve}} function and to identify metabolic days with the \code{\link{ecometab}} function. 
#' 
#' @format A \code{\link[base]{data.frame}} object with 161 rows and 4 variables:
#' \describe{
#'   \item{\code{station_code}}{chr}
#'   \item{\code{latitude}}{numeric}
#'   \item{\code{longitude}}{numeric}
#'   \item{\code{gmt_off}}{int}
#' }
#' 
#' @seealso \code{\link{ecometab}}, \code{\link{map_reserve}}
"stat_locs"

######
#' Example nutrient data for Apalachicola Bay Cat Point station.
#'
#' An example nutrient dataset for Apalachicola Bay Cat Point station.  The data are a \code{\link{swmpr}} object that have been imported into R from csv files using the \code{\link{import_local}} function.  The raw data were obtained from the CDMO data portal but can also be accessed from a zip file created for this package.  See the source below.  The help file for the \code{\link{import_local}} function describes how the data can be imported from the zip file.  Attributes of the dataset include \code{names}, \code{row.names}, \code{class}, \code{station}, \code{parameters}, \code{qaqc_cols}, \code{date_rng}, \code{timezone}, and \code{stamp_class}. 
#'  
#' @format A \code{\link{swmpr}} object and \code{\link[base]{data.frame}} with 215 rows and 13 variables:
#' \describe{
#'   \item{\code{datetimestamp}}{POSIXct}
#'   \item{\code{po4f}}{num}
#'   \item{\code{f_po4f}}{chr}
#'   \item{\code{nh4f}}{num}
#'   \item{\code{f_nh4f}}{chr}
#'   \item{\code{no2f}}{num}
#'   \item{\code{f_no2f}}{chr}
#'   \item{\code{no3f}}{num}
#'   \item{\code{f_no3f}}{chr}
#'   \item{\code{no23f}}{num}
#'   \item{\code{f_no23f}}{chr}
#'   \item{\code{chla_n}}{num}
#'   \item{\code{f_chla_n}}{chr}
#' }
#' 
#' @source \url{https://s3.amazonaws.com/swmpexdata/zip_ex.zip}
#' 
#' @examples 
#' data(apacpnut)
"apacpnut"

######
#' Example water quality data for Apalachicola Bay Cat Point station.
#'
#' An example water quality dataset for Apalachicola Bay Cat Point station.  The data are a \code{\link{swmpr}} object that have been imported into R from csv files using the \code{\link{import_local}} function.  The raw data were obtained from the CDMO data portal but can also be accessed from a zip file created for this package.  See the source below.  The help file for the \code{\link{import_local}} function describes how the data can be imported from the zip file.  Attributes of the dataset include \code{names}, \code{row.names}, \code{class}, \code{station}, \code{parameters}, \code{qaqc_cols}, \code{date_rng}, \code{timezone}, and \code{stamp_class}. 
#'  
#' @format A \code{\link{swmpr}} object and \code{\link[base]{data.frame}} with 70176 rows and 25 variables:
#' \describe{
#'   \item{\code{datetimestamp}}{POSIXct}
#'   \item{\code{temp}}{num}
#'   \item{\code{f_temp}}{chr}
#'   \item{\code{spcond}}{num}
#'   \item{\code{f_spcond}}{chr}
#'   \item{\code{sal}}{num}
#'   \item{\code{f_sal}}{chr}
#'   \item{\code{do_pct}}{num}
#'   \item{\code{f_do_pct}}{chr}
#'   \item{\code{do_mgl}}{num}
#'   \item{\code{f_do_mgl}}{chr}
#'   \item{\code{depth}}{num}
#'   \item{\code{f_depth}}{chr}
#'   \item{\code{cdepth}}{num}
#'   \item{\code{f_cdepth}}{chr}
#'   \item{\code{level}}{num}
#'   \item{\code{f_level}}{chr}
#'   \item{\code{clevel}}{num}
#'   \item{\code{f_clevel}}{chr}
#'   \item{\code{ph}}{num}
#'   \item{\code{f_ph}}{chr}
#'   \item{\code{turb}}{num}
#'   \item{\code{f_turb}}{chr}
#'   \item{\code{chlfluor}}{num}
#'   \item{\code{f_chlfluor}}{chr}
#' }
#' 
#' @source \url{https://s3.amazonaws.com/swmpexdata/zip_ex.zip}
#'
#' @examples 
#' data(apacpwq)
"apacpwq"

######
#' Example water quality data for Apalachicola Bay Dry Bar station.
#'
#' An example water quality dataset for Apalachicola Bay Dry Bar station.  The data are a \code{\link{swmpr}} object that have been imported into R from csv files using the \code{\link{import_local}} function.  The raw data were obtained from the CDMO data portal but can also be accessed from a zip file created for this package.  See the source below.  The help file for the \code{\link{import_local}} function describes how the data can be imported from the zip file.  Attributes of the dataset include \code{names}, \code{row.names}, \code{class}, \code{station}, \code{parameters}, \code{qaqc_cols}, \code{date_rng}, \code{timezone}, and \code{stamp_class}. 
#'  
#' @format A \code{\link{swmpr}} object and \code{\link[base]{data.frame}} with 70176 rows and 25 variables:
#' \describe{
#'   \item{\code{datetimestamp}}{POSIXct}
#'   \item{\code{temp}}{num}
#'   \item{\code{f_temp}}{chr}
#'   \item{\code{spcond}}{num}
#'   \item{\code{f_spcond}}{chr}
#'   \item{\code{sal}}{num}
#'   \item{\code{f_sal}}{chr}
#'   \item{\code{do_pct}}{num}
#'   \item{\code{f_do_pct}}{chr}
#'   \item{\code{do_mgl}}{num}
#'   \item{\code{f_do_mgl}}{chr}
#'   \item{\code{depth}}{num}
#'   \item{\code{f_depth}}{chr}
#'   \item{\code{cdepth}}{num}
#'   \item{\code{f_cdepth}}{chr}
#'   \item{\code{level}}{num}
#'   \item{\code{f_level}}{chr}
#'   \item{\code{clevel}}{num}
#'   \item{\code{f_clevel}}{chr}
#'   \item{\code{ph}}{num}
#'   \item{\code{f_ph}}{chr}
#'   \item{\code{turb}}{num}
#'   \item{\code{f_turb}}{chr}
#'   \item{\code{chlfluor}}{num}
#'   \item{\code{f_chlfluor}}{chr}
#' }
#' 
#' @source \url{https://s3.amazonaws.com/swmpexdata/zip_ex.zip}
#' 
#' @examples 
#' data(apadbwq)
"apadbwq"

######
#' Example weather data for Apalachicola Bay East Bay station.
#'
#' An example weather dataset for Apalachicola Bay East Bay station.  The data are a \code{\link{swmpr}} object that have been imported into R from csv files using the \code{\link{import_local}} function.  The raw data were obtained from the CDMO data portal but can also be accessed from a zip file created for this package.  See the source below.  The help file for the \code{\link{import_local}} function describes how the data can be imported from the zip file.  Attributes of the dataset include \code{names}, \code{row.names}, \code{class}, \code{station}, \code{parameters}, \code{qaqc_cols}, \code{date_rng}, \code{timezone}, and \code{stamp_class}. 
#'  
#' @format A \code{\link{swmpr}} object and \code{\link[base]{data.frame}} with 70176 rows and 23 variables:
#' \describe{
#'   \item{\code{datetimestamp}}{POSIXct}
#'   \item{\code{atemp}}{num}
#'   \item{\code{f_atemp}}{chr}
#'   \item{\code{rh}}{num}
#'   \item{\code{f_rh}}{chr}
#'   \item{\code{bp}}{num}
#'   \item{\code{f_bp}}{chr}
#'   \item{\code{wspd}}{num}
#'   \item{\code{f_wspd}}{chr}
#'   \item{\code{maxwspd}}{num}
#'   \item{\code{f_maxwspd}}{chr}
#'   \item{\code{wdir}}{num}
#'   \item{\code{f_wdir}}{chr}
#'   \item{\code{sdwdir}}{num}
#'   \item{\code{f_sdwdir}}{chr}
#'   \item{\code{totpar}}{num}
#'   \item{\code{f_totpar}}{chr}
#'   \item{\code{totprcp}}{num}
#'   \item{\code{f_totprcp}}{chr}
#'   \item{\code{cumprcp}}{num}
#'   \item{\code{f_cumprcp}}{chr}
#'   \item{\code{totsorad}}{num}
#'   \item{\code{f_totsorad}}{chr}
#' }
#' 
#' @source \url{https://s3.amazonaws.com/swmpexdata/zip_ex.zip}
#'
#' @examples 
#' data(apaebmet)
"apaebmet"

######
#' Identify metabolic days in a swmpr time series
#'
#' Identify metabolic days in a swmpr time series based on sunrise and sunset times for a location and date.  The metabolic day is considered the 24 hour period between sunsets for two adjacent calendar days.  The function calls the \code{\link[maptools]{sunriset}} function from the maptools package, which uses algorithms from the National Oceanic and Atmospheric Administration (\url{http://www.esrl.noaa.gov/gmd/grad/solcalc/}).
#' 
#' @param dat_in data.frame
#' @param stat_in chr vector of station name including data type
#' 
#' @import maptools
#' 
#' @export 
#' 
#' @details This function is only used within \code{\link{ecometab}} and should not be called explicitly.
#' 
#' @seealso 
#' \code{\link{ecometab}}, \code{\link[maptools]{sunriset}}
#' 
metab_day <- function(dat_in, stat_in){
  
  # station locations
  dat_locs <- get('stat_locs')
  stat_meta <- dat_locs[grep(gsub('wq$', '', stat_in), dat_locs$station_code),]
  
  # all times are standard - no DST!
  gmt_tab <- data.frame(
    gmt_off=c(-4, -5, -6, -8, -9),
    tz = c('America/Virgin', 'America/Jamaica', 'America/Regina',
      'Pacific/Pitcairn', 'Pacific/Gambier'),
    stringsAsFactors = F
    )
  
  # get sunrise/sunset times using functions from maptools
  # adapted from streammetabolism sunrise.set function
  lat <- stat_meta$latitude
  long <- stat_meta$longitude
  gmt_off <- stat_meta$gmt_off
  tz <- gmt_tab[gmt_tab$gmt_off == gmt_off, 'tz']
  start_day <- format(
    dat_in$datetimestamp[which.min(dat_in$datetimestamp)] - (60 * 60 * 24), 
    format = '%Y/%m/%d'
    )
  tot_days <- 1 + length(unique(as.Date(dat_in$datetimestamp)))
  lat.long <- matrix(c(long, lat), nrow = 1)
  sequence <- seq(
    from = as.POSIXct(start_day, tz = tz), 
    length.out = tot_days, 
    by = "days"
    )
  sunrise <- sunriset(lat.long, sequence, direction = "sunrise", 
      POSIXct = TRUE)
  sunset <- sunriset(lat.long, sequence, direction = "sunset", 
      POSIXct = TRUE)
  ss_dat <- data.frame(sunrise, sunset)
  ss_dat <- ss_dat[, -c(1, 3)]
  colnames(ss_dat) <- c("sunrise", "sunset")
  
  # remove duplicates, if any
  ss_dat <- ss_dat[!duplicated(strftime(ss_dat[, 1], format = '%Y-%m_%d')), ]
  ss_dat <- data.frame(
    ss_dat,
    metab_date = as.Date(ss_dat$sunrise, tz = tz)
    )
  ss_dat <- reshape2::melt(ss_dat, id.vars = 'metab_date')
  if(!"POSIXct" %in% class(ss_dat$value))
    ss_dat$value <- as.POSIXct(ss_dat$value, origin='1970-01-01',tz=tz)
  ss_dat <- ss_dat[order(ss_dat$value),]
  ss_dat$day_hrs <- unlist(lapply(
    split(ss_dat, ss_dat$metab_date),
    function(x) rep(as.numeric(x[2, 'value'] - x[1, 'value']), 2) 
    ))
  names(ss_dat)[names(ss_dat) %in% c('variable', 'value')] <- c('solar_period', 'solar_time')
  
  # matches is vector of row numbers indicating starting value that each
  # unique datetimestamp is within in ss_dat
  # output is meteorological day matches appended to dat_in
  matches <- findInterval(dat_in$datetimestamp, ss_dat$solar_time)
  out <- data.frame(dat_in, ss_dat[matches, ])
  return(out)
      
  }

######
#' Calculate oxygen mass transfer coefficient
#' 
#' Calculate oxygen mass transfer coefficient using equations in Thiebault et al. 2008.  Output is used to estimate the volumetric reaeration coefficient for ecosystem metabolism.
#'
#' @param temp numeric for water temperature (C)
#' @param sal numeric for salinity (ppt)
#' @param atemp numeric for air temperature (C)
#' @param wspd numeric for wind speed (m/s)
#' @param bp numeric for barometric pressure (mb)
#' @param height numeric for height of anemometer (meters)
#'
#' @import oce
#' 
#' @export
#' 
#' @details
#' This function is used within the \code{\link{ecometab}} function and should not be used explicitly.
#' 
#' @references
#' Ro KS, Hunt PG. 2006. A new unified equation for wind-driven surficial oxygen transfer into stationary water bodies. Transactions of the American Society of Agricultural and Biological Engineers. 49(5):1615-1622.
#' 
#' Thebault J, Schraga TS, Cloern JE, Dunlavey EG. 2008. Primary production and carrying capacity of former salt ponds after reconnection to San Francisco Bay. Wetlands. 28(3):841-851.
#' 
#' @seealso 
#' \code{\link{ecometab}}
#' 
calckl <- function(temp, sal, atemp, wspd, bp, height = 10){
  
  #celsius to kelvin conversion
  CtoK <- function(val) val + 273.15 
    
  Patm <- bp * 100; # convert from millibars to Pascals
  zo <- 1e-5; # assumed surface roughness length (m) for smooth water surface
  U10 <- wspd * log(10 / zo) / log(height / zo)
  tempK <- CtoK(temp)
  atempK <- CtoK(atemp)
  sigT <- swSigmaT(sal, temp, 10) # set for 10 decibars = 1000mbar = 1 bar = 1atm
  rho_w <- 1000 + sigT #density of SW (kg m-3)
  Upw <- 1.002e-3 * 10^((1.1709 * (20 - temp) - (1.827 * 10^-3 * (temp - 20)^2)) / (temp + 89.93)) #dynamic viscosity of pure water (sal + 0);
  Uw <- Upw * (1 + (5.185e-5 * temp + 1.0675e-4) * (rho_w * sal / 1806.55)^0.5 + (3.3e-5 * temp + 2.591e-3) * (rho_w * sal / 1806.55))  # dynamic viscosity of SW
  Vw <- Uw / rho_w  #kinematic viscosity
  Ew <- 6.112 * exp(17.65 * atemp / (243.12 + atemp))  # water vapor pressure (hectoPascals)
  Pv <- Ew * 100 # Water vapor pressure in Pascals
  Rd <- 287.05  # gas constant for dry air ( kg-1 K-1)
  Rv <- 461.495  # gas constant for water vapor ( kg-1 K-1)
  rho_a <- (Patm - Pv) / (Rd * atempK) + Pv / (Rv * tempK)
  kB <- 1.3806503e-23 # Boltzman constant (m2 kg s-2 K-1)
  Ro <- 1.72e-10     #radius of the O2 molecule (m)
  Dw <- kB * tempK / (4 * pi * Uw * Ro)  #diffusivity of O2 in water 
  KL <- 0.24 * 170.6 * (Dw / Vw)^0.5 * (rho_a / rho_w)^0.5 * U10^1.81  #mass xfer coef (m d-1)
  
  return(KL)
  
  }

#' @importFrom stats var
NULL