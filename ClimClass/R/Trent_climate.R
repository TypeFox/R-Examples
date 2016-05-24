 #' Data set of Trentino climate
#'
#' A group of datasets used consistently throughout the ClimClass manual and examples. 
#' It is used as reference definition of the climate for the Trentino region, Italy. 
#' It includes monthly series of temperature and precipitation, and reference tables for 
#' the definition of aridity and continentality/oceanicity
#' 
#' Series like "Txxxx" were supplied by the Autonomous Province of Trento - Meteotrentino (I). 
#' Series like "FEMxx" were supplied by Fondazione Edmund Mach, San Michele all'Adige (I). 
#' 
#' @docType data
#' @keywords datasets
#' @name Trent_climate
#' @usage data(Trent_climate)
NULL
#' Dataset of meteorological measures
#' 
#' A list of 40 data frames (one for each station of the meteoroglical network), with 
#' monthly time series of precipitation and temperature (minimum and maximum).
#' 
#' @docType data
#' @keywords datasets
#' @name lista_cli
#' @usage data(Trent_climate)
#' @format list of 28 elements, each is a data frame of 5 variables and 636 observations
NULL
#' Precipitation
#' 
#' the daily data frame of precipitation for a number of stations. It is used in function  
#' \code{\link{oiv_ind}}
#' 
#' @docType data
#' @keywords datasets
#' @name P
#' @usage data(Trent_climate)
#' @format data frame: 19358 obervations of 15 variables (stations)
NULL
#' Mean daily temperature
#' 
#' the daily data frame of mean daily temperature for a number of stations, used in function 
#'  \code{\link{oiv_ind}}
#'  
#'  
#' @docType data
#' @keywords datasets
#' @name Tm
#' @usage data(Trent_climate)
#' @format data frame: 19358 obervations of 15 variables (stations)
NULL
#' Minimum daily temperature
#' 
#' the daily data frame of minimum daily temperature for a number of stations, used in function 
#'  \code{\link{oiv_ind}}
#'  
#'  
#' @docType data
#' @keywords datasets
#' @name Tn
#' @usage data(Trent_climate)
#' @format data frame: 19358 obervations of 15 variables (stations)
NULL
#' Maximum daily temperature
#' 
#' the daily data frame of maximum daily temperature for a number of stations, used in function 
#'  \code{\link{oiv_ind}}
#'  
#'  
#' @docType data
#' @keywords datasets
#' @name Tx
#' @usage data(Trent_climate)
#' @format data frame: 19358 obervations of 15 variables (stations)
NULL
#' Water balance
#' 
#' is the first list (\code{W_balance}) in \code{thornt_lst} organized according to stations. 
#' See Examples in function \code{\link{thornthwaite}} for its construction.
#'  
#' @docType data
#' @keywords datasets
#' @name W_balance
#' @usage data(Trent_climate)
#' @format list of 28 elements, each is a data frame of 5 variables and 636 observations
NULL
#' Aridity index
#' 
#' Used for reference in aridity indices assessment (see function \code{\link{arid}} and 
#' references for data sources).
#'  
#' @docType data
#' @keywords datasets
#' @name arid_ind_tables
#' @usage data(Trent_climate)
#' @format  list formed by six data frames. 
NULL
#' climatic normals of precipitation and temperatures
#' 
#' climatic normals of precipitation and temperature (minimum, maximum, and mean) for the 
#' climatic period 1981 - 2010. It has been calculated by function \code{\link{climate}}.
#'  
#' @docType data
#' @keywords datasets
#' @name clima_81_10
#' @usage data(Trent_climate)
#' @format  a list (one table for each station) of 40 monthly climatic normals
NULL
#' radiative energy coefficients
#' 
#' "radiative energy coefficients" for Hargreaves' equation, corresponding to the 
#' daily extra-atmospheric solar radiation energy. It is the output of function 
#' \code{\link{ExAtRa}}.
#'  
#' @docType data
#' @keywords datasets
#' @name coeff_rad
#' @usage data(Trent_climate)
#' @format  an array of 12 numerics
NULL
#' continentality/oceanicity indices
#' 
#' Used for reference in continentality/oceanicity indices assessment (see function 
#' \code{\link{contin}} and references for data sources).
#'  
#' @docType data
#' @keywords datasets
#' @name continental_ind_tables
#' @usage data(Trent_climate)
#' @format  list of 4 data frames
NULL
#' geographical position for each meteorological station
#' coordinates and elevation for each station in the data set. Coordinates are geographical 
#' and elevation is measured in meters above mean sea level.
#' 
#' @docType data
#' @keywords datasets
#' @name coord_elev
#' @usage data(Trent_climate)
#' @format  data frame of 40 observations of 4 variables
NULL
#' Monthly quantiles of the meteorological varaibales
#' 
#' the second list (\code{quantiles}) in \code{thornt_lst} organized according to stations. 
#' See Examples in function \code{\link{thornthwaite}} for its construction.
#' 
#' @docType data
#' @keywords datasets
#' @name quantiles
#' @usage data(Trent_climate)
#' @format  data frame of 40 observations of 4 variables
NULL
#' input for the Thornthwaite function
#' 
#' For every station, the first element (a list, \code{Thornth._W._bal}) reports the monthly 
#' series of water balance quantities for the station, each in one data frame (see function 
#' \code{\link{thornthwaite}} for details). The second list (\code{quantiles}) reports the 
#' monthly quantiles for the same quantities.
#' 
#' @docType data
#' @keywords datasets
#' @name thornt_lst
#' @usage data(Trent_climate)
#' @format S3 object: a "hyperlist" (list of lists of lists), one list of lists for each station
NULL
