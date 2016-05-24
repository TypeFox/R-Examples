#' @name weather_SouthAsia
#' @title Weather series for southern Asia from NASA POWER agroclimatology
#' @description This contemporary daily climate dataset for South Asia
#' covers the period 1st January 1997 to 31 December 2008
#' with 12 complete years of data with precipitation.
#' It cover a part of South Asia (North-East of India, Bangladesh, Myanmar, Neapal) with an elevation less than 2500m.
#' The dataset was extracted from the NASA Langley Research Center POWER Project which provide
#' agroclimatology dataset (Chandler et al., 2004). 
#' It was funded through the NASA Earth Science Directorate Applied Science Program
#' This climate dataset contains daily estimates of precipitation, mean, minimum and
#' maximum temperature, relative humidity, dew point, solar radiation
#' and wind speed with global coverage at one degree resolution
#' (approximately 111 km at the equator).
#' The NASA POWER agroclimatology data are derived from various
#' sources: solar radiation from satellite observations, meteorological
#' data from the Goddard Earth Observing System global assimilation
#' model version 4 (GEOS-4), and precipitation from the Global
#' Precipitation Climate Project and Topical Rainfall Measurement
#' Mission. A full description can be found at 
#' \url{http://power.larc.nasa.gov/common/MethodologySSE6/POWER_Methodology_Content.html}
#' Elevation (Altitude) were retrive from Aster Global Digital Elevation Model
#' by using the Webservice api.geonames.org/astergdem?
#' Sample are: ca 30m x 30m, between 83N and 65S latitude.
#' Result : a single number giving the elevation in meters according to aster gdem, ocean areas have been masked as "no data" and have been assigned a value of -9999
#' Example http://api.geonames.org/astergdem?lat=50.01&lng=10.2&username=demo 
#' @docType data
#' @usage weather_SouthAsia
#' @format a \code{RangedData} instance, 1 row per day.
#' @source \url{http://power.larc.nasa.gov/} and \url{http://asterweb.jpl.nasa.gov/gdem.asp} and \url{http://www.geonames.org/about.html}
NULL
