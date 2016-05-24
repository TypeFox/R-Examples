#' @name weather_EuropeEU
#' @title Weather serie for Europe EU from NASA POWER agroclimatology
#' @description This contemporary daily climate dataset for Europe
#' covers the period 1st January 2001 to 31 December 2010
#' with 10 complete years of data.
#' It cover a part of Europe) with an elevation less than 500m.
#' The dataset was was extrated from the NASA Langley Research Center POWER Project which provide
#' agroclimatology dataset (Chandler et al., 2004). 
#' It was funded through the NASA Earth Science Directorate Applied Science Program
#' This climate datasetcontains daily estimates of precipitation, mean, minimum and
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
#' @usage weather_EuropeEU
#' @format a \code{RangedData} instance, 1 row per day.
#' SRAD     daily Insolation Incident On A Horizontal Surface (MJ/m^2/day) 
#' T2M      Average Air Temperature At 2 m Above The Surface Of The Earth (degrees C) 
#' TMIN     Minimum Air Temperature At 2 m Above The Surface Of The Earth (degrees C) 
#' TMAX     Maximum Air Temperature At 2 m Above The Surface Of The Earth (degrees C) 
#' RH2M     Relative Humidity At 2 m (%) 
#' TDEW     Dew/Frost Point Temperature At 2 m (degrees C) 
#' RAIN     Average Precipitation (mm/day) 
#' WIND     Wind Speed At 10 m Above The Surface Of The Earth (m/s)
#' @source \url{http://power.larc.nasa.gov/} and \url{http://asterweb.jpl.nasa.gov/gdem.asp} and \url{http://www.geonames.org/about.html}
NULL
