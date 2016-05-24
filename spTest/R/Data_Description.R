#' Average temperature from WRFG RCM output from NARCCAP
#'
#' A list containing coordinates and averages of yearly average temperatures from a 24-year period of WRFG-NCEP RCM output on 14606 grid boxes over North America (Weather Research and Forecasting - Grell scheme [WRFG] with boundary conditions provided by the National Center for Environmental Prediciton [NCEP] Regional Climate Model [RCM]).
#' 
#' @format A list containing five items:
#' \describe{
#'	\item{lon}{longitude of the grid boxes}
#' 	\item{lat}{latitude of the grid boxes}
#' 	\item{xc}{x coordinates of the grid boxes (Lambert Conformal projection)}
#'	\item{yc}{y coordinates of the grid boxes (Lambert Conformal projection)}
#'	\item{WRFG.NCEP.tas}{matrix of average temperatures}
#'}
#' @details RCM output is available through the North American Regional Climate Change Assessment Program (NARCCAP). Data used for this example are from the WRFG model with NCEP boundary conditions. For the years 1981-2004, grid box temperatures (given by the variable 'tas') are averaged over the entire year. Yearly average temperatures for each grid box were then averaged over the 24 year period to create the matrix WRFG.NCEP.tas.
#'
#' @source \url{https://www.earthsystemgrid.org/project/narccap.html}
#'
#'@examples
#' data(WRFG)
#' library(fields)
#' image.plot(WRFG$lon-360, WRFG$lat, WRFG$WRFG.NCEP.tas)
#' world(add = TRUE)
"WRFG"
