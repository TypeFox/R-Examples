#' APSIM: A general utility package for the Agricultural Production Systems 
#' Simulator
#' 
#' The APSIM package provides a number of general purpose utilities designed to 
#' simplify the importation of simulation output files and to assist in the 
#' building of weather (met) files. Future development will include the ability 
#' to manipulate .apsim files in order to create large factorial simulations.
#' 
#' @section Import output files: \code{\link{loadApsim}} A general purpose 
#'   function to load multiple APSIM results into a single data frame or data 
#'   table. This function can also be used to read a single output file.
#'   
#' @section Manipulate and create met files: \code{\link{prepareMet}} can be
#'   used to generate an APSIM formatted met file from CSV, Excel or netCDF
#'   formats. \code{\link{loadMet}} will import an APSIM formatted met file into
#'   a custom metFile object. \code{\link{writeMetFile}} will write a completed
#'   metFile object to a .met file ready for use in APSIM.
#'   
#' @docType package
#' @name apsim
#' @importFrom data.table rbindlist
#' @import lubridate
#' @import methods
#' @importFrom plyr ddply
#' @import stringr
#' @import sirad
#' @import utils
NULL

#' Weather data from Kingsthorpe, Queensland, Australia.
#' 
#' A dataset containing a year of weather data from Kingsthorpe.
#' 
#' @format A data frame containing 365 rows and 10 columns.
"kingsData"

#' A complete metFile example using data from Dalby, Queensland, Australia.
#' 
#' Three years of data from from Dalby including location coordinates, tav and
#' amp data and units for the weather data
#' 
#' @format See \code{\link{metFile}} for more information on the metFile class.
"met"