#' Documentation
#' 
#' @description translateSPSS2R provides a set of functions with translated SPSS commands. The usage is oriented on the handling of the SPSS-Syntax. 
#' 
#' @details Mainly the package has two purposes: 
#' \enumerate{
#' \item It facilitates SPSS-Users to change over to R. 
#' \item It facilitates migration projects from SPSS to R.
#' }. For automatic translation a webtool (\url{http://www.eoda.de/en/translate2R.html}) is provided by eoda GmbH. The tool translates SPSS-Syntax to R-Code by the use of translateSPSS2R functions.    
#' 
#' @docType package
#' @name translateSPSS2R-package
#' @param Package: translateSPSS2R
#' @param Type: Package
#' @param Version: 1.0.0
#' @param Date: 2015-23-06
#' @param Imports: car ,data.table, e1071, foreign, Hmisc, plyr, stringr, tidyr, zoo
#' @param  License: GPL-3
#' @author Andreas Wygrabek <Andreas.Wygrabek@@eoda.de> and Bastian Wiessner <Bastian.Wiessner@@eoda.de>
#' @param Maintainer: Andreas Wygrabek
#' @references
#' \url{https://service.eoda.de/translater/?lang=en}
#' \cr \url{http://www.eoda.de/en}
#'
NULL 
#' Sample dataset
#' 
#' Sampledata imported by xpssFrame(). The dataset contains 20 different cars from 3 continents.
#'
#' The variables are as follows:
#'
#' \itemize{
#' \item{V1} {Manufacturer. name of the manufacturer (Audi, BMW, Chevrolet,...)}
#' \item{V2} {Model. name of the model (A8, 328i, Malibu,...)}
#' \item{V3} {Country. numeric indicator for the country (1 = Germany, 2 = US, 3= Japan)}
#' \item{V4} {Car-Type. numeric indicator for the car-type (1 = PKW,2 = SUV)}
#' \item{V5} {Sales volume in thousand. (0.954-230.902)}
#' \item{V5_kl2} {Sales volume grouped. (1 = Until 30.000, 2 = More than 30.000)}
#' \item{V6} {Purchase price in thousand. (15.35-85.50)}
#' \item{V6_kl3} {Purchase price in three groups. (1 = Until 20.000, 2 = more than 20.000 until 30.000, 3 = more than 30.000)}
#' \item{V7_1} {PS. Amount of PS (135-310)}
#' \item{V7_2} {Weight. Weight of the car (94.5-138.60)}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name fromXPSS
#' @usage data(fromXPSS)
#' @format A data.frame with 20 rows and 10 variables.
NULL 
