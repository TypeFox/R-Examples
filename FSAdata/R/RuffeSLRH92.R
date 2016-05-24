#' @title Biological data for Ruffe captured from the St. Louis River in 1992.
#' 
#' @description Biological data for Ruffe (\emph{Gymnocephalus cernuus}) captured in the St. Louis River, Lake Superior in 1992.
#' 
#' @name RuffeSLRH92
#' 
#' @docType data
#' 
#' @format A data frame of 738 observations on the following 11 variables:
#' 
#' \describe{
#'  \item{fish.id}{A unique fish identification number (across all years, most of which are not shown in this file}
#'  \item{month}{Month of capture} 
#'  \item{day}{Day of capture}
#'  \item{year}{Year of capture} 
#'  \item{indiv}{A unique fish identification number within the year}
#'  \item{location}{Grid location of capture}
#'  \item{length}{Total length (mm)}
#'  \item{weight}{Weight (g)} 
#'  \item{sex}{Sex factor}
#'  \item{maturity}{Maturity stage factor}
#'  \item{age}{Age (yrs) from scales}
#' }
#' 
#' @section Topic(s):
#'  \itemize{
#'    \item Length Frequency 
#'    \item Weight-Length 
#'    \item Growth
#'    \item von Bertalanffy
#'    \item Maturity
#'  }
#'  
#' @concept Growth 'von Bertalanffy' 'Length Frequency' 'Weight-Length' Maturity
#' 
#' @source personal collection by the United States Geological Survey, Lake Superior Biological Station, Ashland, WI.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(RuffeSLRH92)
#' str(RuffeSLRH92)
#' head(RuffeSLRH92)
#' op <- par(mfrow=c(2,2),pch=19)
#' hist(RuffeSLRH92$length,main="")
#' hist(RuffeSLRH92$age,main="")
#' plot(weight~length,data=RuffeSLRH92)
#' plot(length~age,data=RuffeSLRH92)
#' par(op)
#' xtabs(~age,data=RuffeSLRH92)
#' xtabs(~sex,data=RuffeSLRH92)
#' xtabs(~maturity,data=RuffeSLRH92)
#' 
NULL