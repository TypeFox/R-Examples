#' @title Stock and recruitment data for Exmouth Gulf Tiger Prawn, 1970-83.
#' 
#' @description Stock and recruitment data for Exmouth Gulf Tiger Prawn (\emph{Panaeus esculentus}), 1970-1983.
#' 
#' @name TPrawnsEG
#' 
#' @docType data
#' 
#' @format A data frame with 14 observations on the following 5 variables.
#'  \describe{ 
#'    \item{year}{a numeric vector of years}
#'    \item{stock}{a numeric vector giving the index of spawning stock fish}
#'    \item{recruits}{a numeric vector containing the index of recruits}
#'    \item{cycloneJan}{a numeric vector containing the relative rainfal in January as an index of cyclonic activity}
#'    \item{cycloneFeb}{a numeric vector containing the relative rainfal in February as an index of cyclonic activity}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Stock-Recruit
#'    \item Recruitment
#'  }
#' 
#' @concept 'Stock-Recruit' Recruitment
#' 
#' @source From table 9.1 in Haddon, M. 2000.  Modeling and Quantitative Methods in Fisheries, CRC Press.  Originally from Penn, J. W., and Caputi, N.  1986.  Spawning stock-recruitment relationships and environmental influences on the tiger prawn (\emph{Penaeus esculentus}) fishery in Exmouth Gulf, Western Australia.  Australian Journal of Marine and Freslzwater Research 37:491-505.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(TPrawnsEG)
#' str(TPrawnsEG)
#' head(TPrawnsEG)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=TPrawnsEG,type="l")
#' plot(recruits~stock,data=TPrawnsEG)
#' par(op)
#' 
#' 
NULL