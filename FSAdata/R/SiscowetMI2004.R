#' @title Ages (subsample) and lengths (all fish) for male and female Siscowet Lake Trout captured at four locations in Michigan waters of Lake Superior.
#' 
#' @description Ages (subsample) and lengths (all fish) for male and female Siscowet Lake Trout captured at four locations in Michigan waters of Lake Superior.
#' 
#' @name SiscowetMI2004
#' 
#' @docType data
#' 
#' @format A data frame with 780 observations on the following 8 variables.
#'  \describe{
#'    \item{locID}{Locations (\code{Blind Sucker},\code{Deer Park},\code{Grand Marais},\code{Little Lake Harbor})}
#'    \item{pnldep}{Depth of gillnet panel in which the fish was captured}
#'    \item{mesh}{Gillnet stretch mesh measure}
#'    \item{fishID}{Unique fish identification code}
#'    \item{sex}{Sex (\code{F} and \code{M})}
#'    \item{age}{Assigned ages (yrs; from otoliths)}
#'    \item{len}{Total length (mm)}
#'    \item{wgt}{Weight (g)}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Age-Length Key
#'    \item Growth
#'  }
#'  
#' @concept 'Age-Length Key' 'Growth'
#' 
#' @note Used in the \href{http://derekogle.com/IFAR/}{Introductory Fisheries Analyses with R} book.
#' 
#' @source Obtained directly from the U.S. Fish and Wildlife Service via Michael Seider.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(SiscowetMI2004)
#' str(SiscowetMI2004)
#' head(SiscowetMI2004)
#' xtabs(~age+locID,data=SiscowetMI2004)
#' op <- par(mfrow=c(2,2),pch=19)
#' plot(len~age,data=SiscowetMI2004,subset=locID=="Blind Sucker",main="Blind Sucker")
#' plot(len~age,data=SiscowetMI2004,subset=locID=="Grand Marais",main="Grand Marais")
#' plot(len~age,data=SiscowetMI2004,subset=locID=="Little Lake Harbor",main="Little Lake Harbor")
#' par(op)
#' 
NULL