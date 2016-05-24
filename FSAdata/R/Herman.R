#' @title Lengths for Walleye, Yellow Perch, Black Crappie, and Black Bullheads from Lake Herman, SD.
#' 
#' @description Total lengths of Walleye (\emph{Sander vitreus}), Yellow Perch (\emph{Perca flavescens}), Black Crappie (\emph{Pomoxis nigromaculatus}), and Black Bullheads (\emph{Ameiurus melas}) for four years in Lake Herman, SD.
#' 
#' @details Lake Herman was sampled on June 20-22, 2005 with four overnight gillnet sets and 10 overnight trapnet sets. The trapnets were constructed with 19-mm (0.75 in) bar-mesh netting, 0.9 m high x 1.5 m wide (3 ft high x 5 ft wide) frames and 18.3 m (60 ft) long leads. The gillnets were 45.7 m long x 1.8 m deep (150 ft long x 6 ft deep) with one 7.6 m (25 ft) panel each of 13, 19, 25, 32, 38 and 51-mm (0.5, 0.75, 1, 1.25, 1.5, and 2 in) bar-mesh monofilament netting.
#' 
#' @name Herman
#' 
#' @docType data
#' 
#' @format A data frame of 5931 observations on the following 3 variables:
#'  \describe{
#'    \item{tl}{Total lengths (cm).} 
#'    \item{spec}{Species codes (\code{wae}=walleye, \code{yep}=yellow perch, \code{bkc}=black crappie, and \code{bbh}=black bullhead).} 
#'    \item{yr}{Capture years.}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Length Frequency
#'    \item Size Structure
#'    \item PSD
#'  }
#'  
#' @concept 'Length Frequency' 'Size Structure' PSD
#' 
#' @source From a South Dakota Fish and Game report that was (does not appear to be there (or anywhere) now) at http://www.sdgfp.info/Wildlife/fishing/SELakes/Herman05.pdf.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Herman)
#' str(Herman)
#' head(Herman)
#' op <- par(mfrow=c(2,2),pch=19)
#' ### Four (of 16 possible) examples
#' with(subset(Herman,spec=="bbh" & yr==2003),hist(tl,main="Black Bullhead, 2003"))
#' with(subset(Herman,spec=="bkc" & yr==2001),hist(tl,main="Black Crappie, 2001"))
#' with(subset(Herman,spec=="yep" & yr==2003),hist(tl,main="Yellow Perch, 2003"))
#' with(subset(Herman,spec=="wae" & yr==1999),hist(tl,main="Walleye, 1999"))
#' par(op)
#' 
NULL