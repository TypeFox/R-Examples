#' @title Catch-at-age of Great Slave Lake Whitefish (commercial) by area.
#' 
#' @description Age composition of commercial Whitefish (\emph{Coregonus clupeaformis}) catches for five areas of Great Slave Lake.
#' 
#' @note All data are from summer samples except for Area.IW which is a winter sample
#' 
#' @name WhitefishGSL
#' 
#' @docType data
#' 
#' @format A data frame with 16 observations on the following 6 variables:
#'  \describe{
#'    \item{age}{Assigned ages.} 
#'    \item{area.IE}{Catches for area IE.}
#'    \item{area.II}{Catches for area II.} 
#'    \item{area.IV}{Catches for area IV.}
#'    \item{area.V}{Catches for area V.} 
#'    \item{area.IW}{Catches for area IW.} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Mortality 
#'    \item Catch curve
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#' 
#' @source From Table 19 of Mosenko, R.W., and G. Low.  1980.  Data from the commercial fishery for lake whitefish \emph{Coregonus clupeaformis} (Mitchill), on Great Slave Lake, Northwest Territories, 1979.  Canadian Data Report of Fisheries And Aquatic Sciences, No. 194.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WhitefishGSL)
#' str(WhitefishGSL)
#' head(WhitefishGSL)
#' op <- par(mfrow=c(3,2),pch=19)
#' plot(log(area.IE)~age,data=WhitefishGSL)
#' plot(log(area.II)~age,data=WhitefishGSL)
#' plot(log(area.IV)~age,data=WhitefishGSL)
#' plot(log(area.V)~age,data=WhitefishGSL)
#' plot(log(area.IW)~age,data=WhitefishGSL)
#' par(op)
#' 
#' # can be reshaped to 'long' format with
#' \dontrun{
#' library(reshape)
#' WhitefishGSL1 <- melt(WhitefishGSL,id.vars="age")
#' names(WhitefishGSL1) <- c("age","area","number")
#' }
#' 
NULL