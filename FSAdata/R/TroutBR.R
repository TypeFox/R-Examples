#' @title Ages and lengths of migratory Brown and Rainbow Trout.
#' 
#' @description Total lengths (inches) and ages (from scales) of Brown Trout (\emph{Salmo trutta}) and Rainbow Trout (\emph{Oncorhynchus mykiss}) migrating upstream on the Bois Brule River, WI in 1978 and 1979.
#' 
#' @name TroutBR
#' 
#' @docType data
#' 
#' @format A data frame with 851 observations on the following 3 variables:
#'  \describe{
#'    \item{tl}{Measured total length (inches).} 
#'    \item{age}{Assigned age (from scales).} 
#'    \item{species}{Species (\code{Brown} and \code{Rainbow}).} 
#'  }
#'  
#' @section Topic(s):
#'  \itemize{ 
#'    \item Growth
#'    \item von Bertalanffy 
#'  }
#'  
#' @concept Growth 'von Bertalanffy'
#' 
#' @source From tables 18 and 19 of Scholl, D.K., P.J. Peeters, and S.T. Schram.  1984.  Migratory brown trout and rainbow trout populations of the Brule River, Wisconsin.  Wisconsin Department of Natural Resources, Fish Management Report No. 123.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(TroutBR)
#' str(TroutBR)
#' head(TroutBR)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(tl~age,data=TroutBR,subset=species=="Brown",main="Brown Trout")
#' plot(tl~age,data=TroutBR,subset=species=="Rainbow",main="Rainbow Trout")
#' par(op)
#' 
NULL