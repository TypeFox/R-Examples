#' @title Summarized multiple mark-recapture data for YOY walleye.
#' 
#' @description The numbers of young-of-year walleye (\emph{Sander vitreus}) that were captured, found to have previous marks (i.e., recaptured), and were newly marked on several sampling occasions in 1959, 1960, 1961, and 1962 in Pike Lake, Wisconsin.
#' 
#' @name WalleyePL
#' 
#' @docType data
#' 
#' @format A data frame with 33 observations on the following 5 variables:
#'  \describe{
#'    \item{year}{Sampling year}
#'    \item{t}{Sampling occasion within each year}
#'    \item{caught}{Number of walleye captured}
#'    \item{recaptures}{Number of marked walleyes captured}
#'    \item{newmarks}{Number of unmarked walleyes that were captured, marked, and returned to the population}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{
#'    \item Population Size
#'    \item Abundance
#'    \item Mark-Recapture
#'    \item Capture-Recapture
#'    \item Schnabel
#'    \item Schumacher-Eschmeyer
#'  }
#'  
#' @concept Abundance 'Population Size' 'Mark-Recapture' 'Capture-Recapture' 'Schnabel'
#' 
#' @source From table 3 of Mraz, D.  1968.  Recruitment, growth, exploitation, and management of walleyes in a southeastern Wisconsin lake.  Wisconsin Department of Natural Resources Technical Bulletin 40.  38 pages.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyePL)
#' str(WalleyePL)
#' WalleyePL
#' subset(WalleyePL,year==1960)
#' 
NULL