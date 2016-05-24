#' @title Stock and recruitment data for Pacific Halibut, 1929-1991.
#' 
#' @description Pacific Halibut (\emph{Hippoglossus stenolepis}) stock, recruitment, landings, and fishing mortality by year, 1929-1991.
#' 
#' @name HalibutPAC
#' 
#' @docType data
#' 
#' @format A data frame of 63 observations on the following 5 variables:
#'  \describe{
#'    \item{year}{Year of data}
#'    \item{ssb}{Spawning stock biomass (tonnes)}
#'    \item{rec}{Recruits (thousands)}
#'    \item{land}{Landings (in millions of pounds)}
#'    \item{fmort}{Fishing related mortality}
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
#' @source From Anonymous. 1989. International Pacific Halibut Commission.  Annual Report 1989.; Myhre, Gordon, Peltonen, Stpierre, Skud, and Walden.  International Pacific halibut Commission. Technical Report No. 14.; Pers. Comm. with Pat Sullivan.  Obtained from Ransom Myers online database which was (is?) at http://ram.biology.dal.ca/~myers/data.html.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(HalibutPAC)
#' str(HalibutPAC)
#' head(HalibutPAC)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(rec~year,data=HalibutPAC,type="l")
#' plot(rec~ssb,data=HalibutPAC)
#' par(op)
#' 
NULL