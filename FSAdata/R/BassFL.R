#' @title Catch-at-age for Suwanee and Largemouth Bass.
#' 
#' @description Catch-at-age for Suwanee (\emph{Micropterus notius}) and Largemouth Bass (\emph{Micropterus salmoides}) collected from several lakes in Florida, 2001-2002.
#' 
#' @name BassFL
#' 
#' @docType data
#' 
#' @format A data frame with 39 observations on the following 5 variables.
#'  \describe{
#'    \item{species}{Species of bass (Suwanee and Largemouth)}
#'    \item{loc}{Location (SantaFe, Wacissa, Withlacoochee, Ochlockonee)}
#'    \item{year}{Year (2001, 2002)} 
#'    \item{num}{Number of fish captured}
#'    \item{age}{Age of fish at capture}
#'  }
#'  
#' @section Topic(s):
#'  \itemize{ 
#'    \item Total mortality 
#'    \item Catch curve 
#'  }
#'  
#' @concept Mortality 'Catch Curve'
#' 
#' @source From Figure 2 of Bonvechio, T.F., M.S. Allen, and R.L. Cailteux.  2005.  Relative Abundance, Growth, and Mortality of Suwannee Bass  in Four Florida Rivers.  North American Journal of Fisheries Management 25:275-283.  [Was (is?) from http://sfrc.ufl.edu/allenlab/publications/2005bonvechio_allen_cailteux.pdf.]
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BassFL)
#' str(BassFL)
#' head(BassFL)
#' op <- par(mfrow=c(3,2),mar=c(3,3,2,1),mgp=c(1.75,0.5,0),tcl=-0.2,pch=19)
#' plot(log(num)~age,data=BassFL,subset=(loc=="SantaFe" & year==2001 & species=="Suwanee"),
#'   ylim=c(0,max(log(num))),main="Suwanee, Santa Fe")
#' points(log(num)~age,data=BassFL,subset=(loc=="SantaFe" & year==2002 & species=="Suwanee"),col="red")
#' legend("topright",legend=c("2001","2002"),col=c("black","red"),pch=19)
#' plot(log(num)~age,data=BassFL,subset=loc=="Wacissa" & year==2002 & species=="Suwanee",
#'  ylim=c(0,max(log(num))),main="Suwanee, Wacissa")
#' plot(log(num)~age,data=BassFL,subset=loc=="Withlacoochee" & year==2002 & species=="Suwanee",
#'  ylim=c(0,max(log(num))),main="Suwanee, Withlacoochee")
#' plot(log(num)~age,data=BassFL,subset=(loc=="SantaFe" & year==2001 & species=="Largemouth"),
#'  ylim=c(0,max(log(num))),main="Largemouth, Santa Fe")
#' points(log(num)~age,data=BassFL,subset=(loc=="SantaFe" & year==2002 & species=="Largemouth"),
#'  col="red")
#' legend("topright",legend=c("2001","2002"),col=c("black","red"),pch=19)
#' plot(log(num)~age,data=BassFL,subset=loc=="Ochlockonee" & year==2001 & species=="Largemouth",
#'  ylim=c(0,max(log(num))),main="Largemouth, Ochlockonee")
#' points(log(num)~age,data=BassFL,subset=(loc=="Ochlockonee" & year==2002 & species=="Largemouth"),
#'  col="red")
#' legend("topright",legend=c("2001","2002"),col=c("black","red"),pch=19)
#' par(op)
#' 
NULL
