#' @title Catch-at-age for Walleye from eight Kansas reservoirs.
#' 
#' @description Catch-at-age for Walleye (\emph{Sander vitreus}) from eight Kansas reservoirs during 1991-1999.
#' 
#' @note The authors used all age-2 and older Walleye to construct the catch curves.
#' 
#' @name WalleyeKS
#' 
#' @docType data
#' 
#' @format A data frame with 66 observations on the following 3 variables.
#'  \describe{
#'    \item{reservoir}{Reservoir (Cedar.Bluff, Cheney, Glen.Elder, Kirwin, Lovewell, Marion, Webster, Wilson)}
#'    \item{age}{Age of fish at capture}
#'    \item{catch}{Number of fish captured}
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
#' @source From Figure 2 of Quist, M.C., J.L. Stephen, C.S. Guy, and R.D. Schultz.  2004.  Age structure and mortality of Walleyes in Kansas reservoirs: Use of mortality caps to establish realistic management objectives.  North American Journal of Fisheries Management 24:990-1002.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(WalleyeKS)
#' str(WalleyeKS)
#' head(WalleyeKS)
#' op <- par(mfrow=c(2,2),mar=c(3,3,2,1),mgp=c(1.75,0.5,0),tcl=-0.2,pch=19)
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Cedar.Bluff"),main="Cedar Bluff")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Cheney"),main="Cheney")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Glen.Elder"),main="Glen Elder")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Kirwin"),main="Kirwin")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Lovewell"),main="Lovewell")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Marion"),main="Marion")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Webster"),main="Webster")
#' plot(log(catch)~age,data=WalleyeKS,subset=(reservoir=="Wilson"),main="Wilson")
#' par(op)
#' 
NULL