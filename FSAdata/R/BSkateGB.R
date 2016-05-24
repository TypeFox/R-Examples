#' @title Stock and recruitment data for Barndoor Skate from Georges Bank, 1966-2007.
#' 
#' @description Stock and recruitment data for Barndoor Skate (\emph{Dipturus laevis}) from Georges Bank for three seasons.
#' 
#' @note Only years within each season where more than one spawner and more than one recruit were captured were recorded.  The authors noted that the Beverton-Holt model could NOT be fit to the winter data.
#' 
#' @name BSkateGB
#' 
#' @docType data
#' 
#' @format A data frame with 31 observations on the following 4 variables.
#'  \describe{
#'    \item{spawners}{a numeric vector giving CPUE of spawning fish}
#'    \item{recruits}{a numeric vector containing the CPUE of recruits}
#'    \item{year}{a numeric vector containing the year of the survey (recruits have been properly lagged (3 years) to match with spawners)}
#'    \item{season}{a factor containing the season of capture (fall, spring, winter)}
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
#' @source estimated (using digitize) from figures 4-6 in Gedamke, T., J.M. Hoenig, W.D. DuPaul, and J.A. Musick.  2009.  Stock-recruitment dynamics and the maximum population growth rate of the barndoor skate on Georges Bank.  North American Journal of Fisheries Management 29:512-526.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(BSkateGB)
#' str(BSkateGB)
#' head(BSkateGB)
#' op <- par(mfrow=c(3,2),pch=19)
#' plot(recruits~year,data=BSkateGB,subset=season=="fall",type="b",main="fall")
#' plot(recruits~spawners,data=BSkateGB,subset=season=="fall",main="fall")
#' plot(recruits~year,data=BSkateGB,subset=season=="spring",type="b",main="spring")
#' plot(recruits~spawners,data=BSkateGB,subset=season=="spring",main="spring")
#' plot(recruits~year,data=BSkateGB,subset=season=="winter",type="b",main="winter")
#' plot(recruits~spawners,data=BSkateGB,subset=season=="winter",main="winter")
#' par(op)
#' 
NULL
