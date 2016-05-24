#' @title Stock and recruitment data for Icelandic summer spawning Herring, 1946-1996.
#' 
#' @description Icelandic summer spawning Herring (\emph{Clupeaformis harengus}) stock, recruitment, landings, and fishing mortality by year, 1946-1996.
#' 
#' @name HerringISS
#' 
#' @docType data
#' 
#' @format A data frame of 51 observations on the following 6 variables:
#'  \describe{
#'    \item{year}{Year of data}
#'    \item{ssb}{Spawning stock biomass (tonnes)}
#'    \item{rec}{Recruits -- i.e., 1-year olds (thousands)}
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
#' @source From the Report of the Atlanto-Scandian Herring and Capelin Working Group. ICES Doc. C.M. 1995. Assess: 9.; Anon. 1986. Report of the herring assessment working group for the area south of 62 degrees N. ICES Doc. C.M. 1986. Assess: 19.; Anon. 1991. Report of the Atlanto-Scandian Herring and Capelin Working Group. ICES Doc. C.M. 1991. Assess: 17.; State of marine stocks and environmental conditions in Icelandic waters 1989 Fishing prospects 1990, August 1990. Hafranns\'oknastofnun Fj\"olrit NR. 19. Marine Research Institute, Iceland. Mon Jun 20, 1988.; Report of the Northern Pelagic and Blue Whiting Fisheries Working Group. ICES Doc. C.M. 1997. Assess:14.  Obtained from Ransom Myers online database which was (is?) at http://ram.biology.dal.ca/~myers/data.html.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(HerringISS)
#' str(HerringISS)
#' head(HerringISS)
#' op <- par(mfrow=c(1,2))
#' plot(rec~year,data=HerringISS,type="l")
#' plot(rec~ssb,data=HerringISS)
#' par(op)
#' 
NULL