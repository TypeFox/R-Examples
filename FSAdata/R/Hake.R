#' @title Stock and recruitment data for Hake, 1982-1996.
#' 
#' @description Stock and recruitment data for Hake (\emph{Merluccius merluccius}), 1982-1996.
#' 
#' @name Hake
#' 
#' @docType data
#' 
#' @format A data frame with 15 observations on the following 3 variables.
#'  \describe{ 
#'    \item{year}{a numeric vector of years 1982-1996}
#'    \item{recruits}{a numeric vector of the number of recruits in millions} 
#'    \item{spawn.biomass}{a numeric vector of spawning biomass in thousand tonnes} 
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
#' @source Cadima, E. 2003.  Fish Stock Assessment Manual, FAO Fisheries Department.  131 pp.
#' 
#' @keywords datasets
#' 
#' @examples
#' data(Hake)
#' str(Hake)
#' head(Hake)
#' op <- par(mfrow=c(1,2),pch=19)
#' plot(recruits~year,data=Hake,type="l")
#' plot(recruits~spawn.biomass,data=Hake)
#' par(op)
#' 
NULL