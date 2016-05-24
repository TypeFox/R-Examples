#' @name CountyMonthBirthRate
#' @aliases CountyMonthBirthRate2005Version CountyMonthBirthRate2014Version 
#' @docType data
#' @title Monthly Growth Fertility Rates (GFR) for 12 urban Oklahoma counties
#' @description Monthly Growth Fertility Rates (GFR) for 12 urban counties in Oklahoma
#' between January 1990 and December 1999.  The GFR is defined as the number of births divided 
#' by the number of females (ages 15-44), multiplied by 1,000.
#' 
#' There are two datasets in this package that are almost identical.  The 2014 version is better suited for substantive researchers
#' in the areas of fertility and traumatic cultural events.  The 2005 version recreates the 2005 article
#' and, therefore is better suited for the graphical aims of the 2014 manuscript.
#' 
#' The difference is that the 2005 version uses constant estimate for a county population 
#' --specifically the US Census 1990 estimates.  The 2014 version uses different estimates
#' for each month --specificallly the US intercensal annual estimates, with linear interpolation
#' for February through December of each year.
#' 
#' @format A data frame with 1,440 observations on the following 11 variables.
#' \describe{
#'    \item{Fips}{The county's 5-digit value according to the \emph{F}ederal \emph{I}nformation \emph{P}rocessing \emph{S}tandards.  \code{integer}}
#'    \item{CountyName}{The lower case name of the county. \code{character}} 
#'    \item{Year}{The year of the record, ranging from 1990 to 1999. \code{integer}}
#'    \item{Month}{The month of the record, ranging from 1 to 12. \code{integer}} 
#'    \item{FecundPopulation}{The number of females in the county, ages of 
#'    15 to 44. \code{numeric}} 
#'    \item{BirthCount}{The number of births in a county for the given month. 
#'    \code{integer}}
#'    \item{Date}{The year and month of the record, with a date of the 15th. Centering the date within the month makes the value a little more representative and the graphs a little easier. \code{date}} 
#'    \item{DaysInMonth}{The number of days in the specific month. \code{integer}} 
#'    \item{DaysInYear}{The number of days in the specific years \code{integer}} 
#'    \item{StageID}{The `Stage' of the month.  The pre-bombing records are `1' (accounting for 9 months of gestation); the post-bombing months are `2'. \code{integer}} 
#'    \item{BirthRate}{The Growth Fertility Rate (GFR). \code{numeric}} 
#' }
#' @details 
#' <<Joe, can you please finish/edit this sentence?>>
#' The monthly birth counts were copied from county records by Ronnie Coleman during the 
#' summer of 2001 from state vital statistics records.  It was collected
#' for \href{http://www.ncbi.nlm.nih.gov/pubmed/16463916}{Rodgers, St. John, & Coleman (2005)}.
#' 
#' The US Census' intercensal estimates are used for the January values of
#' \code{FecundPopluation}.  Values for February-December are interpolated using
#' \href{http://stat.ethz.ch/R-manual/R-devel/library/stats/html/approxfun.html}{\code{approx}}.
#' 
#' The datasets were manipulated to produce this data frame by the two R files
#' \href{https://github.com/OuhscBbmc/Wats/blob/master/UtilityScripts/IsolateCensusPopsForGfr.R}{IsolateCensusPopsForGfr.R}
#' and \href{https://github.com/OuhscBbmc/Wats/blob/master/UtilityScripts/CalculateGfr.R}{CalculateGfr.R}.
#' 
#' @author Will Beasley
#' @references
#' Rodgers, J. L., St. John, C. A. & Coleman R. (2005).  
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/16463916}{Did Fertility Go Up after the Oklahoma City Bombing?  An Analysis of Births in Metropolitan Counties in Oklahoma, 1990-1999.}  
#' \emph{Demography, 42}, 675-692.
#' 
#' \href{http://www.census.gov/popest/data/intercensal/st-co/characteristics.html}{Intercensal
#' estimates for 199x.}
#' 
#' \href{http://www.census.gov/popest/data/intercensal/county/county2010.html}{Intercensal
#' estimates for 200x.}
#' 
#' @keywords datasets
#' @examples 
#' library(ggplot2) 
#' 
#' ##2005 Version (see description above)
#' ds2005 <- CountyMonthBirthRate2005Version
#' ggplot(ds2005, aes(x=Date, y=BirthRate, color=factor(Fips))) + 
#' geom_line() +
#' labs(title="County Fertility - Longitudinal") 
#' 
#' ggplot(ds2005, aes(x=BirthRate, color=factor(Fips))) + 
#' geom_density() +
#' labs(title="Distributions of County Fertility")
#' 
#' ##2014 Version (see description above)
#' ds2014 <- CountyMonthBirthRate2014Version
#' ggplot(ds2014, aes(x=Date, y=BirthRate, color=factor(Fips))) + 
#' geom_line() +
#' labs(title="County Fertility - Longitudinal") 
#' 
#' ggplot(ds2014, aes(x=BirthRate, color=factor(Fips))) + 
#' geom_density() +
#' labs(title="Distributions of County Fertility")
NULL