##' @name UNdata
##' @title Data from 2005 UN Human Development Report.
##' 
##' @description
##' These data are from the Human Development report of 2005.  Note however that the
##' values given are from 2003 except for CO2 which is 2002.
##'
##' Table 9.13 in the 4th Edition.
##'
##' @format \Sexpr[stage=build,results=rd]{data(UNdata); smss:::describe_df(UNdata)}
##' \describe{
##' \item{\code{HDI}}{HDI value}
##' \item{\code{Fert}}{Total Fertility rate (births/woman) }
##' \item{\code{Cont}}{Contraceptive prevalence rate (\%) }
##' \item{\code{Cell}}{Cellular subscribers (per 1000 people) }
##' \item{\code{Inter}}{Internet users (per 1000 people) }
##' \item{\code{GDP}}{GDP per capita (US$) }
##' \item{\code{CO2}}{Carbon dioxide emissions per capita (metric tons)}
##' \item{\code{Life}}{Life expectancy at birth, female (years) }
##' \item{\code{Liter}}{Adult literacy rate (female rate \% ages 15 and above)}
##' \item{\code{FemEc}}{Female economic activity rate (\% of male rate, ages 15 and above)}
##' }
##' 
##' @source
##' \url{http://www.stat.ufl.edu/~aa/social/data.html}
##'
##' @examples
##' data(UNdata)
##' summary(UNdata)
##' @docType data
NULL
