##' @name statewide.crime.2
##' 
##' @title Statewide Crime 2
##' 
##' @description The data are from Statistical Abstract of the United States and most variables were measured in 1993. Table 9.1 of the 3rd edition, and some appear in Table 9.1 of the 4th edition.
##" 
##' @docType data
##' 
##' @format \Sexpr[stage=build,results=rd]{data(statewide.crime.2); smss:::describe_df(statewide.crime.2)}
##' \describe{
##' \item{\code{State}}{U.S. State}
##' \item{\code{VR}}{violent crime rate (per 100,000 people in population)}
##' \item{\code{MR}}{murder rate (per 100,000 people in population)}
##' \item{\code{M}}{percent in metropolitan areas}
##' \item{\code{W}}{percent white}
##' \item{\code{H}}{percent high school graduates}
##' \item{\code{P}}{percent with income below the poverty level}
##' \item{\code{S}}{percent of families headed by a single parent. }
##' }
##'
##' @examples
##' data(statewide.crime.2)
##' summary(statewide.crime.2)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL
