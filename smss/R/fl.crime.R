##' @name fl.crime
##' @title Florida Crime Data
##'
##' @description Crime data for florida counties. Table 9.16 or 9.17 from the 4th edition.
##' Source: Dr. Larry Winner, University of Florida.
##'
##' @format \Sexpr[stage=build,results=rd]{data(fl.crime); smss:::describe_df(fl.crime)}
##' \describe{
##'    \item{\code{County}}{county name}
##'    \item{\code{C}}{crime rate}
##'    \item{\code{I}}{median income}
##'    \item{\code{HS}}{percent completing high school}
##'    \item{\code{U}}{percent urban}
##' }
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
##' @examples
##' data(fl.crime)
##' summary(fl.crime)
##' @docType data
NULL
