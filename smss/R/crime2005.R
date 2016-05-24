##' @name crime2005
##' @title 2005 Statewise Crime
##' @description Data from the Statistical Abstract of the United States
##' 
##' @format \Sexpr[stage=build,results=rd]{data(crime2005); smss:::describe_df(crime2005)}
##' \describe{
##' \item{\code{STATE}}{U.S. state}
##' \item{\code{VI}}{violent crime rate (number of violent crimes per 100,000 population)}
##' \item{\code{VI2}}{violent crime rate (number of violent crimes per 10,000 population)}
##' \item{\code{MU}}{murder rate}
##' \item{\code{ME}}{percent in metropolitan areas}
##' \item{\code{WH}}{percent white}
##' \item{\code{HS}}{percent high school graduates}
##' \item{\code{PO}}{percent below the poverty level.}
##' }
##' @examples
##' data(crime2005)
##' summary(crime2005)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL
