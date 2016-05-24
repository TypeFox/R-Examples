##' @name birth.rates
##' 
##' @title Birth Rate Data
##'
##' @description Birth rates of several nations. Table 9.13 of 3rd edition.
##'
##' @docType data
##' @format \Sexpr[stage=build,results=rd]{data(birth.rates); smss:::describe_df(birth.rates)}
##' \describe{
##' \item{\code{B}}{crude birth rate (number of births per 1000 population size)}
##' \item{\code{W}}{women's economic activity (female labor force as percent of male)}
##' \item{\code{C}}{percent women using contraception}
##' \item{\code{LI}}{female adult literacy rate}
##' \item{\code{LE}}{female life expectancy}
##' \item{\code{HDI}}{human development index (which has components referring to life expectancy at birth, educational attainment, and income per capita)}
##' \item{\code{GNP}}{gross national product (per capita, in thousands of dollars)}
##' \item{\code{N}}{daily newspaper circulation per 100 people}
##' \item{\code{T}}{number of televisions per 100 people.}
##' }
##'
##' @examples
##' data(birth.rates)
##' summary(birth.rates)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL
