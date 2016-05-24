##' @name oecd.data
##'
##' @title OECD Data
##'
##' @description The OECD (Organization for Economic Cooperation and Development) is an international
##' organized that consists of developed countries.  This dataset has UN data on OECD countries.
##'
##' See Table 3.11 and Exercise 3.6.
##'
##' @docType data
##'
##' @format \Sexpr[stage=build,results=rd]{data(oecd.data); smss:::describe_df(oecd.data)}
##' \describe{
##' \item{\code{nation}}{Country name}
##' \item{\code{GDP}}{GDP per capita (in US dollars)}
##' \item{\code{Unemploy}}{Percent unemployed}
##' \item{\code{Inequal}}{A measure of inequality that compares the wealth of the richest 10\% to the poorest 10\%}
##' \item{\code{Health}}{Public expenditure on health (as a percent of GDP)}
##' \item{\code{Phys}}{Number of physicians per 100,000 people}
##' \item{\code{CO2}}{Carbon dioxide emissions (per capita, in metric tons)}
##' \item{\code{Parlia}}{The percentage of seats in parliament held by women.}
##' \item{\code{FemEcon}}{Female economic activity as a percentage of the male rate}
##' }
##'
##' @examples
##' data(oecd.data)
##' summary(oecd.data)
##' @source \url{http://www.stat.ufl.edu/~aa/social/data.html}
NULL
