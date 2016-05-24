#' Maddison Project Dataset
#'
#' This package contains the \href{http://www.ggdc.net/maddison/maddison-project/data.htm}{Maddison Project Dataset}
#' with estimates of GDP per capita for all countries in the world between AD 1
#' and 2010 in a format amenable to analysis in R.
#'
#' @format A data frame with nine variables:
#' \describe{
#' \item{\code{year}}{Year of estimate}
#' \item{\code{country_original}}{Country name in the original form used in the data base (e.g. "(Centre-North) Italy" instead of simply "Italy" as in variable \code{country})}
#' \item{\code{gdp_pc}}{Estimated GDP per capita in 1990 international Geary-Khamis dollar}
#' \item{\code{country}}{Full country name}
#' \item{\code{iso2c}}{Country iso2c code}
#' \item{\code{iso3c}}{Country iso3c code}
#' \item{\code{continent}}{Country continent}
#' \item{\code{region}}{Country region}
#' \item{\code{aggregate}}{Whether or not the "country" is an aggregate (e.g. "Total world")}
#' }
#'
#' The database was last updated in January 2013.
#'
#' As per instructions on the Maddison Project website, please site the data as
#' follows:
#'
#'\tabular{ll}{
#' When using the data \tab The Maddison-Project, <http://www.ggdc.net/maddison/maddison-project/home.htm>, 2013 version.\cr
#' When refering to underlying methodology and main results \tab Bolt, J. and J. L. van Zanden (2014). The Maddison Project: collaborative research on historical national  accounts. The Economic History Review, 67 (3): 627-651. When using individual country  data \cr
#' When using individual country data \tab See country-source references in the appendix of Bolt and van Zanden (2014).\cr
#' }
#'
#' The package is not affiliated with, nor endorsed by, the Maddison Project.
"maddison"