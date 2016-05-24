#' Mapped plant community time series, Dubois, ID
#'
#' A list with two datasets containing presence-absence and environmental data
#' for a plant community of sagebrush steppe in Dubois, Idaho, USA
#'
#' A historical dataset consisting of a series of permanent 1-\eqn{m^2} quadrats
#' located on the sagebrush steppe in eastern Idaho, USA, between 1923 and 1973.
#' It also contains records of monthly precipitation, mean temperature and
#' snowfall. Total precipitation, total snowfall, and mean annual temperature
#' have been calculated from the original data.
#'
#' @format A list with 2 dataframes, one corresponding to the presence-absence
#'   data and the other to the environmental variables. The first dataframe has
#'   in columns: \describe{ \item{quad}{Name of the quadrat surveyed}
#'   \item{species}{Name of the species found} \item{Presence-absence
#'   data}{Several columns with the year in which the surveys were conducted}}
#'   The second dataframe has the following columns: \describe{ \item{YEAR}{Year
#'   in which surveys were conducted} \item{Environmental variables}{Data of the
#'   recorded environmental variables in the form XXX.YYY, where XXX denotes a
#'   month (or a total) and YYY can refer to snow (in inches), temperature
#'   (fahrenheit degrees) or precipitation (in inches)}}
#'
#' @note Only quadrats Q1, Q2, Q3, Q4, Q5, Q6, Q25 and Q26 are included here.
#'   The surveys were conducted annually from 1932 to 1955 with some gaps for
#'   the quadrats included here.
#'
#' @source
#' \url{https://knb.ecoinformatics.org/#view/doi:10.5063/AA/lzachmann.6.36}
#'
#' @name idaho
NULL
