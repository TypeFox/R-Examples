#' Pantheon dataset
#'
#' Dataframe of globally famous people according to MIT's Pantheon 1.0.
#' Dataset includes the number of globally famous people for a \strong{sample} 
#' of 10 countries and 53 different occupations. 
#' The \strong{complete} dataset is described in [Yu et al., 2015].
#'
#' @format A dataframe with the variables:
#' \describe{
#' \item{Country}{Name of the country}
#' \item{Occupation}{Occupation according to the taxonomy of Pantheon}
#' \item{Value}{ Quantity of globally famous people that were born in that country} 
#' }
#'
#' @source \url{http://pantheon.media.mit.edu/}
#' @references 
#' A. Z. Yu, S. Ronen, K. Hu, T. Lu, and C. A. Hidalgo, 'Pantheon: A Dataset for the Study of Global Cultural Production', arXiv:1502.07310 [physics], Feb. 2015.
#' @keywords dataset
#' @examples
#' data(pantheon)
#' str(pantheon)
#' summary(pantheon)
#' pantheon[pantheon$Country=="Chile",]
"pantheon"

#' Geese dataset
#'
#' A matrix of species of geese.
#' The dataset includes the quantity of 4 species of geese observed by year in the Netherlands.
#' The data comes from the Dutch bird protection organisation Sovon.
#'
#' @format A matrix with the variables:
#' \describe{
#' \item{Columns}{Year of observation}
#' \item{Rows}{Species}
#' }
#'
#' @source \url{https://www.sovon.nl/en}
#' @source \url{http://www.compass-project.eu/applets/3/index_EN.html}
#' @keywords dataset
#' @examples
#' str(geese)
#' summary(geese)
#' geese[,"2000"]
#' geese["Mute Swan",]
"geese"

#' Scidat dataset
#'
#' A matrix of the number of papers authored by 10 countries in 27 areas 
#' of science in the year 2013. Data was retrieved and aggregated from SCimago.
#'
#' @format A matrix with the variables:
#' \describe{
#' \item{Columns}{Areas of Science according SCimago}
#' \item{Rows}{Name of the country}
#' }
#'
#' @source Raw data before the aggregation was queried from 
#' \url{http://www.scimagojr.com/} in 2014.
#' @references SCImago. (2007). SJR-SCImago Journal & Country Rank.
#' @keywords dataset
#' @examples
#' str(scidat)
#' summary(scidat)
#' scidat["United States",]
#' scidat[,"Chemistry"]
"scidat"