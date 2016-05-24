#' U.S. Major League Batting Average Data: 2002-2012
#'
#' Data frame consisting of the following variables:
#'
#' Data is aggregated into half seasons: so season indicates whether
#' the observation is in the first or second half of the season of a
#' given year. Only players who have more than 10 at bats in any half
#' season are included, and only players who have more than three
#' half seasons are represented. The transformed batting average is
#' \eqn{arcsin(sqrt((H + 1/4)/(AB + 1/2)))}. Only regular seasons data are included. 
#' R programs to extract the data from the original sources are available on request. 
#' \itemize{
#' 	\item Name
#' 	\item IdNum
#' 	\item Year
#' 	\item Halfseason
#' 	\item Pitcher
#' 	\item HA transformed batting average; 
#' 	\item AB at bats 
#' 	\item H hits 
#' 	\item BB walks 
#' 	\item YOB Year of Birth; 
#' 	\item age age of the player
#' 	\item agesq age squared
#' }
#' @source  ESPN Website: \url{http://espn.go.com/mlb/statistics}
#' @references Gu, Jiaying and Roger Koenker (2015) Empirical Bayesball Remixed,
#' preprint.
#' @keywords data
#' @name bball
NULL
