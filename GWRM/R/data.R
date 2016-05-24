#' Goals scored by footballers in the first division of the Spanish league 
#' 
#' The response variable \code{goals}, is the number of goals scored by the footballers (excluding goalkeepers) in the first division of the Spanish league from the 2000/2001 to the 2006/2007 seasons. Since there are footballers who played more than one season, the season in which each one has played more matches has been selected. The covariates considered are the final classification of the team in each season, the position in the field (forward, midfielder and defender) and the number of matches played.
#' 
#' \itemize{
#'   \item{\code{clasif}} { a numeric vector}
#'   \item{\code{position}} { a factor with levels \code{Defender} \code{Forward} \code{Midfielder}} 
#'   \item{\code{played}} { a numeric vector}
#'   \item{\code{goals}} { a numeric vector}
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name goals
#' @usage data(goals)
#' @source MARCA sports paper
#' @format A data frame with 1224 observations on the following 4 variables
#' @references Rodriguez-Avi, J., Conde-Sanchez, A., Saez-Castillo, A. J., Olmo-Jimenez, M. J. and Martinez Rodriguez, A. M.(2009). A generalized Waring regression model for count data. Computational Statistics and Data Analysis, 53, pp. 3717-3725.
NULL