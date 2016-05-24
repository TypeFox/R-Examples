#' @name rfPermute
#' @title Estimate Permutation p-values for Random Forest Importance Metrics.
#' @description Estimate significance of importance metrics for a Random Forest 
#'   model by permuting the response variable. Produces null distribution of 
#'   importance metrics for each predictor variable and p-value of observed.
#'
#' @param x,y,formula,data,subset,na.action,\dots See \code{\link{randomForest}} 
#'   for definitions.
#' @param nrep Number of permutation replicates to run to construct 
#'   null distribution and calculate p-values (default = 100).
#' @param num.cores Number of CPUs to distribute permutation results over.
#'
#' @details All other parameters are as defined in \code{randomForest.formula}. 
#'   A Random Forest model is first created as normal to calculate the observed 
#'   values of variable importance. The response variable is then permuted 
#'   \code{nrep} times, with a new Random Forest model built for each 
#'   permutation step. 
#'
#' @return An \code{rfPermute} object which contains all of the components of a 
#'   \code{randomForest} object plus:
#'   \item{null.dist}{A list containing two three-dimensional arrays of null 
#'     distributions for \code{unscaled} and \code{scaled} importance measures.} 
#'   \item{pval}{A three dimensional array containing permutation p-values for 
#'     \code{unscaled} and \code{scaled} importance measures.}
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @keywords tree classif regression
#' 
#' @seealso 
#' \code{\link{plot.rfPermute}} for plotting null distributions from the \code{rfPermute} objects. \cr
#' \code{\link{rp.importance}} for extracting importance measures. \cr
#' \code{\link{rp.combine}} for combining multiple \code{rfPermute} objects.\cr
#' \code{\link{randomForest}}
#'
#' @examples
#'   # A regression model using the ozone example
#'   data(airquality)
#'   ozone.rfP <- rfPermute(Ozone ~ ., data = airquality, ntree = 500, na.action = na.omit, nrep = 50)
#'   
#'   # Plot the null distributions and observed values.
#'   layout(matrix(1:6, nrow = 2))
#'   plot(ozone.rfP) 
#'   layout(matrix(1))
#'   
#'   # Plot the unscaled importance distributions and highlight significant predictors
#'   plot(rp.importance(ozone.rfP, scale = FALSE))
#'   
#'   # ... and the scaled measures
#'   plot(rp.importance(ozone.rfP, scale = TRUE))
#'
#'
#' @importFrom randomForest randomForest
#' @export
#' 
rfPermute <- function(x, ...) UseMethod("rfPermute")
