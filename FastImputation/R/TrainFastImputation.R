#' Learn from the training data so that later you can fill in missing data
#'
#' Like Amelia, FastImputation assumes that the columns of the data are
#' multivariate normal or can be transformed into approximately
#' multivariate normal.
#' 
#' @param x Dataframe containing training data. Can have incomplete rows.
#' @param constraints A list of constraints.  See the examples below for formatting details.
#' @return An object of class 'FastImputationPatterns' that contains
#'   information needed later to impute on a single row.
#' @export
#' @seealso \code{\link{FastImputation}}
#' @references
#' \url{http://gking.harvard.edu/amelia/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#'
#' data(FItrain)   # provides FItrain dataset
#' patterns <- TrainFastImputation(FItrain)
#'
#' patterns.with.constraints <- TrainFastImputation(
#'   FItrain,
#'   constraints=list(list(1, list(set=0:1)),
#'                    list(2, list(lower=0)),
#'                    list(3, list(lower=0)),
#'                    list(4, list(lower=0)),
#'                    list(5, list(lower=0)),
#'                    list(6, list(lower=0, upper=1)),
#'                    list(7, list(lower=0)),
#'                    list(8, list(lower=0))))
TrainFastImputation <-
function(
  x,
  constraints=list()
) {
  # TODO:
  # - add idvars parameter such that: a vector of column numbers or 
  #   column names that indicates identification variables.  These 
  #   will be dropped from the analysis but copied into the imputed 
  #   datasets.

  if( "data.frame" != class(x) ) stop("Training data must be in a data.frame")
  
  x <- UnfactorColumns(x)  # unfactor the columns
  
  # Fill the constraints so there is a constraint entry for each column
  if( 0==length(constraints) ) {
    filled.constraints <- replicate(ncol(x), list())
  } else {
    filled.constraints <- sapply(1:ncol(x), function(i.col) {
      is.each.constraint.for.this.col <- sapply(constraints, function(this.cons) {
        return( this.cons[[1]] == i.col )
      })

      if( 0 == sum(is.each.constraint.for.this.col) ) {
        return( list() )
      } else if( 1 == sum(is.each.constraint.for.this.col) ) {
        return( constraints[[which(is.each.constraint.for.this.col)]][[2]] )
      } else {
        return( constraints[[max(which(is.each.constraint.for.this.col))]][[2]] )
      }
    })
  }
  
  # Tally the columns with each type of constraint
  cols.bound.to.intervals <- which(sapply(filled.constraints, function(this.cons) 
    !(is.null(this.cons$upper) && is.null(this.cons$lower))))
  cols.bound.to.sets <- which(sapply(filled.constraints, function(this.cons) 
    !is.null(this.cons$set) ))
  
  # Normalize variables bounded to an interval
  for(this.col in cols.bound.to.intervals) {
    x[,this.col] <- NormalizeBoundedVariable(x[,this.col], constraints=filled.constraints[[this.col]])
  }
  
  FastImputationMeans <- colMeans(x, na.rm=TRUE)
  FastImputationCovariance <- CovarianceWithMissing(x)
  
  patterns <- list(FI.means=FastImputationMeans, 
    FI.covariance=FastImputationCovariance, 
    FI.constraints=filled.constraints, 
    FI.cols.bound.to.intervals=cols.bound.to.intervals,
    FI.cols.bound.to.sets=cols.bound.to.sets)
  class(patterns) <- "FastImputationPatterns"
  return( patterns )
}
