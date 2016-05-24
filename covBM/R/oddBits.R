.onUnload <- function(libpath) {
  library.dynam.unload("covBM", libpath)
}

#' Serial CD4 counts in children with HIV.
#'
#' Dataset used in \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}
#' by Andrew Gelman and Jennifer Hill (Cambridge University Press, 2006). Rows with missing
#' values of 'CD4CNT', 'visage' or 'baseage' have been removed.
#'
#' @format A data frame with 976 rows and 11 variables:
#' \describe{
#'   \item{newpid}{Patient ID code.}
#'   \item{t}{Time, in years from first visit.}
#'   \item{sqrtcd4}{Square root of CD4 count.}
#'   \item{treatmnt}{Indicator variable for treatment, 1 represents control group and 2 indicates zinc treatment group.}
#'   \item{CD4CNT}{CD4 count on original (untransformed) scale.}
#'   \item{baseage}{Age of child in years at initial visit.}
#'   \item{visage}{Age of child in years at given visit.}
#' }
#' @source \url{http://www.stat.columbia.edu/~gelman/arm/examples/cd4/}
"cd4"