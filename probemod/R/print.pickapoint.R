#' Print Function For Pick-A-Point
#'
#' Print function for objects of class \code{"pickapoint"}
#'
#' @param x An object of class \code{"pickapoint"}.
#' @param \dots Additional arguments (not supported yet).
#'
#' @return none
#'
#' @examples
#' \dontrun{
#' myModel <- lm('dv ~ iv + mod', data=someData)
#' papresults <- pickapoint(myModel, dv='DV', iv='IV', mod='MOD')
#' papresults
#' }
#' @rdname print.pickapoint
#' @export

print.pickapoint <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("Conditional effects of ", x$iv, " on ", x$dv, " at values of ", paste(x$mod,sep=',') ,"\n")
  print(x$outcome)
  if(x$method == 'percentiles'){
    cat('\nValues for quantitative moderators are 10th, 25th, 50th, 75th, and 90th percentiles\n')
  } else if(x$method == 'meansd'){
    cat('\nValues for quantitative moderators are the mean and plus/minus one SD from the mean\n')
  }
  cat('Values for dichotomous moderators are the two values of the moderator')
}
