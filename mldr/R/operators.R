#' Checks if two mldr objects have the same structure
#'
#' @param mldr1 First \code{mldr} object to compare
#' @param mldr2 Second \code{mldr} object to compare
#' @return \code{TRUE} if the two mldr objects have the same structure, \code{FALSE} otherwise
#' @export
"==.mldr" <- function(mldr1, mldr2) {
  length(mldr1$attributes) == length(mldr2$attributes) &&
    names(mldr1$attributes) == names(mldr2$attributes) &&
    mldr1$attributes == mldr2$attributes
}

#' Generates a new mldr object joining the rows
#' in the two mldrs given as input
#'
#' @param mldr1 First \code{mldr} object to join
#' @param mldr2 Second \code{mldr} object to join
#' @return a new \code{mldr} object with all rows in the two parameters
#' @export
"+.mldr" <- function(mldr1, mldr2) {
  # Check the two mldr's structure
  if(mldr1 == mldr2)
    mldr_from_dataframe(
      rbind(subset(mldr1$dataset, select = 1:mldr1$measures$num.attributes),
            subset(mldr2$dataset, select = 1:mldr1$measures$num.attributes)),
      mldr1$labels$index,
      mldr1$name
    )
  else
    stop(paste(substitute(mldr1), " and ", substitute(mldr2), " don't have the same structure"))
}
