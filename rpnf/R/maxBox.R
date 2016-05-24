#' Returns the maximum box number in given column
#'
#' @param redData Data to consider
#' @param column Column to consider
maxBox <- function(redData,column) {
  # TODO we have to make a difference between X and O columns!
  # max will be inaccurate for O columns
  max(redData$boxnumber[redData$column==column])
}
