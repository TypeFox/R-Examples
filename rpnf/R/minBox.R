#' Returns the minimum box number in given column
#'
#' @param redData Data to consider
#' @param column Column to consider
minBox <- function(redData,column){
  # TODO we have to make a difference between X and O columns!
  # min will be inaccurate for X columns
  min(redData$boxnumber[redData$column==column])
}
