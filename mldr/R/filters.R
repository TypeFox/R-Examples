#' @title Filter rows in a\code{mldr} returning a new \code{mldr}
#' @description Generates a new \code{mldr} object containing the selected
#' rows from an existent \code{mldr}
#' @param mldrObject Original \code{mldr} object from which some rows are going to be selected
#' @param rowFilter Expression to filter the rows
#' @return A new \code{mldr} object with the selected rows
#' @seealso \code{\link{mldr_from_dataframe}}, \code{\link{==.mldr}}, \code{\link{+.mldr}}
#' @examples
#'
#' library(mldr)
#'
#' highlycoupled <- genbase[.SCUMBLE > 0.05] # Select instances with highly imbalanced coupled labels
#' summary(highlycoupled)   # Compare the selected instances
#' summary(genbase)         # with the traits of the original MLD
#'
#' @export

"[.mldr" <- function(mldrObject, rowFilter = T) {
  rowFilter <- substitute(rowFilter)
  rows <- eval(rowFilter, mldrObject$dataset, parent.frame())
  newDataset <- mldrObject$dataset[rows, 1:mldrObject$measures$num.attributes]

  mldr_from_dataframe(newDataset, labelIndices = mldrObject$labels$index, name = mldrObject$name)
}
