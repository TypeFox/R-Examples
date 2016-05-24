#' @title Add new row to dataframe
#' @description Facilitates insertion of a new row in a data.frame.
#' @param .data The existing data.frame
#' @param newrow The new row to be appended.
#' @param where An integer for the position to add the row, default is at the top.
#'
#' @examples
#' existingDF <- as.data.frame(matrix(seq(20),nrow=5,ncol=4))
#' existingDF
#' r <- 3
#' newrow <- seq(4)
#'
#' insert.row(existingDF, newrow, r)
#'
#' @export
insert.row <- function(.data, newrow, where=1){
  .data[seq(where+1,nrow(.data)+1),] <- .data[seq(where,nrow(.data)),]
  .data[where,] <- newrow
  .data
}
NULL

