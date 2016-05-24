#' Function to insert data in time series data
#'
#' @param dataIn as input time series data
#' @param pos as position at with insert new data patch to be insert
#' @param dataInsert as data to be inserted in time series data "dataIn"
#' @return returns the time series data with inserted data patch
#' @export
#========================================================
# Function "insert_patch()" starts here-----------------
#========================================================
insert_patch <- function(dataIn, pos, dataInsert){
  dots <- list(dataInsert)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(dataIn, cumsum(seq_along(dataIn) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  unlist(result)
}
#========================================================
# Function "insert_patch()" ends here-----------------
#========================================================
