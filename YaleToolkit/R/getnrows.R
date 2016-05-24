#' @title Get the number of rows of the file
#' @description Use iterators to avoid the memory overhead of
#' obtaining the number of rows of a file.
#' @param file the name of a file (possible with a path)
#' @param n the size of the chunks used by the iterator
#' @return an integer
#' @examples
#' data(CO2)
#' write.csv(CO2, "CO2.csv", row.names=FALSE)
#' getnrows("CO2.csv")
#' @export
getnrows <- function(file, n=10000) {
  i <- NULL # To kill off an annoying R CMD check NOTE
  it <- ireadLines(file, n=n)
  return( foreach(i=it, .combine=sum) %do% length(i) )
}