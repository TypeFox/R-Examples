#' Convert a list to a dataframe.
#' 
#' Convert a list of different time series to a dataframe. Usually the list is the output of
#' \code{extractPeriod}
#' NOTE: Since it's dataframe, so the dataframes in the input datalist should have the same 
#' date, if not, please use  \code{extractPeriod} to process.
#'
#' @param datalist A list containing different time series, each sub list has to have the same length.
#' @return The converted dataframe
#' 
#' @examples
#' # open file attached in the package.
#' file <- system.file("extdata", "testdl.txt", package = "hyfo")
#' datalist <- dget(file) # read list file.
#' datalist_new <- extractPeriod(datalist, commonPeriod = TRUE)
#' 
#' dataframe <- list2Dataframe(datalist_new)
#' 
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' @export
list2Dataframe <- function(datalist) {
  
  data <- lapply(datalist, function(x) x[, 2:ncol(x)])
  names <- lapply(datalist, function(x) colnames(x)[2:ncol(x)])
  names <- do.call('cbind', names)
  Date <- datalist[[1]][, 1]
  data <- data.frame(data)
  colnames(data) <- names
  data <- data.frame(cbind(Date, data))
  
  return(data)
}