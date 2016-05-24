#' Performs a Cross Join of Unique combinations
#' 
#' This function makes use of \code{\link{CJ}} function of the data.table package to perform a
#' cross join. The function makes sure that the combinations are unique and removes NAs before
#' joining. doUniqueCJ is rather not used as a standalone function but inside \code{\link{computeShares}}. 
#' 
#' @author Matthias Bannert, Gabriel Bucur
#' @param dt data.table
#' @param cols character vector that denotes names of relevant columns
#' @export
doUniqueCJ <- function(dt, cols) {

unique_values <- lapply(cols, function(x) unique(na.omit(dt[[x]])))
names(unique_values) <- cols

do.call("CJ", unique_values)
}
