#' Read in simulated data
#' 
#' Functions to read in simulation data from folder. Can be csv or RData files.
#' 
#' @param path path to the files you want to read in.
#' @param ... arguments passed to \code{\link{read.csv}}
#' @param returnList if \code{TRUE} a list containing the data.frames. Very much
#'   like the output of sim. If \code{FALSE} a single data.frame is returned,
#'   using \code{\link[dplyr]{rbind_all}}
#' 
#' @rdname read_data
#' @export
sim_read_data <- function(path, ..., returnList = FALSE) {
  filesToLoad <- dir(path = path, pattern = "*.csv$", full.names = TRUE)
  dataList <- lapply_remove_null(filesToLoad, file_reader, ...)
  if (returnList) dataList else dplyr::bind_rows(dataList)
}

#' @rdname read_data
#' @export
sim_clear_data <- function(path, ...) {
  # To remove non readable files from folder
  filesToLoad <- dir(path = path, pattern = "*.csv$", full.names = TRUE)
  dataList <- lapply(filesToLoad, file_reader, ...)
  file.remove(filesToLoad[sapply(dataList, is.null)])
}

file_reader <- function(file, ...) {
  # Sometimes strange files may be in a folder - for example they are empty
  out <- try({
    read.csv(file = file, ...)
  })
  if (inherits(out, "try-error")) NULL else out
}

lapply_remove_null <- function(l, f, ...) {
  outList <- lapply(l, f, ...)
  outList[!sapply(outList, is.null)]
}

#' @rdname read_data
#' @export
sim_read_list <- function(path) {
  # this fun reads all RData files and returns its results.
  filesToLoad <- dir(path = path, pattern = "*.RData$", full.names = TRUE)
  lapply_remove_null(filesToLoad, file_reader_rdata)
}

file_reader_rdata <- function(filename) {
  res <- NULL
  try({
    load(file = filename, envir = environment())
  })
  res
}

#' @rdname read_data
#' @export
sim_clear_list <- function(path) {
  # To remove non readable files from folder
  filesToLoad <- dir(path = path, pattern = "*.RData$", full.names = TRUE)
  dataList <- lapply(filesToLoad, file_reader_rdata)
  file.remove(filesToLoad[sapply(dataList, is.null)])
}
