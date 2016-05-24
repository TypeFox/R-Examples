#' Update GIMMS 3G file inventory
#'
#' @description
#' Download the latest version of the GIMMS 3G file inventory from the NASA Ames
#' Ecological Forecasting Lab's FTP server
#' (\url{http://ecocast.arc.nasa.gov/data/pub/gimms/3g.v0/}, accessed on
#' January 15, 2016) or, if the server is not accessible, load local version of
#' the file inventory.
#'
#' @param sort 'logical'. Determines whether or not to sort the list of
#' available GIMMS files in ascending order of time via
#' \code{\link{rearrangeFiles}}.
#' @param quiet 'logical'. If TRUE, information sent to the console is reduced.
#'
#' @return
#' A vector of online filepaths.
#'
#' @author
#' Florian Detsch
#'
#' @seealso
#' \code{\link{rearrangeFiles}}.
#'
#' @examples
#' updateInventory()
#'
#' @export updateInventory
#' @name updateInventory
updateInventory <- function(sort = TRUE, quiet = TRUE) {

  ## available files (online)
  if (!quiet)
    cat("Trying to update GIMMS inventory from server...\n")

  gimms_fls <-
    suppressWarnings(
      try(
        readLines("http://ecocast.arc.nasa.gov/data/pub/gimms/3g.v0/00FILE-LIST.txt"),
        silent = TRUE
      )
    )

  ## remove duplicate entries
  gimms_fls <- gimms_fls[!duplicated(basename(gimms_fls))]

  ## available files (offline)
  if (class(gimms_fls) == "try-error") {
    if (!quiet)
      cat("Online update failed. Using local inventory...\n")

    gimms_fls <- readRDS(system.file("extdata", "inventory.rds",
                                     package = "gimms"))
  } else {
    if (!quiet)
      cat("Online update of the GIMMS file inventory successful!\n")
  }

  ## sort files (optional)
  if (sort)
    gimms_fls <- rearrangeFiles(gimms_fls)

  ## return files
  return(gimms_fls)
}
