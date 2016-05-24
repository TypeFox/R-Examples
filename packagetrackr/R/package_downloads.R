#' packagetracker
#'
#' @name packagetracker
#' @docType package
#' @import magrittr
NULL
globalVariables(c("package", "content"))

#' Track package downloads from Rstudio's CRAN mirror
#'
#' Results are cached in a local folder.
#'
#' @param package_name Name of the package to get download statistics for
#' @param start first day of requested download stats
#' @param end last day of requested download stats
#' @param cache_dir cache folder to use, defaults to one given by rappdirs
#' @param force if TRUE, user is not prompted to confirm writing to hard
#' disk (intended for non-interactive use)
#' \code{\link{package_download_cache_dir}} (via rappdirs)
#'
#' @importFrom dplyr bind_rows arrange
#' @importFrom utils read.csv
#' @examples \dontrun{package_download("package_downloads", start = as.Date("2015-09-01"))}
#' @export
package_downloads <- function(package_name,
                              start = as.Date("2012-10-01"),
                              end = Sys.Date()-1,
                              cache_dir = package_download_cache_dir(package_name),
                              force = FALSE) {

  if (!force && "y" != readline(paste("Please confirm that package 'packagetrackr'",
                                      "will write downloaded cache data to you home",
                                      "folder under", cache_dir, "[y/N]:"))) {
    return()
  }

  seq(start, end, by = 'day') %>%
    as.character() %>%
    setdiff(dates_in_cache(cache_dir)) %>%
    lapply(download_log, package_name = package_name, cache_dir = cache_dir) %>%
    dplyr::bind_rows() %>%
    dplyr::bind_rows(list.files(cache_dir, full.names = TRUE) %>%
                       lapply(utils::read.csv, stringsAsFactors = HELLNO) %>%
                       dplyr::bind_rows())

}

#' Download and unzip a package download log .csv
#'
#' Format is assumed to be as provided by cran-logs.rstudio.com
#'
#' @param day day of the log to download
#' @param package_name package name if not NULL, downloads are directly filtered for this package
#' @param cache_dir if not NULL, results are cached in this directory
#'
#' @importFrom httr GET
#' @importFrom utils read.csv write.csv
#'
#' @export
download_log <- function(day, package_name = NULL, cache_dir = NULL) {

  tmp_file <- tempfile()

  day %>%
    to_log_url() %>%
    httr::GET() %$%
    writeBin(content, tmp_file)

  tmp_file %>%
    gunzip() %>%
    textConnection() %>%
    utils::read.csv(stringsAsFactors = HELLNO) %>%
    subset(is.null(package_name) | package == package_name)                -> df

  if(!is.null(cache_dir)) {
    utils::write.csv(df, file.path(cache_dir, paste0(day, ".csv")))
  }

  unlink(tmp_file)

  return(df)

}

gunzip <- function(filename) {
  system(paste("gunzip", "-c", filename), intern = TRUE) %>%  ## This will be Linux only
    paste(collapse = "\n")
}

dates_in_cache <- function(cache_dir) {
  list.files(cache_dir) %>%
    tools::file_path_sans_ext() %>%
    basename()
}

to_log_url <- function(days, logs_root = "http://cran-logs.rstudio.com/") {

  year <- as.POSIXlt(days)$year + 1900

  paste0(logs_root,  year, "/", days, ".csv.gz")

}

#' Canonical download directory for package download logs
#'
#' @param package_name name of package
#'
#' @importFrom rappdirs app_dir
#'
#' @export
package_download_cache_dir <- function(package_name) {
  the_dir <- rappdirs::app_dir(paste0("R_package_downloads"))$cache() %>%
    file.path(package_name)
  if (!file.exists(the_dir)) {
    dir.create(the_dir, recursive = TRUE)
  }
  return(the_dir)
}

#' @rdname package_download_cache_dir
#' @export
remove_package_download_cache_dir <- function(package_name) {
  package_name %>%
    package_download_cache_dir() %>%
    unlink(recursive = TRUE)
}

HELLNO <- FALSE

