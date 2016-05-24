if ( !isGeneric("downloadGimms") ) {
  setGeneric("downloadGimms", function(x, ...)
    standardGeneric("downloadGimms"))
}
#' Download GIMMS 3G data
#'
#' @description
#' Download GIMMS 3G binary data for a given time span from the NASA Ames
#' Ecological Forecasting Lab's FTP server
#' (\url{http://ecocast.arc.nasa.gov/data/pub/gimms/3g.v0/}, accessed on
#' January 15, 2016).
#'
#' @param x If 'Date', start date for download (e.g. "2000-01-01"). If
#' 'numeric', start year for download (e.g. 2000). If 'character', a vector of
#' full online filepath(s) to download, typically returned from
#' \code{\link{updateInventory}}. If not supplied, download will start from the
#' oldest file available.
#' @param y If 'Date', end date for download. If 'numeric', end year for
#' download. If not supplied, download will stop with the latest file available.
#' @param dsn 'character'. Destination folder for file download. If not supplied,
#' all downloaded files will be stored in the current working directory.
#' @param overwrite Logical. If \code{TRUE}, already downloaded files in 'dsn'
#' will be overwritten.
#' @param quiet Logical. If \code{TRUE}, information sent to the console is
#' reduced.
#' @param mode See \code{\link{download.file}}.
#' @param cores Integer. Number of cores for parallel computing.
#' @param ... Further arguments passed on to \code{\link{download.file}}, e.g.
#' 'method'.
#'
#' @return
#' A vector of filepaths.
#'
#' @author
#' Florian Detsch
#'
#' @seealso
#' \code{\link{download.file}}
#'
#' @examples
#' \dontrun{
#' ## Destination folder for data download
#' gimms_dir <- paste0(getwd(), "/data")
#'
#' ## 'Date' method
#' gimms_files_date <- downloadGimms(x = as.Date("2000-01-01"),
#'                                   y = as.Date("2000-06-30"),
#'                                   dsn = gimms_dir)
#'
#' ## 'numeric' method, i.e. full years
#' gimms_files_year <- downloadGimms(x = 2000, y = 2002, dsn = gimms_dir)
#'
#' ## 'character' method, i.e. file names
#' gimms_files <- updateInventory()
#' gimms_files <- gimms_files[grep("geo00", gimms_files)]
#' gimms_files_char <- downloadGimms(x = gimms_files, dsn = gimms_dir)
#'
#' ## 'missing' method, i.e. entire collection
#' gimms_files_full <- downloadGimms(dsn = gimms_dir)
#' }
#'
#' @export downloadGimms
#' @name downloadGimms

################################################################################
### function using 'Date' input ################################################
#' @aliases downloadGimms,Date-method
#' @rdname downloadGimms
setMethod("downloadGimms",
          signature(x = "Date"),
          function(x, y,
                   dsn = getwd(), overwrite = FALSE, quiet = TRUE,
                   mode = "wb", cores = 1L, ...) {

            ## check 'cores'
            cores <- checkCores(cores)

            ## jump to downloadGimms,missing-method if neither 'x' nor 'y' is specified
            if (missing(x) & missing(y))
              downloadGimms(dsn = dsn, overwrite = overwrite, quiet = quiet,
                            mode = mode, cores = cores, ...)

            ## available files
            gimms_fls <- updateInventory()

            ## extract timestamps
            gimms_dts_chr <- monthlyIndices(gimms_fls, timestamp = TRUE,
                                            format = "%Y-%m-%d")
            gimms_dts <- as.Date(gimms_dts_chr)

            ## start (finish) with the first (last) timestamp available if 'x'
            ## ('y') is not specified
            if (missing(x)) x <- gimms_dts[1]
            if (missing(y)) y <- gimms_dts[length(gimms_dts)]

            ## create subset of available files covering user-defined temporal range
            user_dts <- seq(x, y, 1)
            user_dts <- user_dts[which(substr(user_dts, 9, 10) %in% c("01", "15"))]
            gimms_fls <- gimms_fls[gimms_dts %in% user_dts]

            ## download
            gimms_out <- downloader(gimms_fls, dsn = dsn, overwrite = overwrite,
                                    quiet = quiet, mode = mode, cores = cores, ...)

            ## return local filenames
            return(gimms_out)
          })


################################################################################
### function using numeric input (i.e. years) ##################################
#' @aliases downloadGimms,numeric-method
#' @rdname downloadGimms
setMethod("downloadGimms",
          signature(x = "numeric"),
          function(x, y,
                   dsn = getwd(), overwrite = FALSE, quiet = TRUE,
                   mode = "wb", cores = 1L, ...) {

            ## check 'cores'
            cores <- checkCores(cores)

            ## jump to downloadGimms,missing-method if neither 'x' nor 'y' is specified
            if (missing(x) & missing(y))
              downloadGimms(dsn = dsn, overwrite = overwrite, quiet = quiet,
                            mode = mode, cores = cores, ...)

            ## available files
            gimms_fls <- updateInventory()

            ## if specified, subset available files by time frame
            gimms_bsn <- basename(gimms_fls)
            gimms_yrs_chr <- substr(gimms_bsn, 4, 5)
            gimms_yrs_num <- as.numeric(gimms_yrs_chr)

            id_old <- gimms_yrs_num >= 81
            id_new <- gimms_yrs_num <= 80
            gimms_yrs_chr[id_old] <- paste0("19", gimms_yrs_chr[id_old])
            gimms_yrs_chr[id_new] <- paste0("20", gimms_yrs_chr[id_new])
            gimms_yrs_num <- as.numeric(gimms_yrs_chr)

            ## start (finish) with the first (last) year available if 'x' ('y') is not
            ## specified
            if (missing(x)) x <- gimms_yrs_num[1]
            if (missing(y)) y <- gimms_yrs_num[length(gimms_yrs_num)]

            ## subset files
            gimms_fls <- gimms_fls[gimms_yrs_chr %in% seq(x, y)]

            ## download
            gimms_out <- downloader(gimms_fls, dsn = dsn, overwrite = overwrite,
                                    quiet = quiet, mode = mode, cores = cores, ...)

            ## return local filenames
            return(gimms_out)
          })


################################################################################
### function using character input (i.e. files) ################################
#' @aliases downloadGimms,character-method
#' @rdname downloadGimms
setMethod("downloadGimms",
          signature(x = "character"),
          function(x, dsn = getwd(), overwrite = FALSE, quiet = TRUE,
                   mode = "wb", cores = 1L, ...) {

            ## check 'cores'
            cores <- checkCores(cores)

            ## download
            gimms_out <- downloader(x, dsn = dsn, overwrite = overwrite,
                                    quiet = quiet, mode = mode, cores = cores, ...)

            ## return vector with output files
            return(gimms_out)
          })


################################################################################
### function using no input (i.e. download entire collection) ##################
#' @aliases downloadGimms,missing-method
#' @rdname downloadGimms
setMethod("downloadGimms",
          signature(x = "missing"),
          function(dsn = getwd(), overwrite = FALSE, quiet = TRUE,
                   mode = "wb", cores = 1L, ...) {

            ## check 'cores'
            cores <- checkCores(cores)

            ## available files
            gimms_fls <- updateInventory()

            ## download
            gimms_out <- downloader(gimms_fls, dsn = dsn, overwrite = overwrite,
                                    quiet = quiet, mode = mode, cores = cores, ...)

            ## return local filenames
            return(gimms_out)
          })
