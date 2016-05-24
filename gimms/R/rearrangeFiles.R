#' Rearrange GIMMS files by date
#'
#' @description
#' Rearrange GIMMS-related files in ascending order of time.
#'
#' @param x Character. Vector of (local or online) filepaths. If \code{NULL},
#' \code{dsn} will be searched for available files via pattern matching.
#' @param dsn Character. Path to look for GIMMS data. If not supplied and 'x' is
#' missing, this defaults to the current working directory.
#' @param pattern Character. A regular expression passed on to
#' \code{\link{list.files}}.
#' @param pos Integer. The start positions of year, month and period ('a' or
#' 'b') in the target GIMMS files. Unless modified, this usually defaults to
#' \code{c(4, 6, 11)} (see 'References').
#' @param ... Further arguments passed on to \code{\link{list.files}}.
#'
#' @return
#' A vector of filepaths arranged in ascending order of time.
#'
#' @author
#' Florian Detsch
#'
#' @references
#' \url{http://ecocast.arc.nasa.gov/data/pub/gimms/3g.v0/00READMEgeo.txt}
#' (accessed on January 15, 2016).
#'
#' @seealso
#' \code{\link{list.files}}
#'
#' @examples
#' ## latest version of files inventory
#' gimms_files <- updateInventory(sort = FALSE)
#' head(gimms_files)
#'
#' ## re-arrange vector with available files according to date
#' gimms_files_arr <- rearrangeFiles(gimms_files)
#' head(gimms_files_arr)
#'
#' @export rearrangeFiles
#' @name rearrangeFiles
rearrangeFiles <- function(x = NULL,
                           dsn = getwd(),
                           pattern = "^geo.*.VI3g$",
                           pos = c(4, 6, 11),
                           ...) {

  if (length(pos) != 3)
    stop("'pos' must be a vector of length 3 (i.e., start position of year, month and day); see ?rearrangeFiles. \n")

  ## if `is.null(fls)`, apply pattern matching in 'dsn'
  if (is.null(x))
    x <- list.files(dsn, pattern = pattern, ...)

  ## vector to data.frame
  gimms_df <- data.frame(file = x, stringsAsFactors = FALSE)

  ## switch current locale time to us standard
  systime_locale <- Sys.getlocale(category = "LC_TIME")

  if (Sys.info()[["sysname"]] == "Windows") {
    invisible(Sys.setlocale(category = "LC_TIME", locale = "C"))
  } else {
    invisible(Sys.setlocale(category = "LC_TIME", locale = "en_US.UTF-8"))
  }

  ## create columns 'year', 'month' and 'day'
  gimms_df <- transform(gimms_df,
                        "year" = substr(basename(file), pos[1], pos[1] + 1),
                        "month" = substr(basename(file), pos[2], pos[2] + 2),
                        "day" = ifelse(substr(basename(file), pos[3], pos[3]) == "a", 1, 15))

  ## create column 'date'
  gimms_df$date <- as.Date(paste0(gimms_df$day, gimms_df$month, gimms_df$year),
                           format = "%d%b%y")

  ## re-arrange rows by 'date'
  gimms_df <- gimms_df[order(gimms_df$date), ]

  ## revoke locale time adjustment
  Sys.setlocale(category = "LC_TIME", locale = systime_locale)

  ## return rearranged files
  gimms_fls <- gimms_df$file
  return(gimms_fls)
}
