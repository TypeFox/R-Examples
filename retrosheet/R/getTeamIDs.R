#' @name getTeamIDs
#'
#' @title Retrieve team IDs for event files
#'
#' @description This function retrieves the team ID needed for the
#' \code{team} argument of \code{getRetrosheet("play", year, team)}.
#'
#' @param year A single valid four-digit numeric year.
#' @param quiet logical. Passed to \code{\link[utils]{download.file}}
#' @param ... further arguments passed to \code{\link[utils]{download.file}}.
#'
#' @return If the file exists, a named vector of IDs for the given year.
#' Otherwise \code{NA}.
#'
#' @details All currently available years can be retrieved with
#' \code{type.convert(substr(getFileNames()$event, 1L, 4L))}
#'
#' @importFrom RCurl url.exists
#' @importFrom data.table fread
#' @export
#'
#' @examples getTeamIDs(2010)
#'
getTeamIDs <- function(year, quiet = TRUE, ...) {
    stopifnot(is.numeric(year), length(year) == 1L)
    path <- sprintf("http://www.retrosheet.org/events/%deve.zip", year)
    if (url.exists(path)) {
        tmp <- tempfile()
        on.exit(unlink(tmp))
        download.file(path, destfile = tmp, quiet = quiet, ...)
    } else {
        available <- grep(year, getFileNames()$event)
        if(!length(available)) {
            return(NA)
        }
    }
    fname <- paste0("TEAM", year)
    unzip(tmp, files = fname)
    on.exit(unlink(fname), add = TRUE)
    read <- suppressWarnings(fread(fname, header = FALSE, drop = 2:3))
    out <- structure(read[[1L]], .Names = read[[2L]])
    out
}
