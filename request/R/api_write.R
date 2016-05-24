#' Write helper
#'
#' @export
#' @param .data Result of a call to \code{api}
#' @param file (character) Full file path to write to
#' @param overwrite (logical) Will only overwrite existing path if \code{TRUE}
#' @param ... ignored for now
#' @examples \dontrun{
#' ## write to disk
#' ff <- tempfile(fileext = ".json")
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues) %>%
#'   api_write(ff)
#' jsonlite::fromJSON(ff)
#' }
api_write <- function(.data, file, overwrite = FALSE, ...){
  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  modifyList(.data, list(
    write = write_disk(path = file, overwrite = overwrite))
  )
}
