#' Curl settings
#'
#' @export
#' @family dsl
#' @param .data Result of a call to \code{api}
#' @param ...	Comma separated list of unquoted variable names
#' @examples \dontrun{
#' # Config handler
#' api('http://api.crossref.org/works') %>%
#'  api_config(verbose(), progress()) %>%
#'  peep()
#'
#' xx <- api('http://api.crossref.org') %>%
#'   api_path(works, 10.3897/zookeys.515.9459) %>%
#'   api_config(verbose())
#' }
api_config <- function(.data, ...) {
  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  tmp <- list(...)
  tmp <- if (length(tmp) == 1) tmp[[1]] else do.call('c', tmp)
  modifyList(.data, list(config = tmp))
}
