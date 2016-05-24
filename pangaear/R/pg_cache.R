#' @export
#' @rdname pg_data
pg_cache_clear <- function(doi = NULL, prompt = TRUE) {
  if (is.null(doi)) {
    files <- list.files(env$path, full.names = TRUE)
    resp <- if (prompt) readline(sprintf("Sure you want to clear all %s files? [y/n]:  ", length(files))) else "y"
    if (resp == "y") unlink(files, force = TRUE) else NULL
  } else {
    files <- file.path(env$path, rdoi(doi))
    unlink(files, force = TRUE)
  }
}

#' @export
#' @rdname pg_data
pg_cache_list <- function() list.files(env$path)
