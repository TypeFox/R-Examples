#' Loads available packages from the given repository.
#' @param repos character vector, the base URLs of the repositories to use
#' @return data.frame of available packages
load_available_packages <- function(repos) {
  con <- url(file.path(contrib.url(repos), 'PACKAGES'))
  available <- as.data.frame(read.dcf(con))
  close(con)
  available
}

