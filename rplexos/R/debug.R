# Custom function to print debug messages
rplexos_message <- function(...) {
  if (is_debug_rplexos()) {
    message("*** rplexos debug: ", ...)
  }
}

#' Enable or disable debug mode
#'
#' The debug mode will print progress on screen and additional information to 
#' help diagnose problems.
#'
#' @export
start_debug_rplexos <- function() {
  options(rplexos.debug = TRUE)
  check_debug_rplexos()
}

#' @rdname start_debug_rplexos
#' @export
stop_debug_rplexos <- function() {
  options(rplexos.debug = FALSE)
  check_debug_rplexos()
}

#' @rdname start_debug_rplexos
#' @export
check_debug_rplexos <- function() {
  out <- is_debug_rplexos()
  
  if (out) {
    cat("rplexos debug mode is enabled\n")
  } else {
    cat("rplexos debug mode is disabled\n")
  }
  
  invisible(out)
}

#' @rdname start_debug_rplexos
#' @export
is_debug_rplexos <- function() {
  getOption("rplexos.debug", FALSE)
}

#' Shortcut functions for rplexos sample files
#'
#' These functions return the folder where the PLEXOS sample files
#' are located. They are uses in different examples and the vignettes.
#'
#' @export
location_solution_rplexos <- function() {
  system.file("extdata", "solution", package = "rplexos")
}

#' @rdname location_solution_rplexos
#' @importFrom utils unzip
#' @export
location_input_rplexos <- function() {
  out <- system.file("extdata", "database", package = "rplexos")
  
  # File comes compressed by default. Decompress it
  xml.path <- file.path(out, "three_nodes.xml")
  zip.path <- system.file("extdata", "three_nodes.zip", package = "rplexos")
  if (!file.exists(xml.path))
    unzip(zip.path, exdir = out)
  
  # Return folder location
  out
}
