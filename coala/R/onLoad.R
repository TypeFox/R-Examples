.onLoad <- function(libname, pkgname) {
  # scrm should always be available
  activate_scrm()

  # Silently create the simulators for which a binary is available
  suppressMessages({
    tryCatch(activate_ms(), error = function(e) {}) #nolint
    tryCatch(activate_msms(), error = function(e) {}) #nolint
    tryCatch(activate_seqgen(), error = function(e) {}) #nolint
  })

  invisible(NULL)
}
