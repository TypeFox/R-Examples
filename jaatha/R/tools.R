# --- Argument checking --------------------------------------------------------
is_single_numeric <- function(value) is.numeric(value) && length(value) == 1
is_single_logical <- function(value) is.logical(value) && length(value) == 1


# --- Seeds --------------------------------------------------------------------
sample_seed <- function(n = 1) sample.int(2 ^ 20, n)


# --- Check for suggested packages
require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install package '", pkg, "'", call. = FALSE)
  }
  invisible(TRUE)
} 
