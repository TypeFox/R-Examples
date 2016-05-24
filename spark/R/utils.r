
#' @importFrom magrittr is_less_than

which_max <- function(x, tolerance = .Machine$double.eps^0.5) {
  abs(max(x) - x) %>%
    is_less_than(tolerance) %>%
    which()
}

#' @importFrom utils installed.packages

is_installed <- function(pkg) {
  pkg %in% rownames(installed.packages())
}


import <- function(pkg, name) {
  getExportedValue(pkg, name)
}
