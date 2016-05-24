tempfile <- function(name="unnamed") {
  base::tempfile(paste0("coala-", Sys.getpid(), "-", name, "-"))
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Please install package '", pkg, "'", call. = FALSE)
  }
  invisible(TRUE)
}
