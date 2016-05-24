
isa.status.default <- function(...) {
  if (isa.option("verbose")) {
    args <- list(...)
    cat(args[[1]], "\n")
  }
}

isa.status <- function(...) {
  func <- isa.option("status.function")
  func(...)
}
