print.metaplus <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "metaplus"))
    stop("Use only with 'metaplus' objects.\n")
  print(summary(x))
  invisible()
}