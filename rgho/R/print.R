to_data_frame <- function(x, n = length(x)) {
  data.frame(
    Label = head(attr(x, "labels"), n = n),
    ID = head(as.vector(x), n = n)
  )
}

#' @export
print.gho <- function(x, n = 6, ...) {
  if (n == Inf) {
    n <- length(x)
  }

  cat(sprintf("A 'GHO' object of %i elements.\n\n", length(x)))

  print(to_data_frame(x, n = n), ...)

  cat(sprintf("...\n\n(Printing %i first elements.)\n", n))

  if (! is.null(attr(x, "attrs"))) {
    names_attr <- names(attr(x, "attrs"))
    names_attr <- names_attr[names_attr != "code"]
    cat("\nAttributes:\n\n")
    cat(names_attr, sep = "\n")
  }
}
