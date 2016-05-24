printObject = function(x, type) {
  info = x$info()
  cat(
    sprintf("%s on path %s", type, info$path),
    sprintf("  %-9s : %s", "extension", info$extension),
    sprintf("  %-9s : %s", "use.cache", info$use.cache),
    sprintf("  %-9s : %s", "simplify", info$simplify),
    sprintf("  %-9s : %i", "items", length(x$ls())),
    sprintf("  %-9s : %i", "cached", length(x$cached())),
    sprintf("  %-9s : %s", "functions", collapse(names(x))),
    sprintf("  %-9s : %s", "methods", collapse(sub("\\.fail$", "", methods(class = "fail")), ", ")),
    sep = "\n"
  )
}

#' @export
print.fail = function(x, ...) {
  printObject(x, "File Abstraction Interface Layer")
}

#' @export
print.sail = function(x, ...) {
  printObject(x, "Source Abstraction Interface Layer")
}

#' @export
as.list.fail = function(x, ...) {
  x$as.list(...)
}

#' @export
as.list.sail = function(x, ...) {
  x$as.list(...)
}
