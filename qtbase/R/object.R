## methods on RQtObject, which wraps every object created by Smoke

setOldClass("RQtObject")

print.RQtObject <- function(x, ...) {
  cat(head(class(x), 1), "instance\n")
}

print.RQtInvalid <- function(x, ...) {
  cat("**INVALID** reference to a", class(x)[2], "instance\n")
}
