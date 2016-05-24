qap <- function(A, B, method = NULL, ...) {
  if(is.null(method)) method <- "SA"

  methods <- c("SA")
  method <- methods[pmatch(tolower(method), tolower(methods))]
  if(is.na(method)) stop("Unknown method. Available methods are: ",
    paste(methods, collapse = ", "))

  if(method == "SA") qapSA(A, B, ...)
  else stop("Unknown method. Available methods are: ",
    paste(methods, collapse = ", "))
}

qap.obj <- function(A, B, o) {
  sum(diag(A%*%B[o,o]))
}
