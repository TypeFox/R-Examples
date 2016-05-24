
#' @S3method names mutaframe
names.mutaframe <- function(x, ...) attr(x, "col.names")

#' @S3method dimnames mutaframe
dimnames.mutaframe <- function(x, ...) {
  list(attr(x, "row.names"), attr(x, "col.names"))
}

#' @S3method "names<-" mutaframe
`names<-.mutaframe` <- function(x, ..., value) {
  attr(x, "col.names") <- value
  x
}

#' @S3method "dimnames<-" mutaframe
`dimnames<-.mutaframe` <- function(x, ..., value) {
  attr(x, "row.names") <- value[[1]]
  attr(x, "col.names") <- value[[2]]
  x
}

#' Make valid variable names
#' @param var_names variable names
variable_names <- function(var_names) {
  no_name <- is.na(var_names) | nzchar(var_names) == 0
  var_names[no_name] <- paste("X", seq_len(sum(no_name)), sep = "")

  make.names(var_names, unique = T)
}
