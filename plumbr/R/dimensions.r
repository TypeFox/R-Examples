#' @S3method nrow mutaframe
nrow.mutaframe <- function(x) length(attr(x, "row.names")) %||% 0
#' @S3method ncol mutaframe
ncol.mutaframe <- function(x) length(attr(x, "col.names")) %||% 0
#' @S3method dim mutaframe
dim.mutaframe <- function(x) c(nrow.mutaframe(x), ncol.mutaframe(x))
