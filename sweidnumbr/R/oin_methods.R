
#' @export
print.oin <- function (x,...) 
{
  print(as.character(x),...)
  cat("Organizational identity number(s)")
}

#' @export
`[.oin` <- create_s3_method("[")

#' @export
rep.oin <- create_s3_method("rep")

#' @export
`[<-.oin` <- function(x, ..., value){
  value <- as.oin(value)
  NextMethod()
}
