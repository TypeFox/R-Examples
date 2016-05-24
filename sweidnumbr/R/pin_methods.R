#' @export
print.pin <- function (x,...) 
{
  print(as.character(x),...)
  cat("Personal identity number(s)")
}

#' @export
`[.pin` <- create_s3_method("[")
#' @export
rep.pin <- create_s3_method("rep")

#' @export
`[<-.pin` <- function(x, ..., value){
  value <- as.pin(value)
  NextMethod()
}

