varcheck <- function(x, vars, msg) {

  check <- pmatch(vars, names(x))

  if(any(is.na(check))) {
    if(missing(msg))
      stop("Missing or misnamed variables:",
        paste(vars[is.na(check)], collapse = ", "))
    else
      stop(msg)
  }
}
