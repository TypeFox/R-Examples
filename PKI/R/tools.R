raw2hex <- function(what, sep, upper=FALSE)
  .Call(PKI_raw2hex, what, if (missing(sep)) NULL else sep, upper)
