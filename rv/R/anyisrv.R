

anyisrv <- function (...) {
  any(unlist(lapply(list(...), is.rv)))
}




