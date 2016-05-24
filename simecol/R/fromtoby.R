## the fromtoby-vector form is intended to give users
## a compact representation of time steps without the need to apply
## time series classes as "ts" or "zoo"
## "fromtoby" vectors are optional, ordinary sequences are possible too

## check for from-to-by structure of vector
isfromtoby <- function(x) {
  (sum(names(x) %in% c("from","to","by"))==3)
}

hasfromtoby <- function(x) {
  any(names(x) %in% c("from","to","by"))
}

## expand from-to-by-vector to sequence
fromtoby <- function(times) {
  if (isfromtoby(times)) {
    times <- seq(times["from"], times["to"], times["by"])
  } 
  times
}
