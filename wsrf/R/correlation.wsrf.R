correlation <- function(object, ...) UseMethod("correlation")

correlation.wsrf <- function(object, ...) {

  object[[.CORRELATION_IDX]]

}