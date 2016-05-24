roc.nsc <-function(object) {
###Computes the roc curve for a nsc model
  nonzero <- object$nonzero^(1/4)
  errors <- object$errors
  if(is.null(errors))
    stop("No errors component")
  n <- length(errors)
  heights <- (errors[1:(n - 1)] + errors[2:n])/2
  bases <- diff(nonzero)
  area <- sum((nonzero[-1] + nonzero[-n]) * heights * bases) /
    (-2 * diff(range(nonzero)))
  area
}
