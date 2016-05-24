alternate.order <-
function(input.vector,order.vector) {
  vector.length <- length(input.vector)
  if (vector.length != length(order.vector) || vector.length != max(order.vector)) {
    return(NULL)
  }
  output.vector <- vector(length=vector.length)
  assign.value <- function(x) {
    output.vector[order.vector[x]] <<- input.vector[x]
  }
  sapply(seq(vector.length),assign.value)
  return(output.vector)
}
