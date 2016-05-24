#' @export 
rev.ff <- function(x){
  y <- clone(x)
  N <- length(x) + 1L
  for (i in chunk(x)){
    y[i] <- x[N-i]
  }
  y
}

# #' quick testing
# x <- ff(1:100)
# rev(x)
