up.sample <- function(x, f, y=NA) {
  n <- length(x)
  as.vector(rbind(x, matrix(rep(y, (f-1)*n), nrow=f-1)))
}
