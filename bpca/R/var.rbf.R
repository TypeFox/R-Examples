var.rbf <- function(x)
{
  lv <- function(x) sqrt(t(x) %*% x)  # length of vector
  l  <- apply(x, 1, lv)
  n  <- nrow(x)
  var.rb <- diag(n)

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {

      cost <- (t(x[i,]) %*%
               x[j,]) /
      (l[i]*l[j])

      var.rb[j,i] <- cost         # fill lower.tri
      var.rb[i,j] <- var.rb[j,i]  # fill upper.tri
    }
  }

  dimnames(var.rb) <- list(dimnames(x)[[1]],
                           dimnames(x)[[1]])
  return(var.rb)
}
