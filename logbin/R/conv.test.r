conv.test <- function(theta1, theta2, epsilon) {
  diffvec <- theta1 - theta2
  diff <- sqrt(sum(diffvec**2))/sqrt(sum(theta1**2))
  (diff < epsilon)
}