flowsRandom <-
function(k=300, n=150) {
  # Input:
  # k = number of flows
  # n = number of base pairs
  # Output:
  # P = percentage of sequences that are completely covered
  # Q = percentage of sequences that are not completely covered
  # delta = 1 - (P + Q), should be a small number when calculation is accurate
  # N = number of summands for calculation of P
  # M = number of summands for calculation of Q
  # summands = N + M = total number of summands
  prob <- c(1/4, 1/4, 1/4, 1/4)
  N <- 0
  M <- 0
  P <- 0
  Q <- 0
  for (n1 in 0:n) {
    cat(".")
    for (n2 in 0:(n-n1)) {
      for (n3 in 0:(n-n1-n2)) {
        n0 <- n - n1 - n2 - n3
        l <- n1 + 2*n2 + 3*n3 + 1
        if (l <= k) {
          N <- N + 1
          p <- dmultinom(c(n0, n1, n2, n3), prob=prob)
          P <- P + p
        }
        else {
          M <- M + 1
          q <- dmultinom(c(n0, n1, n2, n3), prob=prob)
          Q <- Q + q
        }
      }
    }
  }
  cat("\n")
  result <- vector()
  result["P"] <- P
  result["Q"] <- Q
  result["delta"] <- 1 - (P +  Q)
  result["N"] <- N
  result["M"] <- M
  result["summands"] <- N + M
  return(result)      
}
