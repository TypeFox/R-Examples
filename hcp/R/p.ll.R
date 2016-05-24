p.ll <-
function(n, j, k, s2, t2){
 q1 <- n * log(sqrt(2 * pi))
 q2 <- 0.5 * (n - k + j) * (1 + log(s2))
 q3 <- 0.5 * (k - j) * (1 + log(t2))
 - (q1 + q2 + q3)
}
