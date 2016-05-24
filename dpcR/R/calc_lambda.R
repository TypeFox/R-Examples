#calculate values of lambda (average number of molecules per partition)
calc_lambda <- function(k, n) {
  #   if (length(k) > 1) {
  #     k_n <- get_k_n(k)
  #     k <- k_n[1]
  #     n <- k_n[2]
  #   }
  p <- k/n
  p_conf <- sqrt(p * (1 - p)/n) #dube
  u_lambda <- sqrt(p/(n * (1 - p))) #bhat
  l <- fl(p)
  lower <- c(fl(p - qnorm(0.975)*p_conf), l - u_lambda)
  upper <- c(fl(p + qnorm(0.975)*p_conf), l + u_lambda)
  res <- data.frame(method = c(rep("dube", length(p)), rep("bhat", length(p))), 
                    lambda = c(l,l), lambda.low = lower, lambda.up = upper, m = l*n,
                    m.low = lower*n, m.up = upper*n, 
                    k = rep(k, 2), n = rep(n, 2), row.names = 1L:(length(l)*2))
  #reorder
  res[unlist(lapply(1L:(nrow(res)/2), function(i) c(i, i + nrow(res)/2))), ]
}
