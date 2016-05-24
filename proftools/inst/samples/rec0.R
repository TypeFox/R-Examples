f <- function(n) { n <- n - 1; g(n) }
g <- function(n) { if (n > 0) f(n) else 0 }
