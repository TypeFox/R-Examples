o <- c(315, 102, 108, 31)
n <- sum(o)
e <- n * c(9, 3, 3, 1) / 16; e
G <- 2 * sum( o * log(o/e) ); G
pval <- 1-pchisq(G, 3); pval
