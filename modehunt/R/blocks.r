blocks <- function(n, m0 = 10, fm = 2){

# generate range of number of observations for each block
n_blocks <- floor(log(n / m0) / log(fm))

l1 <- floor(m0 * fm ^ (0:n_blocks))
l2 <- floor(m0 * fm ^ (1:(n_blocks + 1)) - 1)
res <- cbind(l1, l2)
res <- res[length(res[, 1]):1, ]

return(matrix(res, ncol = 2))
}
