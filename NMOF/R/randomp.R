## -*- truncate-lines: t; -*-
## create random portfolios

## n .. number of assets available
## k .. number of assets in solution
## sum_w .. sum of weights (typically 1)
## abssum_w .. sum of absolute weights (useful for long/short portfolios)
## wmin ..
## wmax ..

long.only <- function(n, budget = 1) {
    ans <- runif(n)
    budget * ans/sum(ans)
}

long.only.k <- function(n, k, budget = 1) {
    ans <- numeric(n)
    i <- sample.int(n, k)
    ans[i] <- runif(length(i))
    budget*ans/sum(ans)
}

long.short <- function(n, budget = 1) {
    ans <- runif(n)
    while (all(ans >= 0) || all(ans < 0)) {
        tmp <- runif(n) < 0.3
        ans[tmp] <- -ans[tmp]
    }
    budget*ans/sum(ans)
}
