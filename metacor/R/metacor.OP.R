`metacor.OP` <-
function(r, n, labels, alpha = 0.05, plot = TRUE, xlim = c(-1, 1)){
G <- r * (1 + ((1 - r^2) / (2 * (n - 1 - 3))))
names(G) <- labels
G.var <- G^2 - 1 + (((n - 3) * (1 - r^2) * hyperg_2F1(1, 1, n / 2, 1 - r^2)) / (n - 2))
names(G.var) <- labels
G.lower <- G - qnorm(alpha / 2) * sqrt(G.var)
G.upper <- G + qnorm(alpha / 2) * sqrt(G.var)
N <- sum(n)
k <- length(r)
G.mean <- sum(n * G) / N
G.mean.se <- (1 - G.mean^2) / sqrt(N - k)
g <- G.mean / G.mean.se
G.mean.lower <- G.mean - qnorm(alpha / 2) * G.mean.se
G.mean.upper <- G.mean + qnorm(alpha / 2) * G.mean.se
p <- pnorm(abs(g), lower.tail = F)
if (plot) metaplot(G, sqrt(G.var), labels = labels, xlab = quote("Unique minimum variance unbiased estimator"~italic(G)), ylab = "", summn = G.mean, sumse = G.mean.se, sumnn = G.mean.se^-2, xlim = c(-1, 1))

res <- list()
res$G <- G
res$G.var <- G.var
res$G.lower <- G.lower
res$G.upper <- G.upper
res$G.mean <- G.mean
res$G.mean.se <- G.mean.se
res$G.mean.lower <- G.mean.lower
res$G.mean.upper <- G.mean.upper
res$p <- p
invisible(res)



}

