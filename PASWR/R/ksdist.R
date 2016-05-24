ksdist <-
function (n = 10, sims = 10000, alpha = 0.05)
{
    Dn <- replicate(sims, ks.test(rnorm(n), pnorm)$statistic)
    cv <- quantile(Dn, 1 - alpha)
    plot(density(Dn), col = "blue", lwd = 2, main = "",
    xlab = paste("Simulated critical value =", round(cv,3),
    "for n =", n, "when the alpha value =", alpha))
    title(main = list(expression(paste("Simulated Sampling Distribution of " ,
    D[n]))))
}

