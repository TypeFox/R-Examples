ksLdist <-
function (n = 10, sims = 1000, alpha = 0.05)
{
    Dn <- c()
    DnL <- c()
    for (i in 1:sims) {
        x <- rnorm(n)
        mu <- mean(x)
        sig <- sd(x)
        Dn[i] <- ks.test(x, pnorm)$statistic
        DnL[i] <- ks.test(x, pnorm, mean = mu, sd = sig)$statistic
    }
    ys <- range(density(DnL)$y)
    xs <- range(density(Dn)$x)
    cv <- quantile(Dn, 1 - alpha)
    cvp <- quantile(DnL, 1 - alpha)
    plot(density(Dn, bw = 0.02), col = "blue", lwd = 2, ylim = ys, xlim = xs,
        main = "", ,xlab="", sub = paste("Simulated critical value =",
        round(cv, 3), "(simple hypothesis) and ", round(cvp, 3),
            "(composite hypothesis)\n for n =", n, "when the alpha value =",
            alpha))
    title(main = list(expression(paste("Simulated Sampling Distribution of " ,
    D[n]))))
    lines(density(DnL, bw = 0.02), col = "red", lwd = 2, lty = 2)
    legend(mean(xs), max(ys), legend = c("Simple Hypothesis",
        "Composite Hypothesis"), col = c("blue", "red"), xjust = 0,
        text.col = c("black", "black"), lty = c(1, 2), bg = "gray95",
        cex = 1, lwd = 2)
    box()
    abline(h = 0)
}

