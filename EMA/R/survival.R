km <- function (time, status, group = NULL, xlab="Time (years)", ylab="", ...) {
    if (!is.null(group)) {
        group <- factor(group)
        n <- nlevels(group)
    }
    else {
        group <- factor(rep(1, length(time)))
        n <- 1
    }
    if (n == 1) {
        lr <- NULL
        p.lr <- NULL
    }
    else {
        lr <- survdiff(Surv(time, status) ~ group, na.action = na.omit)
        p.lr <- 1 - pchisq(lr$chisq, df = n - 1)
    }
    
    fit <- survfit(Surv(time, status) ~ group, na.action = na.omit)
    par(font.lab=2)
    plot(fit, col = 1:n, xlab = xlab, ylab = ylab, ...)
    if (n != 1) {
        legend("topright", inset = 0.01, legend = names(summary(group))[names(summary(group)) != 
            "NA's"], bty = "n", lwd = rep(1, n), col = 1:n, cex = 0.7, 
            text.col = "gray50")
        text(x = 0.08 * max(time), y = 0.3, paste("p=", round(p.lr, 4), sep = ""))
    }
    
    return(list(fit.km = fit, lr = lr, p.lr = p.lr))
}
