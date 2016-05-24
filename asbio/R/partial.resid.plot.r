partial.resid.plot <- function (x, smooth.span = 0.8, lf.col = 2, sm.col = 4,...) 
{
    p <- ncol(x$model)
    X <- as.matrix(x$model[,2:p])
    Y <- x$model[,1]
    resid <- matrix(ncol = ncol(X), nrow = nrow(X))
    other.resid <- matrix(ncol = ncol(X), nrow = nrow(X))
    for (i in 1:ncol(X)) {
        resid[, i] <- lm(Y ~ X[, -i])$residuals
        other.resid[, i] <- lm(X[, i] ~ X[, -i])$residuals
    }
    for (i in 1:ncol(X)) {
        plot(other.resid[, i], resid[, i], xlab = bquote(paste(hat(epsilon), 
            "(", .(colnames(X)[i]), " | model without ", .(colnames(X)[i]), 
            ")")), ylab = bquote(paste(hat(epsilon), "(Y | model without ", 
            .(colnames(X)[i]), ")")),...)
        l1 <- lm(resid[, i] ~ other.resid[, i])
        abline(as.numeric(l1$coefficients[1]), as.numeric(l1$coefficients[2]), 
            col = lf.col, lty = 1)
        lines(lowess(other.resid[, i], resid[, i], f = smooth.span), 
            lty = 2, col = sm.col)
        readline("Press return for next plot")
    }
}