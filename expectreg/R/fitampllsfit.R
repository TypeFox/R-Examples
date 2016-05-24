fitampllsfit <-
function (y, B, b, p, cc, DD, lambda, nb) 
{
    m <- nrow(B)
    n <- ncol(B)
    np <- length(p)
    a <- B %*% b
    BtB <- 0
    Bty <- 0
    w = matrix(NA, nrow = length(y), ncol = np)
    for (j in 1:np) {
        w[, j] = p[j]
        w[!(y >= cc[j] * a), j] = 1 - p[j]
        WB <- cc[j] * as.vector(w[, j]) * B
        BtB <- BtB + cc[j] * t(WB) %*% B
        Bty <- Bty + t(WB) %*% y
    }
    lambda = c(rep(0, times = n - sum(nb)), rep(lambda, times = nb))
    P <- sqrt(lambda) * t(DD) %*% DD
    model <- lsfit(x = BtB + P, y = Bty, intercept = FALSE)
    sigma = 0.5 * sum((Bty - (BtB + P) %*% model$coef)^2, na.rm = TRUE)/(length(Bty) - 
        sum(hat(model$qr)[1:length(Bty)]))
    return(list(b = model$coef, hatma = hat(model$qr)[1:length(Bty)], 
        weight = w, sig = sigma))
}
