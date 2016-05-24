Fiellers <- function(model,
                     Which = 1:length(which(b)),
                     sample = "Sample",
                     factor = paste0("factor(", sample, ")"),
                     independent = "Z",
                     df = summary(model)$df[2],
                     alpha = 0.05) {
    comp.p <- 1 - alpha / 2
    namesCoefs <- names(coef(model))
    b <- (substr(namesCoefs, 1, nchar(factor)) == factor) &
           unlist(lapply(strsplit(namesCoefs, ":"), length)) < 2
    i <- which(b)[Which]
    j <- which(namesCoefs == independent)
    alpha <- coef(model)[i]
    beta <- coef(model)[j]
    s11 <- diag(vcov(model))[i]
    s12 <- vcov(model)[i, j]
    s22 <- vcov(model)[j, j]
    m <- alpha / beta
    g <- qt(comp.p, df)^2 * s22 / beta^2
    Z <- s11 - 2 * m * s12 + m^2 * s22 - g * (s11 - s12^2 / s22)
    term1 <- g * s12 / s22
    term2 <- qt(comp.p, df) / beta * sqrt(Z)
    namesX <- c("Alpha", "Beta", "s11", "s12", "s22")
    namesF <- c("Lower", "Rho", "Upper")
    if (length(i) == 1) {
        X <- c(alpha, beta, s11, s12, s22)
        names(X) <- namesX
        Fiellers <- c(m, (m - term1 + c(-1, 1) * term2) / (1-g))[c(2, 1, 3)]
        names(Fiellers) <- namesF
        Fiellers <- c(Fiellers, X)
    } else {
        X <- cbind(alpha, beta, s11, s12, s22)
        dimnames(X)[[2]] <- namesX
        Fiellers <-
            cbind((m - term1 - term2) / (1-g), m,
                  (m - term1 + term2) / (1-g))
        dimnames(Fiellers)[[2]] <- namesF
        Fiellers <- cbind(Fiellers, X)
    }
    return(Fiellers)
}
