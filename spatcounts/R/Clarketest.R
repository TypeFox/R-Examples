`Clarketest` <-
function (LogLike1, LogLike2, alpha = 0.05, p = NULL, q = NULL, correction = TRUE){
    if (is.null(p) == TRUE) 
        p <- LogLike1$Coef
    if (is.null(q) == TRUE) 
        q <- LogLike1$Coef
    clarketab <- matrix(NA, 1, dim(LogLike1$ll)[2])
    for (i in 1:dim(LogLike1$ll)[2]) {
        ll1 <- LogLike1$ll[, i]
        ll2 <- LogLike2$ll[, i]
        n <- length(ll1)
        m.i <- ll1 - ll2 - ifelse(correction, p/2 * log(n) - 
            q/2 * log(n), 0)/n
        B <- sum(m.i > 0)
        if (B < n/2) {
            p.value <- pbinom(B, n, 0.5)
            p.value2 <- ifelse(p.value < 10^(-16), "<2e-16", 
                as.character(formatC(p.value, ifelse(p.value < 
                  10^(-4), 3, 4), format = ifelse(p.value < 10^(-4), 
                  "g", "f"))))
        }
        if (B >= n/2) {
            p.value <- 1 - pbinom(B - 1, n, 0.5)
            p.value2 <- ifelse(p.value < 10^(-16), "<2e-16", 
                as.character(formatC(p.value, ifelse(p.value < 
                  10^(-4), 3, 4), format = ifelse(p.value < 10^(-4), 
                  "g", "f"))))
        }
        rownames(clarketab) <- c("Favour model")
        clarketab[1, i] <- as.numeric(B/n)
    }
    up <- qbinom(0.975, n, 0.5)/n
    lo <- qbinom(0.025, n, 0.5)/n
    cat("Favour model 1  ", length(clarketab[clarketab > up])/dim(LogLike1$ll)[2])
    cat("\n")
    cat("No decision     ", length(clarketab[clarketab > lo & 
        clarketab < up])/dim(LogLike1$ll)[2])
    cat("\n")
    cat("Favour model 2  ", length(clarketab[clarketab < lo])/dim(LogLike1$ll)[2])
    cat("\n")
    return(clarketab)
}

