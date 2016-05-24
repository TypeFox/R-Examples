`Vuongtest` <-
function (LogLike1, LogLike2, alpha = 0.05, p = NULL, q = NULL, correction = TRUE){
    if (is.null(p) == TRUE) 
        p <- LogLike1$Coef
    if (is.null(q) == TRUE) 
        q <- LogLike1$Coef
    vuongtab <- matrix(NA, 3, dim(LogLike1$ll)[2])
    for (i in 1:dim(LogLike1$ll)[2]) {
        ll1 <- LogLike1$ll[, i]
        ll2 <- LogLike2$ll[, i]
        n <- length(ll1)
        m.i <- ll1 - ll2
        nu <- (sqrt(n) * 1/n * (sum(m.i) - ifelse(correction, 
            p/2 * log(n) - q/2 * log(n), 0)))/(sqrt((n - 1)/n * 
            var(m.i)))
        if (abs(nu) < qnorm(1 - alpha/2)) {
            fav <- 0
        }
        if (nu >= qnorm(1 - alpha/2)) {
            fav <- 1
        }
        if (nu <= -qnorm(1 - alpha/2)) {
            fav <- 2
        }
        rownames(vuongtab) <- c("nu", "Favour model", 
            "P-value")
        vuongtab[1, i] <- as.numeric(formatC(nu, 3, format = "g"))
        vuongtab[2, i] <- as.integer(fav)
        vuongtab[3, i] <- as.numeric(2 * pnorm(-abs(nu)))
        rm(ll1, ll2)
    }
    cat("Favour model 1  ", length(vuongtab[2, ][vuongtab[2, 
        ] == 1])/dim(LogLike1$ll)[2])
    cat("\n")
    cat("No decision     ", length(vuongtab[2, ][vuongtab[2, 
        ] == 0])/dim(LogLike1$ll)[2])
    cat("\n")
    cat("Favour model 2  ", length(vuongtab[2, ][vuongtab[2, 
        ] == 2])/dim(LogLike1$ll)[2])
    cat("\n")
    return(vuongtab)
}

