`DIC` <-
function (Yin, fm.X, region, regmodel, burnin = 1){
    phi <- regmodel$phi
    omega <- regmodel$omega
    r <- regmodel$r
    be <- regmodel$beta
    ga <- regmodel$gamma
    t.i <- regmodel$t.i
    if (fm.X == ~1) 
        xb <- matrix(1, length(Yin), 1)
    else xb <- model.matrix(fm.X)
    y <- as.matrix(Yin)
    fm.X <- cbind(region, xb)
    regindex <- fm.X[, 1]
    if (is.null(r) & is.null(omega) == TRUE) {
        n <- length(phi)
        result <- matrix(0, 1, (n - burnin + 1))
        for (i in burnin:n) {
            mu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% be[, i] + 
                ga[regindex, i])
            result[i - burnin + 1] <- -2 * sum(log(mu.i) + (y - 
                1) * log(mu.i + (phi[i] - 1) * y) - y * log(phi[i]) - 
                1/phi[i] * (mu.i + (phi[i] - 1) * y) - lgamma(y + 
                1))
        }
        meandeviance <- mean(result)
        pphi <- mean(phi[burnin:n])
        pbe <- matrix(0, 1, dim(be)[1])
        pga <- matrix(0, 1, dim(ga)[1])
        for (i in 1:dim(be)[1]) {
            pbe[i] <- sum(be[i, burnin:n])/(n - burnin + 1)
        }
        for (i in 1:dim(ga)[1]) {
            pga[i] <- sum(ga[i, burnin:n])/(n - burnin + 1)
        }
        mu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% t(pbe) + pga[regindex])
        devmean <- -2 * sum(log(mu.i) + (y - 1) * log(mu.i + 
            (pphi - 1) * y) - y * log(pphi) - 1/pphi * (mu.i + 
            (pphi - 1) * y) - lgamma(y + 1))
        effpar <- meandeviance - devmean
        dicresult <- meandeviance + effpar
        cat("DIC            ", dicresult)
        cat("\n")
        cat("mean deviance  ", meandeviance)
        cat("\n")
        cat("p.D            ", effpar)
        cat("\n")
    }
    if (is.null(r) == FALSE) {
        n <- length(r)
        result <- matrix(0, 1, (n - burnin + 1))
        for (i in burnin:n) {
            mu.i <- t.i * exp(xb %*% be[, i] + ga[regindex, i])
            result[i - burnin + 1] <- -2 * sum(lgamma(r[i] + 
                y) - lgamma(r[i]) - lgamma(y + 1) + r[i] * log(r[i]) + 
                y * log(mu.i + (y == 0)) - (r[i] + y) * log(r[i] + 
                mu.i))
        }
        meandeviance <- mean(result)
        pr <- mean(r[burnin:n])
        pbe <- matrix(0, 1, dim(be)[1])
        pga <- matrix(0, 1, dim(ga)[1])
        for (i in 1:dim(be)[1]) {
            pbe[i] <- sum(be[i, burnin:n])/(n - burnin + 1)
        }
        for (i in 1:dim(ga)[1]) {
            pga[i] <- sum(ga[i, burnin:n])/(n - burnin + 1)
        }
        mu.i <- t.i * exp(xb %*% t(pbe) + pga[regindex])
        devmean <- -2 * sum(lgamma(pr + y) - lgamma(pr) - lgamma(y + 
            1) + pr * log(pr) + y * log(mu.i + (y == 0)) - (pr + 
            y) * log(pr + mu.i))
        effpar <- meandeviance - devmean
        dicresult <- meandeviance + effpar
        cat("DIC            ", dicresult)
        cat("\n")
        cat("mean deviance  ", meandeviance)
        cat("\n")
        cat("p.D            ", effpar)
        cat("\n")
    }
    if (is.null(omega) == FALSE) {
        n <- length(phi)
        result <- matrix(0, 1, (n - burnin + 1))
        for (i in burnin:n) {
            mu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% be[, i] + 
                ga[regindex, i])
            result[i - burnin + 1] <- -2 * sum(ifelse(y == 0, 
                1, 0) * (log(omega[i] + (1 - omega[i]) * exp(-1/phi[i] * 
                mu.i))) + ifelse(y > 0, 1, 0) * (log(1 - omega[i]) + 
                log(mu.i) + (y - 1) * log(mu.i + (phi[i] - 1) * 
                y) - y * log(phi[i]) - 1/phi[i] * (mu.i + (phi[i] - 
                1) * y) - lgamma(y + 1)))
        }
        meandeviance <- mean(result)
        pphi <- mean(phi[burnin:n])
        pomega <- mean(omega[burnin:n])
        pbe <- matrix(0, 1, dim(be)[1])
        pga <- matrix(0, 1, dim(ga)[1])
        for (i in 1:dim(be)[1]) {
            pbe[i] <- sum(be[i, burnin:n])/(n - burnin + 1)
        }
        for (i in 1:dim(ga)[1]) {
            pga[i] <- sum(ga[i, burnin:n])/(n - burnin + 1)
        }
        pmu.i <- t.i * exp(xb[, 1:dim(be)[1]] %*% t(pbe) + pga[regindex])
        devmean <- -2 * sum(ifelse(y == 0, 1, 0) * (log(pomega + 
            (1 - pomega) * exp(-1/pphi * pmu.i))) + ifelse(y > 
            0, 1, 0) * (log(1 - pomega) + log(pmu.i) + (y - 1) * 
            log(pmu.i + (pphi - 1) * y) - y * log(pphi) - 1/pphi * 
            (pmu.i + (pphi - 1) * y) - lgamma(y + 1)))
        effpar <- meandeviance - devmean
        dicresult <- meandeviance + effpar
        cat("DIC            ", dicresult)
        cat("\n")
        cat("mean deviance  ", meandeviance)
        cat("\n")
        cat("p.D            ", effpar)
        cat("\n")
    }
    DIC <- list(DIC = dicresult, mean.defiance = meandeviance, 
        p.D = effpar)
    return(DIC)
}

