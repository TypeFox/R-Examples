mlpm <-
function (x, p, prob, nsim, smean, svar, fre = 365, outer = 0, 
    plot = F, rain = T, over = T, estimate = "ols", CCFlag = 20, 
    Plag = 20, lsign = 0.05, des = T) 
{
    options(object.size = 1e+08, warn = -1)
    x <- data.matrix(x)
    n <- nrow(x)
    k <- ncol(x)
    y <- matrix(0, nrow(x), ncol(x))
    x1 <- x
    if (des == T) {
        cat("Deseasonalization", "\n")
        x1 <- multides(x, smean, svar, fre = fre, outer)
    }
    for (h in 1:k) y[, h] <- (x1[, h] - mean(x1[, h]))
    cat("Estimation", "\n")
    if (estimate == "ols") {
        y1 <- ar.ols(x1, aic = F, order.max = p, demean = T, 
            intercept = F)
        if (p == 1) 
            rr <- y1$asy.se.coef$ar
        else rr <- varcoeff(k, p, y1$asy.se.coef)
        if (over) {
            cat("Restrictions", "\n")
            rst <- restrict(varcoeff(k, p, y1), rr, n, prob)
            c <- count(sync(rst))
            cat("Parameters estimated", c, "out of", p * k^2, 
                "\n")
            y1 <- ar.egls(x1, rst, order.max = p)
            if (p == 1) 
                rr <- y1$asy.se.coef$ar
            else rr <- varcoeff(k, p, y1$asy.se.coef)
        }
    }
    if (estimate == "burg") 
      if ((k == 1)) y1 <- ar.burg(as.vector(x1), aic = F, order.max = p)
    else y1 <- ar(ts(x1,class="mts"), aic = F, order.max = p, method="burg")
    
    if (estimate == "yw") 
        y1 <- ar.yw(x1, aic = F, order.max = p)
    if (k == 1) 
        res1 <- matrix(0, length(na.exclude(y1$resid)), 1)
    else res1 <- matrix(0, nrow(as.matrix(na.exclude(y1$resid))), 
        ncol(na.exclude(y1$resid)))
    for (q in 1:ncol(as.matrix(na.exclude(y1$resid)))) res1[, 
        q] <- as.matrix(na.exclude(y1$resid))[, q]
    sim <- 0
    if (nsim > 0) {
        cat("Simulation", "\n")
        if ((p == 1)) 
            sim <- simvar(varcoeff(k, p, y1), acf(x1, type = "covariance", 
                plot = F)$acf[1, , ], res1, n, k, p, nsim, 1)
        else sim <- simvar(varcoeff(k, p, y1), sigmaZ(x1, p), 
            res1, n, k, p, nsim, 1)
        for (r in 1:nsim) {
            for (s in 1:k) sim[[r]][, s] <- (sim[[r]][, s] + 
                mean(x1[, s]))
            if (des == T) {
                cat("Seasonalization sim", r, "\n")
                sim[[r]] <- multistag(sim[[r]], x, smean, svar, 
                  fre = fre, outer)
            }
            if (rain) {
                cat("Rain adaptor sim", r, "\n")
                sim[[r]] <- multipro(x, sim[[r]])
            }
        }
    }
    nr <- nrow(x1)
    if (k == 1) {
        sk = var(res1)
        som = var(y)
    }
    else {
        skm <- determinant(acf(res1, type = "covariance", plot = F)$acf[1, 
            , ], logarithm = F)
        sk <- skm$sign * skm$modulus[1]
        som <- determinant(acf(y, type = "covariance", plot = F)$acf[1, 
            , ], logarithm = F)
        so <- som$sign * som$modulus[1]
    }
    np <- p * k^2
    if ((estimate == "ols") && (over)) 
        np <- c
    Qst <- Portmanteau(res1, Plag)
    cat("Portmanteau Test :", Qst, qchisq((1 - lsign), ((k^2) * 
        Plag - np)), "\n")
    aic <- nr * log(sk) + 2 * np
    sbc <- nr * log(sk) + (np * log(nr))
    cat("AIC :", aic, "\n")
    cat("SBC :", sbc, "\n")
    result <- list()
    result$coeff <- varcoeff(k, p, y1)
    if (estimate == "ols") 
        result$coeffstd <- rr
    if ((estimate == "ols") && (over)) 
        result$struct <- rst
    result$res <- res1
    if (plot) {
        plot(acf(x, lag.max = CCFlag, type = "partial"), ylab = "PACF Series", 
            max.mfrow = 2)
        dev.new()
        plot(acf(x, lag.max = CCFlag), ylab = "CCF Series", max.mfrow = 2)
        dev.new()
        plot(acf(res1, lag.max = CCFlag), ylab = "CCF Residuals", 
            max.mfrow = 2)
        if (nsim > 0) {
            dev.new()
            plot(acf(sim[[1]], lag.max = CCFlag), ylab = "CCF Simulated Series", 
                max.mfrow = 2)
            dev.new()
            plot(acf(sim[[1]], lag.max = CCFlag, type = "partial"), 
                ylab = "PACF Simulated Series", max.mfrow = 2)
        }
    }
    result$PACRx <- acf(x, lag.max = CCFlag, , plot = F)
    result$CCRx <- acf(x, lag.max = CCFlag, plot = F)
    result$CCRres <- acf(res1, lag.max = CCFlag, plot = F)
    if (nsim > 0) {
        result$CCRsim1 <- acf(sim[[1]], lag.max = CCFlag, plot = F)
        result$PACRsim1 <- acf(sim[[1]], lag.max = CCFlag, type = "partial", 
            plot = F)
    }
    result$fit <- y1
    result$aic <- aic
    result$Qst <- Qst
    result$sim <- sim
    return(result)
    options(warn = -1)
}

