"iARMAapprox" <- function(phi = numeric(0), theta = numeric(0), seas = F) {
    maxlag <- 2048
    if (!(InvertibleQ(phi) & InvertibleQ(theta))) {
        return(NULL)
    }
    p <- length(phi)
    q <- length(theta)
    unames <- character(0)
    vnames <- character(0)
    if (p > 0) {
        psisphi <- psiwtsAR(phi, maxlag = maxlag)
        vv <- numeric(p)
        vv[1] <- crossprod(psisphi)
        if (p > 1) 
            vv[2:p] <- sapply(1:(p - 1), function(i) crossprod(psisphi[-(1:i)], psisphi[-((maxlag - 
                i + 1):maxlag)]))
        vvmatrix <- toeplitz(vv)
        imatrix <- vvmatrix
        if (!seas) 
            vnames <- paste("phi(", 1:p, ")", sep = "") else vnames <- paste("phiseas(", 1:p, ")", sep = "")
    }
    if (q > 0) {
        psisth <- -psiwtsAR(theta, maxlag = maxlag)
        uu <- numeric(q)
        uu[1] <- crossprod(psisth)
        if (q > 1) 
            uu[2:q] <- sapply(1:(q - 1), function(i) crossprod(psisth[-(1:i)], psisth[-((maxlag - 
                i + 1):maxlag)]))
        uumatrix <- toeplitz(uu)
        imatrix <- uumatrix
        if (!seas) 
            unames <- paste("theta(", 1:q, ")", sep = "") else unames <- paste("thetaseas(", 1:q, ")", sep = "")
    }
    if (p > 0 && q > 0) {
        uvmatrix <- matrix(numeric(1), nrow = p, ncol = q)
        uvmatrix[1, 1] <- crossprod(psisphi, psisth)
        if (q > 1) 
            uvmatrix[1, -1] <- sapply(1:(q - 1), function(i) crossprod(psisphi[-(1:i)], 
                psisth[-((maxlag - i + 1):maxlag)]))
        if (p > 1) 
            uvmatrix[-1, 1] <- sapply(1:(p - 1), function(i) crossprod(psisth[-(1:i)], psisphi[-((maxlag - 
                i + 1):maxlag)]))
        if (p > 1 && q > 1) {
            for (i in 2:p) {
                for (j in 2:q) {
                  if (i == j) 
                    uvmatrix[i, j] <- uvmatrix[1, 1] else if ((k = i - j) < 0) 
                    uvmatrix[i, j] <- crossprod(psisphi[-(1:(-k))], psisth[-((maxlag + k + 
                      1):maxlag)]) else uvmatrix[i, j] <- crossprod(psisth[-(1:k)], psisphi[-((maxlag - k + 
                    1):maxlag)])
                }
            }
        }
        imatrix <- cbind(rbind(vvmatrix, t(uvmatrix)), rbind(uvmatrix, uumatrix))
    }
    inames <- c(vnames, unames)
    dimnames(imatrix) <- list(inames, inames)
    return(imatrix)
}
"iARMASeasMultapprox" <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), period = 0) {
    if (period == 0) {
        warning("Seasonal call, but period = 0; reverting to nonseasonal model")
        return(iARMAapprox(phi = phi, theta = theta))
    }
    if (!(InvertibleQ(phi) & InvertibleQ(theta) & InvertibleQ(phiseas) & InvertibleQ(thetaseas))) {
        return(NULL)
    }
    maxlag <- 2048
    maxlag1 <- maxlag * period
    imatrix <- iarma <- isarma <- block <- ppsmatrix <- qpsmatrix <- pqsmatrix <- qqsmatrix <- NULL
    p <- length(phi)
    q <- length(theta)
    ps <- length(phiseas)
    qs <- length(thetaseas)
    
    if (p + q > 0) 
        iarma <- iARMAapprox(phi = phi, theta = theta)
    if (ps + qs > 0) 
        isarma <- iARMAapprox(phi = phiseas, theta = thetaseas, seas = T)
    inames <- c(dimnames(iarma)[[1]], dimnames(isarma)[[1]])
    
    if (p > 0) 
        psiph <- psiwtsAR(phi, maxlag1)
    if (q > 0) 
        psith <- -psiwtsAR(theta, maxlag1)
    if (ps > 0) 
        psips <- shift(psiwtsAR(phiseas, maxlag), period)
    if (qs > 0) 
        psits <- -shift(psiwtsAR(thetaseas, maxlag), period)
    
    if (p > 0 && ps > 0) {
        ppsmatrix <- matrix(0, nrow = p, ncol = ps)
        ppsmatrix[1, 1] <- crossprod(psiph[-(1:(period - 1))], psips[-((maxlag1 - period + 
            2):maxlag1)])
        if (ps > 1) 
            ppsmatrix[1, -1] <- sapply(2:ps, function(i) crossprod(psiph[-(1:(period * i - 
                1))], psips[-((maxlag1 - i * period + 2):maxlag1)]))
        if (p > 1) 
            ppsmatrix[-1, 1] <- sapply(2:p, function(i) crossprod(psiph[-(1:(period + i - 
                1))], psips[-((maxlag1 - i - period + 2):maxlag1)]))
        if (p > 1 && ps > 1) {
            for (i in 2:p) {
                for (j in 2:ps) {
                  ppsmatrix[i, j] <- crossprod(psiph[-(1:(period * j + i - 1))], psips[-((maxlag1 - 
                    i - period * j + 2):maxlag1)])
                }
            }
        }
    }
    if (q > 0 && ps > 0) {
        qpsmatrix <- matrix(0, nrow = q, ncol = ps)
        qpsmatrix[1, 1] <- crossprod(psith[-(1:(period - 1))], psips[-((maxlag1 - period + 
            2):maxlag1)])
        if (ps > 1) 
            qpsmatrix[1, -1] <- sapply(2:ps, function(i) crossprod(psith[-(1:(period * i - 
                1))], psips[-((maxlag1 - i * period + 2):maxlag1)]))
        if (q > 1) 
            qpsmatrix[-1, 1] <- sapply(2:q, function(i) crossprod(psith[-(1:(period + i - 
                1))], psips[-((maxlag1 - i - period + 2):maxlag1)]))
        if (q > 1 && ps > 1) {
            for (i in 2:q) {
                for (j in 2:ps) {
                  qpsmatrix[i, j] <- crossprod(psith[-(1:(period * j + i - 1))], psips[-((maxlag1 - 
                    i - period * j + 2):maxlag1)])
                }
            }
        }
    }
    if (p > 0 && qs > 0) {
        pqsmatrix <- matrix(0, nrow = p, ncol = qs)
        pqsmatrix[1, 1] <- crossprod(psiph[-(1:(period - 1))], psits[-((maxlag1 - period + 
            2):maxlag1)])
        if (qs > 1) 
            pqsmatrix[1, -1] <- sapply(2:qs, function(i) crossprod(psiph[-(1:(period * i - 
                1))], psits[-((maxlag1 - i * period + 2):maxlag1)]))
        if (p > 1) 
            pqsmatrix[-1, 1] <- sapply(2:p, function(i) crossprod(psiph[-(1:(period + i - 
                1))], psits[-((maxlag1 - i - period + 2):maxlag1)]))
        if (p > 1 && qs > 1) {
            for (i in 2:p) {
                for (j in 2:qs) {
                  pqsmatrix[i, j] <- crossprod(psiph[-(1:(period * j + i - 1))], psits[-((maxlag1 - 
                    i - period * j + 2):maxlag1)])
                }
            }
        }
    }
    
    if (q > 0 && qs > 0) {
        qqsmatrix <- matrix(0, nrow = q, ncol = qs)
        qqsmatrix[1, 1] <- crossprod(psith[-(1:(period - 1))], psits[-((maxlag1 - period + 
            2):maxlag1)])
        if (qs > 1) 
            pqsmatrix[1, -1] <- sapply(2:qs, function(i) crossprod(psith[-(1:(period * i - 
                1))], psits[-((maxlag1 - i * period + 2):maxlag1)]))
        if (q > 1) 
            qqsmatrix[-1, 1] <- sapply(2:q, function(i) crossprod(psith[-(1:(period + i - 
                1))], psits[-((maxlag1 - i - period + 2):maxlag1)]))
        if (q > 1 && qs > 1) {
            for (i in 2:q) {
                for (j in 2:qs) {
                  qqsmatrix[i, j] <- crossprod(psith[-(1:(period * j + i - 1))], psits[-((maxlag1 - 
                    i - period * j + 2):maxlag1)])
                }
            }
        }
    }
    block <- rbind(cbind(ppsmatrix, pqsmatrix), cbind(qpsmatrix, qqsmatrix))
    imatrix <- cbind(iarma, block)
    imatrix <- rbind(imatrix, cbind(if (!is.null(block)) 
        t(block) else NULL, isarma))
    dimnames(imatrix) <- list(inames, inames)
    return(imatrix)
    
}
"iARMASeasMult" <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), period = 0, warn = T) {
    if (period == 0) {
        if (warn) 
            warning("Seasonal call, but period = 0; reverting to nonseasonal model")
        return(iARMA(phi = phi, theta = theta))
    }
    if (!(InvertibleQ(phi) & InvertibleQ(theta) & InvertibleQ(phiseas) & InvertibleQ(thetaseas))) {
        warning("Model is non-causal or non-invertible\n")
        return(NULL)
    }
    imatrix <- iarma <- isarma <- block <- ppsmatrix <- qpsmatrix <- pqsmatrix <- qqsmatrix <- NULL
    p <- length(phi)
    q <- length(theta)
    ps <- length(phiseas)
    qs <- length(thetaseas)
    
    if (p + q > 0) 
        iarma <- iARMA(phi = phi, theta = theta)
    if (ps + qs > 0) 
        isarma <- iARMA(phi = phiseas, theta = thetaseas, period = period)
    inames <- c(dimnames(iarma)[[1]], dimnames(isarma)[[1]])
    phi.s <- shift(c(1, phiseas), period)[-1]
    theta.s <- shift(c(1, thetaseas), period)[-1]
    
    if (p > 0 && ps > 0) {
        ppsmatrix <- matrix(0, nrow = p, ncol = ps)
        tt <- tccfAR(phi, phi.s)
        for (i in 1:p) for (j in 1:ps) ppsmatrix[i, j] <- tt[ps * period + (i - 0) - (j - 
            0) * period]
    }
    if (q > 0 && ps > 0) {
        qpsmatrix <- matrix(0, nrow = q, ncol = ps)
        tt <- -tccfAR(theta, phi.s)
        for (i in 1:q) for (j in 1:ps) qpsmatrix[i, j] <- tt[ps * period + (i - 0) - (j - 
            0) * period]
    }
    if (p > 0 && qs > 0) {
        pqsmatrix <- matrix(0, nrow = p, ncol = qs)
        tt <- -tccfAR(phi, theta.s)
        for (i in 1:p) for (j in 1:qs) pqsmatrix[i, j] <- tt[qs * period + (i - 0) - (j - 
            0) * period]
    }
    
    if (q > 0 && qs > 0) {
        qqsmatrix <- matrix(0, nrow = q, ncol = qs)
        tt <- tccfAR(theta, theta.s)
        for (i in 1:q) for (j in 1:qs) qqsmatrix[i, j] <- tt[qs * period + (i - 0) - (j - 
            0) * period]
    }
    
    
    block <- rbind(cbind(ppsmatrix, pqsmatrix), cbind(qpsmatrix, qqsmatrix))
    imatrix <- cbind(iarma, block)
    imatrix <- rbind(imatrix, cbind(if (!is.null(block)) 
        t(block) else NULL, isarma))
    dimnames(imatrix) <- list(inames, inames)
    return(imatrix)
    
}

"iFARMA" <- function(phi = numeric(0), theta = numeric(0)) {
    
    if (!(InvertibleQ(phi) & InvertibleQ(theta))) {
        warning("Model is non-causal or non-invertible\n")
        return(NULL)
    }
    I22 <- (pi^2)/6
    if ((length(phi) == 0) && (length(theta) == 0)) {
        I22 <- as.matrix(I22)
        dimnames(I22) <- list("d.f", "d.f")
        return(I22)
    }
    I11 <- iARMA(phi = phi, theta = theta)
    J11 <- numeric(0)
    J12 <- numeric(0)
    
    if (length(phi) > 0) 
        J11 <- jFARMA(phi)
    if (length(theta) > 0) 
        J12 <- -jFARMA(theta)
    J <- c(J11, J12)
    I <- rbind(I11, J)
    J <- c(J, I22)
    I <- cbind(I, J)
    inames <- c(dimnames(I11)[[1]], "d.f")
    dimnames(I) <- list(inames, inames)
    return(I)
}
"iFARMAapprox" <- function(phi = numeric(0), theta = numeric(0)) {
    
    if (!(InvertibleQ(phi) & InvertibleQ(theta))) {
        warning("Model is non-causal or non-invertible\n")
        return(NULL)
    }
    I22 <- (pi^2)/6
    if ((length(phi) == 0) && (length(theta) == 0)) {
        I22 <- as.matrix(I22)
        dimnames(I22) <- list("d.f", "d.f")
        return(I22)
    }
    I11 <- iARMAapprox(phi = phi, theta = theta)
    J11 <- numeric(0)
    J12 <- numeric(0)
    
    if (length(phi) > 0) 
        J11 <- jFARMA(phi)
    if (length(theta) > 0) 
        J12 <- -jFARMA(theta)
    J <- c(J11, J12)
    I <- rbind(I11, J)
    J <- c(J, I22)
    I <- cbind(I, J)
    inames <- c(dimnames(I11)[[1]], "d.f")
    dimnames(I) <- list(inames, inames)
    return(I)
}

# write about below

"iARFIMA" <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
    period = 0, dfrac = TRUE, dfs = FALSE, exact = TRUE) {
    if (exact) 
        ans <- iARFIMAexact(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period, d = dfrac, dfs = dfs) else ans <- iARFIMAapprox(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
        period = period, d = dfrac, dfs = dfs)
    ans
}


"iARFIMAexact" <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), thetaseas = numeric(0), 
    period = 0, d = TRUE, dfs = T) {
    
    if (period == 0) {
        if (dfs || length(phiseas) > 0 || length(thetaseas) > 0) 
            stop("0 period, but seasonal parts")
        if (d) 
            return(iFARMA(phi, theta)) else {
            if ((length(phi) == 0) && (length(theta) == 0)) 
                return(numeric(0)) else return(iARMA(phi, theta))
        }
    }
    if (d && dfs) {
        I22 <- (pi^2)/6
        I33 <- (pi^2)/(6 * period)
        I4 <- cbind(c(I22, I33), c(I33, I22))
        if ((length(phi) == 0) && (length(theta) == 0) && (length(phiseas) == 0) && (length(thetaseas) == 
            0)) {
            dimnames(I4) <- list(c("d.f", paste("d.f.", period, sep = "")), c("d.f", paste("d.f.", 
                period, sep = "")))
            return(I4)
        }
        
        I11 <- iARMASeasMult(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        J11 <- numeric(0)
        if (length(phi) > 0) 
            J11 <- jFARMA(phi) * -1
        J12 <- numeric(0)
        if (length(phiseas) > 0) 
            J12 <- jFARMA(phiseas, period = period, fs = F) * -1
        if (length(theta) > 0) 
            J11 <- c(J11, jFARMA(theta))
        if (length(thetaseas) > 0) 
            J12 <- c(J12, jFARMA(thetaseas, period = period, fs = F))
        J <- c(J11, J12)
        K11 <- numeric(0)
        if (length(phi) > 0) 
            K11 <- jFARMA(phi, period = period, fs = T) * -1
        K12 <- numeric(0)
        if (length(phiseas) > 0) 
            K12 <- jFARMA(phiseas) * -1
        if (length(theta) > 0) 
            K11 <- c(K11, jFARMA(theta, period = period, fs = T))
        if (length(thetaseas) > 0) 
            K12 <- c(K12, jFARMA(thetaseas))
        K <- c(K11, K12)
        I <- rbind(I11, J, K)
        I <- cbind(I, rbind(cbind(J, K), I4))
        inames <- c(dimnames(I11)[[1]], "d.f", paste("d.f.", period, sep = ""))
        dimnames(I) <- list(inames, inames)
        return(I)
    } else if (d) {
        I22 <- (pi^2)/6
        if ((length(phi) == 0) && (length(theta) == 0) && (length(phiseas) == 0) && (length(thetaseas) == 
            0)) 
            return(I22)
        I11 <- iARMASeasMult(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        J11 <- numeric(0)
        if (length(phi) > 0) 
            J11 <- jFARMA(phi) * -1
        J12 <- numeric(0)
        if (length(phiseas) > 0) 
            J12 <- jFARMA(phiseas, period = period, fs = F) * -1
        if (length(theta) > 0) 
            J11 <- c(J11, jFARMA(theta))
        if (length(thetaseas) > 0) 
            J12 <- c(J12, jFARMA(thetaseas, period = period, fs = F))
        J <- c(J11, J12)
        I <- rbind(I11, J)
        I <- cbind(I, rbind(cbind(J), I22))
        inames <- c(dimnames(I11)[[1]], "d.f")
        dimnames(I) <- list(inames, inames)
        return(I)
    } else if (dfs) {
        I22 <- (pi^2)/6
        if ((length(phi) == 0) && (length(theta) == 0) && (length(phiseas) == 0) && (length(thetaseas) == 
            0)) {
            I22 <- as.matrix(I22)
            dimnames(I22) <- list(paste("d.f.", period, sep = ""), paste("d.f.", period, 
                sep = ""))
            return(I22)
        }
        I11 <- iARMASeasMult(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        J11 <- numeric(0)
        if (length(phi) > 0) 
            J11 <- jFARMA(phi, period = period, fs = T) * -1
        J12 <- numeric(0)
        if (length(phiseas) > 0) 
            J12 <- jFARMA(phiseas) * -1
        if (length(theta) > 0) 
            J11 <- c(J11, jFARMA(theta, period = period, fs = T))
        if (length(thetaseas) > 0) 
            J12 <- c(J12, jFARMA(thetaseas))
        J <- c(J11, J12)
        I <- rbind(I11, J)
        I <- cbind(I, rbind(cbind(J), I22))
        inames <- c(dimnames(I11)[[1]], paste("d.f.", period, sep = ""))
        dimnames(I) <- list(inames, inames)
        return(I)
    } else {
        I <- iARMASeasMult(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        return(I)
    }
}
"iARFIMAapprox" <- function(phi = numeric(0), theta = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), period = 0, d = T, dfs = T) {
    
    if (period == 0) {
        if (dfs || length(phiseas) > 0 || length(thetaseas) > 0) 
            stop("0 period, but seasonal parts")
        if (d) 
            return(iFARMAapprox(phi, theta)) else {
            if ((length(phi) == 0) && (length(theta) == 0)) 
                return(numeric(0)) else return(iARMAapprox(phi, theta))
        }
    }
    if (d && dfs) {
        I22 <- (pi^2)/6
        I33 <- (pi^2)/(6 * period)
        I4 <- cbind(c(I22, I33), c(I33, I22))
        if ((length(phi) == 0) && (length(theta) == 0) && (length(phiseas) == 0) && (length(thetaseas) == 
            0)) {
            dimnames(I4) <- list(c("d.f", paste("d.f.", period, sep = "")), c("d.f", paste("d.f.", 
                period, sep = "")))
            return(I4)
        }
        I11 <- iARMASeasMultapprox(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        J11 <- numeric(0)
        if (length(phi) > 0) 
            J11 <- jFARMA(phi) * -1
        J12 <- numeric(0)
        if (length(phiseas) > 0) 
            J12 <- jFARMA(phiseas, period = period, fs = F) * -1
        if (length(theta) > 0) 
            J11 <- c(J11, jFARMA(theta))
        if (length(thetaseas) > 0) 
            J12 <- c(J12, jFARMA(thetaseas, period = period, fs = F))
        J <- c(J11, J12)
        K11 <- numeric(0)
        if (length(phi) > 0) 
            K11 <- jFARMA(phi, period = period, fs = T) * -1
        K12 <- numeric(0)
        if (length(phiseas) > 0) 
            K12 <- jFARMA(phiseas) * -1
        if (length(theta) > 0) 
            K11 <- c(K11, jFARMA(theta, period = period, fs = T))
        if (length(thetaseas) > 0) 
            K12 <- c(K12, jFARMA(thetaseas))
        K <- c(K11, K12)
        I <- rbind(I11, J, K)
        I <- cbind(I, rbind(cbind(J, K), I4))
        inames <- c(dimnames(I11)[[1]], "d.f", paste("d.f.", period, sep = ""))
        dimnames(I) <- list(inames, inames)
        return(I)
    } else if (d) {
        I22 <- (pi^2)/6
        if ((length(phi) == 0) && (length(theta) == 0) && (length(phiseas) == 0) && (length(thetaseas) == 
            0)) {
            I22 <- as.matrix(I22)
            dimnames(I22) <- list(paste("d.f.", period, sep = ""), paste("d.f.", period, 
                sep = ""))
            return(I22)
        }
        I11 <- iARMASeasMultapprox(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        J11 <- numeric(0)
        if (length(phi) > 0) 
            J11 <- jFARMA(phi) * -1
        J12 <- numeric(0)
        if (length(phiseas) > 0) 
            J12 <- jFARMA(phiseas, period = period, fs = F) * -1
        if (length(theta) > 0) 
            J11 <- c(J11, jFARMA(theta))
        if (length(thetaseas) > 0) 
            J12 <- c(J12, jFARMA(thetaseas, period = period, fs = F))
        J <- c(J11, J12)
        I <- rbind(I11, J)
        I <- cbind(I, rbind(cbind(J), I22))
        inames <- c(dimnames(I11)[[1]], "d.f")
        dimnames(I) <- list(inames, inames)
        return(I)
    } else if (dfs) {
        I22 <- (pi^2)/6
        if ((length(phi) == 0) && (length(theta) == 0) && (length(phiseas) == 0) && (length(thetaseas) == 
            0)) 
            return(I22)
        I11 <- iARMASeasMultapprox(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        J11 <- numeric(0)
        if (length(phi) > 0) 
            J11 <- jFARMA(phi, period = period, fs = T) * -1
        J12 <- numeric(0)
        if (length(phiseas) > 0) 
            J12 <- jFARMA(phiseas) * -1
        if (length(theta) > 0) 
            J11 <- c(J11, jFARMA(theta, period = period, fs = T))
        if (length(thetaseas) > 0) 
            J12 <- c(J12, jFARMA(thetaseas))
        J <- c(J11, J12)
        I <- rbind(I11, J)
        I <- cbind(I, rbind(cbind(J), I22))
        inames <- c(dimnames(I11)[[1]], paste("d.f.", period, sep = ""))
        dimnames(I) <- list(inames, inames)
        return(I)
    } else {
        I <- iARMASeasMultapprox(phi = phi, theta = theta, phiseas = phiseas, thetaseas = thetaseas, 
            period = period)
        return(I)
    }
}
"jFARMA" <- function(theta, period = 1, fs = F, maxlag = 256) {
    q <- length(theta)
    J <- numeric(q)
    maxlag <- maxlag * period
    maxlagger <- 0:maxlag
    
    psis <- psiwtsAR(theta, maxlag = maxlag)
    
    if (fs) {
        for (k in 1:q) {
            
            a <- shift(psis[-(1:(period - k))], -period)
            J[k] <- sum(a/(1 + maxlagger[1:length(a)]))
        }
    } else {
        for (k in 1:q) {
            J[k] <- sum(psis/((k + maxlagger) * period))
        }
    }
    return(J)
}

"iARMA" <- function(phi = numeric(0), theta = numeric(0), period = 0) {
    ######## information matrix of arma#####################
    if (!(InvertibleQ(phi) & InvertibleQ(theta))) {
        cat("Model is non-causal or non-invertible\n")
        return(NULL)
    }
    p <- length(phi)
    q <- length(theta)
    unames <- character(0)
    vnames <- character(0)
    if (p > 0) {
        if (p > 1) {
            vvmatrix <- (tccfAR(phi, phi)[-(1:(p - 1))])[-(p + 1)]
        } else if (p == 1) {
            vvmatrix <- tccfAR(phi, phi)[-(p + 1)]
        }
        vvmatrix <- toeplitz(vvmatrix)
        imatrix <- vvmatrix
        if (!period) 
            vnames <- paste("phi(", 1:p, ")", sep = "") else vnames <- paste("phi.", period, "(", 1:p, ")", sep = "")
    }
    if (q > 0) {
        if (q > 1) {
            uumatrix <- (tccfAR(theta, theta)[-(1:(q - 1))])[-(q + 1)]
        } else if (q == 1) {
            uumatrix <- tccfAR(theta, theta)[-(q + 1)]
        }
        uumatrix <- toeplitz(uumatrix)
        imatrix <- uumatrix
        if (!period) 
            unames <- paste("theta(", 1:q, ")", sep = "") else unames <- paste("theta.", period, "(", 1:q, ")", sep = "")
    }
    if (p > 0 && q > 0) {
        uvmatrix <- matrix(numeric(1), nrow = p, ncol = q)
        tuv <- -tccfAR(phi, theta)
        for (i in 1:p) {
            for (j in 1:q) {
                uvmatrix[i, j] <- tuv[q + i - j]
            }
        }
        imatrix <- cbind(rbind(vvmatrix, t(uvmatrix)), rbind(uvmatrix, uumatrix))
    }
    
    inames <- c(vnames, unames)
    dimnames(imatrix) <- list(inames, inames)
    return(imatrix)
}

"tccfAR" <- function(phi, theta) {
    # auxilary function used with iarma######### computes the theoretical cross-covariance
    # function of two autoregressions z[t]-phi[1] z_[t-1] --- phi[p] z[t-p] = a[t]
    # z[t]-theta[1] z_[t-1] --- theta[q] z[t-q] = a[t] where p, q are length(phi),
    # length(theta)
    p <- length(phi)
    q <- length(theta)
    if (p == 0 || q == 0) 
        return(numeric(0))
    k <- p + q
    rhs <- c(-1, rep(0, k - 1))
    A <- matrix(numeric(k^2), nrow = k, ncol = k)
    for (i in 1:k) {
        for (j in 1:k) {
            imj <- i - j
            ijq <- i + j - q - 1
            if (i > q) {
                if (i > j && imj <= q) 
                  A[i, j] <- theta[imj] else if (i > q && imj == 0) 
                  A[i, j] <- -1
            } else {
                if (ijq > 0 && ijq <= p) 
                  A[i, j] <- phi[ijq] else if (ijq == 0) 
                  A[i, j] <- -1
            }
        }
    }
    
    return(solve(A, rhs))
} 
