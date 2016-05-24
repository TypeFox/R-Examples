`DLLoglikelihood` <-
function (r, z, useC = TRUE) 
{
    EPS <- .Machine$double.eps # 1+EPS==1, machine epsilon
    if (useC) {
        if (length(r) != length(z)) 
            stop("Error: arguments must be same length")
	stanQ <- 1
        out <- .C("DLResid", as.double(z), res = as.double(numeric(length(z))), 
            as.integer(length(z)), as.double(r), LogL = as.double(1), EPS, as.integer(stanQ),
            fault = as.integer(1), PACKAGE = "ltsa")
        fault <- out$fault
        if (fault == 1) 
            stop("error: sequence not p.d.")
        LogL <- out$LogL
    }
    else {
        n <- length(z)
        if (n != length(r)) 
            stop("arguments have unequal length")
        error <- numeric(n)
        sigmasq <- numeric(n)
        error[1] <- z[1]
        sigmasq[1] <- r[1]
        phi <- r[2]/r[1]
        error[2] <- z[2] - phi * z[1]
        sigmasqkm1 <- r[1] * (1 - phi^2)
        sigmasq[2] <- sigmasqkm1
        logg <- log(r[1]) + log(sigmasqkm1)
        for (k in 2:(n - 1)) {
            if (sigmasqkm1 < 0 || abs(sigmasqkm1) < EPS) 
                stop("r is not a p.d. sequence")
            phikk <- (r[k + 1] - phi %*% rev(r[2:k]))/sigmasqkm1
            sigmasqk <- sigmasqkm1 * (1 - phikk^2)
            phinew <- phi - phikk * rev(phi)
            phi <- c(phinew, phikk)
            sigmasqkm1 <- sigmasqk
            logg <- logg + log(sigmasqk)
            error[k + 1] <- z[k + 1] - crossprod(phi, rev(z[1:k]))
            sigmasq[k + 1] <- sigmasqk
        }
        S <- sum((error * error)/sigmasq)
        LogL <- c(-n/2 * log(S/n) - logg/2)
    }
    LogL
}

