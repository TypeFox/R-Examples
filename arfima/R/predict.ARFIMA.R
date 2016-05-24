

predict.arfima <- function(object, n.ahead = 1, newxreg = NULL, predint = 0.95, bootpred = TRUE, 
    B = if (bootpred) 1000 else 0, trex = FALSE, seed = NA, setmuhat0 = FALSE, cpus = 1, 
    trend = NULL, ...) {
    npar <- bootpred
    if (!is.null(object$s) && any(object$s > 1)) 
        stop("predict only takes static regression parameters, not full transfer functions/dynamic regressions")
    if (!is.null(object$r) && any(object$r > 1)) 
        stop("predict only takes static regression parameters, not full transfer functions/dynamic regressions")
    myNROW <- function(x) if (is.null(x)) 
        0 else nrow(x)
    
    myNCOL <- function(x) if (is.null(x)) 
        0 else ncol(x)
    
    if (length(seed) == 0) 
        seed <- NA
    
    if (bootpred && is.na(seed)) 
        warning("bootstrap predictions and intervals may not be reproducible without setting the seed(s)")
    
    nrxreg <- myNROW(newxreg)
    
    if (nrxreg && !object$strReg)
      stop("Only predict with regular regression at this time (no transfer functions)")
    
    if (nrxreg && is.data.frame(newxreg)) {
      namx <- setdiff(object$namx, "Intercept")
      if(length(union(colnames(newxreg), namx)) != length(colnames(newxreg)))
        stop("Named arguments in xreg and newxreg are different")
      if(length(union(colnames(newxreg), namx)) != length(namx))
        stop("Named arguments in xreg and newxreg are different")
      newxreg <- newxreg[namx]
    }
    if (nrxreg) {
      newxreg <- cbind(rep(1, dim(newxreg[2])), as.matrix(newxreg))
    }

    
    if (nrxreg && is.null(object$xreg)) 
        stop("no xreg in input to arfima, but newxreg input to predict")
    
    if (!nrxreg && !is.null(object$xreg)) 
        stop("xreg in input to arfima, but no newxreg input to predict")
    
    if (myNCOL(object$xreg) != myNCOL(newxreg)) 
        stop("unconformable newxreg and xreg")
    
    if (myNROW(newxreg) && (myNROW(newxreg) != n.ahead)) 
        stop("unconformable 'newxreg' and 'n.ahead'")
    
    
    
    m <- length(object$modes)
    tacvfs <- tacvf(obj = object, xmaxlag = n.ahead, forPred = TRUE, n.ahead = n.ahead + 
        1)
    
    if (nrxreg || (bootpred && npar)) 
        resids <- resid(object)
    if (nrxreg) {
        newXreg <- newxreg
        res <- resids
        for (i in 1:m) res[[i]] <- as.vector(res[[i]]$regResiduals)
        y <- NULL
    } else y <- as.vector(object$z)
    dint <- object$dint
    dseas <- object$dseas
    period <- object$period
    
    ## multicore!!  call function inside function!
    preds <- vector("list", m)
    limiting <- FALSE
    for (i in 1:m) {
        zy <- if (!is.null(y)) 
            y else res[[i]]
        nz <- length(zy)
        zz <- c(zy, rep(0, n.ahead))
        muHat <- tacvfs[[i + 1]]$muHat
        if (object$differencing) {
            if (setmuhat0) 
                muHat <- 0
            icap <- dint + dseas * period
            zinit <- zy[(nz - icap + 1):nz]
            if (dint > 0) 
                zy <- diff(zy, differences = dint)
            if (dseas > 0) 
                zy <- diff(zy, differences = dseas, lag = period)
            if (nrxreg) {
                if (dint > 0) 
                  Xreg <- diff(Xreg, differences = dint)
                if (dseas > 0) 
                  Xreg <- diff(Xreg, differences = dseas, lag = period)
            }
        } else {
            icap <- 0
            zinit <- NULL
        }
        
        # do bootstrap preds in predictWork too!!!  so can eventually parallelize.
        
        znew <- predictWork(y = zy, r = tacvfs[[i + 1]]$tacvf, dint = dint, dseas = dseas, 
            period = period, muHat = muHat, n.ahead = n.ahead, trex = trex)
        nz <- length(zy)
        znew$Forecast <- znew$Forecast
        
        
        if (!is.null(trend)) {
            if (length(trend) < n.ahead) 
                stop("length of trend less than n.ahead")
            if (length(trend) != n.ahead) 
                warning("trend is not the same length as n.ahead")
            meanwx <- trend[1:n.ahead]
        } else meanwx <- rep(0, n.ahead)
        
        
        if (nrxreg) {
            meanwx <- meanwx + as.numeric(newXreg %*% object$modes[[i]]$omega)
        }
        
        znew$Forecast <- znew$Forecast + meanwx
        if (icap > 0) {
            znew$Forecast <- integ(z = znew$Forecast, zinit = zinit, dint = dint, dseas = dseas, 
                period = period)
            znew$Forecast <- znew$Forecast[-c(1:icap)]
        }
        
        if (bootpred && B > 0) {
            maxxer <- function(x) {
                x <- sort(x)[ceiling(length(x) * (1 - (1 - predint)/2))]
                x
            }
            minner <- function(x) {
                x <- sort(x)[floor(length(x) * ((1 - predint)/2))]
            }
            ## muHat correct in below?
            A <- Boot(object$modes[[i]], dint = dint, dseas = dseas, period = period, R = B, 
                n.ahead = n.ahead, n = nz, zinit = zinit, lastpoint = muHat, pred = TRUE, 
                seed = seed)
            up <- as.vector(apply(A, 1, maxxer))
            down <- as.vector(apply(A, 1, minner))
            znew$uppernp <- up + meanwx
            znew$lowernp <- down + as.numeric(meanwx)
            znew$meanvalnp <- as.vector(apply(A, 1, mean)) + meanwx
            
        } else {
            znew$lowernp <- znew$uppernp <- znew$meanvalnp <- znew$lowerp <- znew$upperp <- znew$meanvalp <- NULL
        }
        if (length(tacvfs[[i + 1]]$psis) > 0) {
            limiting <- TRUE
            sigmas <- cumsum(tacvfs[[i + 1]]$psis^2)[1:n.ahead]
            znew$limitVar <- sigmas * tacvfs[[i + 1]]$sigma2
            znew$limitSD <- sqrt(znew$limitVar)
        } else {
            znew$limitVar <- NULL
            znew$limitSD <- NULL
        }
        znew$sigma2 <- tacvfs[[i + 1]]$sigma2
        preds[[i]] <- znew
    }
    preds$z <- object$z
    preds$seed <- seed
    preds$limiting <- limiting
    preds$bootpred <- npar
    preds$B <- B
    preds$predint <- predint
    preds$name <- deparse(substitute(object))
    
    class(preds) <- "predarfima"
    preds
}

predictWork <- function(y, r, dint, dseas, period, muHat, n.ahead, trex = FALSE) {
    
    
    ny <- length(y)
    
    znew <- if (dseas + dint > 0) 
        TrFore(y, r, muHat, maxLead = n.ahead, getV = FALSE) else TrFore(y, r, muHat, maxLead = n.ahead)  #correct to use muHat, as it's from the differenced version
    
    coeffs <- NULL
    if (dint + dseas > 0) {
        if (trex) 
            coeffs <- rev(Z(l = n.ahead, d = dint, ds = dseas, s = period)) else coeffs <- wtsforexact(dint = dint, dseas = dseas, period = period, len = n.ahead)
    }
    
    znew$Forecast <- znew$Forecasts[1, ]
    znew$Forecasts <- NULL
    
    if (length(coeffs) > 0) {
        exactsig <- rep(0, n.ahead)
        P <- exactVals(r = r, n = ny, maxLead = n.ahead)
        for (i in 1:n.ahead) {
            exactsig[i] <- sum(sapply(1:i, function(u) sapply(1:i, function(v) coeffs[i - 
                u + 1] * coeffs[i - v + 1] * (r[abs(u - v) + 1] - P[u, v]))))
        }
    } else exactsig <- as.numeric(znew$SDForecasts[1, ]^2)
    
    znew$exactVar <- exactsig
    znew$exactSD <- sqrt(exactsig)
    
    znew
}


Z <- function(l, d, ds, s) {
    
    if ((d == 0) && (ds == 0)) 
        return(numeric(0))
    if (ds > 0 && s == 0) 
        stop("No period supplied")
    worker1 <- function(m, value, val) {
        if (m > 0) {
            val[m] <- val[m] + value
            if (d > 0) {
                for (i in 1:d) {
                  val <- worker1(m - i, value * choose(d, i) * (-1)^(i + 1), val)
                }
            }
            if (ds > 0) {
                for (j in 1:ds) {
                  val <- worker1(m - j * s, value * choose(ds, j) * (-1)^(j + 1), val)
                }
            }
            if (d > 0 && ds > 0) {
                for (i in 1:d) {
                  for (j in 1:ds) {
                    val <- worker1(m - i - j * s, value * choose(d, i) * choose(ds, j) * 
                      (-1)^(i + j + 1), val)
                  }
                }
            }
        }
        val
    }
    
    val <- numeric(l)
    val <- worker1(l, 1, val)
    
    val
}

print.predARFIMA <- function(x, digits = max(6, getOption("digits") - 3), ...) {
    n <- length(x$z)
    if (is.null(x$meanval)) 
        bootpred <- FALSE else bootpred <- TRUE
    if (bootpred) {
        lower <- x$lower
        upper <- x$upper
        meanpred <- x$meanval
        B <- x$B
        predint <- x$predint
    }
    forecasts <- x$Forecast
    exactSD <- x$exactSD
    approxSD <- x$approxSD
    n.ahead <- length(forecasts)
    ans <- rbind(forecasts, exactSD, approxSD)
    namer <- c("Forecasts", "Exact SD", if (!is.null(approxSD)) "Approximate SD" else NULL)
    rownames(ans) <- namer
    colnames(ans) <- 1:n.ahead
    if (bootpred) {
        ans1 <- rbind(upper, meanpred, lower)
        intt <- paste(round(100 * predint), "%", sep = "")
        namer1 <- c(paste("Upper", intt), "Prediction (Mean)", paste("Lower", intt))
        rownames(ans1) <- namer1
        colnames(ans1) <- 1:n.ahead
        ret <- list(`Forecasts and SDs` = ans, `Bootstrap Replicates` = B, `Bootstrap Predictions and Intervals` = ans1)
    } else ret <- list(`Forecasts and SDs` = ans)
    print(ret, digits = digits, ...)
}


TrFore <- function(z, r, zm, maxLead, getV = TRUE) {
    n <- length(z)
    
    if (length(r) < (n + maxLead - 1)) 
        stop("error: length(r) must be >= ", n + maxLead - 1)
    if (!(maxLead > 0)) 
        stop("maxLead must be > 0")
    
    zk <- z - zm
    zf <- vf <- matrix(numeric(maxLead), ncol = maxLead)
    gk <- t(matrix(r[n + 1 + outer(1:maxLead, 1:n, "-")], ncol = maxLead, byrow = TRUE))
    
    GI <- TrenchInverse(toeplitz(r[1:n]))
    gkGI <- crossprod(t(gk), GI)
    zf[1, ] <- zm + gkGI %*% zk
    if (getV) {
        for (j in 1:maxLead) vf[1, j] <- r[1] - sum(gkGI[j, ] * gk[j, ])
    }
    dimnames(zf) <- dimnames(vf) <- list(n, 1:maxLead)
    if (getV) 
        ans <- list(Forecasts = zf, SDForecasts = sqrt(vf)) else ans <- list(Forecasts = zf)
    ans
}

exactVals <- function(r, n, maxLead, way = 5) {
    if (length(r) < n + maxLead) 
        stop("r is too short")
    GI <- TrenchInverse(toeplitz(r[1:n]))
    vfs <- matrix(numeric(maxLead^2), ncol = maxLead)
    if (way == 1) {
        for (j in 1:(maxLead)) {
            gkj <- matrix(r[(n + j):(j + 1)], nrow = 1)
            for (h in 1:(maxLead)) {
                gkh <- matrix(r[(n + h):(h + 1)], ncol = 1)
                vfs[j, h] <- gkj %*% GI %*% gkh
            }
        }
    } else if (way == 2) {
        for (j in 1:(maxLead)) {
            gkj <- matrix(r[(n + j):(j + 1)], nrow = 1)
            for (h in 1:(maxLead)) {
                gkh <- matrix(r[(n + h):(h + 1)], ncol = 1)
                vfs[j, h] <- vfs[h, j] <- gkj %*% GI %*% gkh
            }
        }
    } else if (way == 3) {
        mat <- mat1 <- NULL
        for (j in 1:(maxLead)) {
            gkj <- matrix(r[(n + j):(j + 1)], nrow = 1)
            mat <- rbind(mat, gkj)
        }
        vfs <- mat %*% GI %*% t(mat)
    } else if (way == 4) {
        mat <- t(n + 1 + outer(1:maxLead, 1:n, "-"))
        gk <- matrix(r[mat], ncol = maxLead)
        vfs <- crossprod(t(crossprod(gk, GI)), gk)
    } else {
        ## seems to be fastest
        mat <- t(n + 1 + outer(1:maxLead, 1:n, "-"))
        gk <- t(matrix(r[mat], ncol = maxLead))
        vfs <- gk %*% GI %*% t(gk)
    }
    
    vfs
} 
