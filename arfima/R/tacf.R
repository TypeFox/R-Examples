mix <- function(x, y) {
    
    n <- 2 * length(x) - 2
    rev(Re(fft(fft(symtacvf(x)) * fft(symtacvf(y)), inverse = TRUE)/n)[(n/2 - 1):(n - 1)])
    
}

tacvf <- function(obj, xmaxlag = 0, forPred = FALSE, n.ahead = 0, ...) {
    if (xmaxlag < 0) 
        stop("xmaxlag must be >= 0")
    if (class(obj) == "arfima") {
        if (!obj$weeded) {
            warning("The object has not been weeded:  there may be duplicate modes, as well as excessive output. \n For more sensible output, please call weed() on the object")
        }
        m <- length(obj$modes)
        
        res <- vector("list", m + 1)
        name <- deparse(substitute(obj))
        res[[1]] <- name
        period <- obj$period
        dint <- obj$dint
        dseas <- obj$dseas
        maxlag <- obj$n - dint - dseas * period + xmaxlag - 1
        for (i in 1:m) {
            phi <- obj$modes[[i]]$phi
            theta <- obj$modes[[i]]$theta
            dfrac <- obj$modes[[i]]$dfrac
            phiseas <- obj$modes[[i]]$phiseas
            thetaseas <- obj$modes[[i]]$thetaseas
            dfs <- obj$modes[[i]]$dfs
            H <- obj$modes[[i]]$H
            Hs <- obj$modes[[i]]$Hs
            alpha <- obj$modes[[i]]$alpha
            alphas <- obj$modes[[i]]$alphas
            sigma2 <- obj$modes[[i]]$sigma2
            rr <- tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, 
                thetaseas = thetaseas, dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, 
                period = period, maxlag = maxlag, useCt = TRUE, sigma2 = sigma2)
            if (!forPred) 
                res[[i + 1]] <- rr else {
                res[[i + 1]]$muHat <- obj$modes[[i]]$muHat
                res[[i + 1]]$tacvf <- rr
                res[[i + 1]]$sigma2 <- obj$modes[[i]]$sigma2
                if (length(H) == 0 && length(Hs) == 0) {
                  ## need to do alpha, too.
                  if (n.ahead == 0) 
                    stop("invalid n.ahead")
                  res[[i + 1]]$psis <- psiwts(phi = phi, theta = theta, phiseas = phiseas, 
                    thetaseas = thetaseas, dfrac = dfrac, dfs = dfs, dint = dint, dseas = dseas, 
                    period = period, len = n.ahead)
                }
            }
        }
    } else if (class(obj) == "ARFIMA") {
        name <- deparse(substitute(obj))
        phi <- obj$phi
        theta <- obj$theta
        dfrac <- obj$dfrac
        phiseas <- obj$phiseas
        thetaseas <- obj$thetaseas
        dfs <- obj$dfs
        H <- obj$H
        Hs <- obj$Hs
        alpha <- obj$alpha
        alphas <- obj$alphas
        period <- obj$period
        maxlag <- length(obj$res) + xmaxlag - 1
        sigma2 <- obj$sigma2
        res <- vector("list", 2)
        res[[1]] <- name
        res[[2]] <- tacvfARFIMA(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, 
            thetaseas = thetaseas, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, 
            maxlag = maxlag, useCt = TRUE, sigma2 = sigma2)
        
    } else stop("The tacvf function is only defined for arfima or ARFIMA objects.")
    class(res) <- "tacvf"
    res
}

print.tacvf <- function(x, ...) {
    m <- length(x) - 1
    if (m > 1) {
        cat("Printing the tacvfs of", x[[1]], "\n")
        for (i in 2:(m + 1)) {
            cat("Mode", i - 1, "\n")
            if (is.vector(x[[i]], mode = "numeric")) {
                print(x[[i]])
            } else print(x[[i]]$tacvf)
        }
    } else {
        cat("Printing the tacvf of", x[[1]], "\n")
        if (is.vector(x[[2]], mode = "numeric")) {
            print(x[[2]])
        } else print(x[[2]]$tacvf)
    }
    invisible(x)
}


plot.tacvf <- function(x, type = "o", pch = 20, xlab = NULL, ylab = NULL, main = NULL, xlim = NULL, 
    ylim = NULL, tacf = FALSE, maxlag = NULL, lag0 = !tacf, ...) {
    if (is.null(ylab)) 
        ylab <- if (!tacf) 
            "tacvf" else "tacf"
    if (is.null(xlab)) 
        xlab <- "lag"
    
    m <- length(x) - 1
    for (i in 2:(m + 1)) {
        if (!is.vector(x[[i]], mode = "numeric")) 
            x[[i]] <- x[[i]]$tacvf
        if (tacf) 
            x[[i]] <- x[[i]]/x[[i]][1]
    }
    if (m > 9) {
        stop("more than nine modes currently not supported: please use the weed() function")
    }
    if (is.null(main)) {
        if (!tacf) {
            if (m > 1) 
                main <- paste("The tacvfs of ", x[[1]], sep = "") else main <- paste("The tacvf of ", x[[1]], sep = "")
        } else {
            if (m > 1) 
                main <- paste("The tacfs of ", x[[1]], sep = "") else main <- paste("The tacf of ", x[[1]], sep = "")
        }
    }
    if (!lag0 && is.null(xlim)) {
        xbds <- if (is.null(maxlag)) 
            1:(length(x[[2]]) - 1) else 1:maxlag
        xlim <- range(xbds)
    } else if (is.null(xlim)) {
        xbds <- if (is.null(maxlag)) 
            0:(length(x[[2]]) - 1) else 0:maxlag
        xlim <- range(xbds)
    } else {
        if (xlim[1] > xlim[2]) {
            warning("upper x limit lower than lower x limit: resetting")
            temp <- xlim[1]
            xlim[1] <- xlim[2]
            xlim[2] <- temp
        }
        
        if (xlim[2] > length(x[[2]]) - 1) {
            warning("the bounds of the x limits (lag) end at n - 1 + xmaxlag; subtracting one from each of the limits")
            xlim <- xlim - 1
        }
        if (xlim[1] < 0) {
            warning("lower x limit (lag) less than 0: resetting to 0")
            xlim[1] <- 0
        }
        xbds <- seq(xlim[1], xlim[2])
    }
    if (is.null(ylim)) {
        minner <- Inf
        maxxer <- -Inf
        for (i in 1:m) {
            minner <- min(minner, x[[i + 1]][xbds + 1])
            maxxer <- max(maxxer, x[[i + 1]][xbds + 1])
        }
        ylim <- c(minner, maxxer)
    } else if (ylim[1] > ylim[2]) {
        warning("upper y limit lower than lower y limit: resetting")
        temp <- ylim[1]
        ylim[1] <- ylim[2]
        ylim[2] <- temp
    }
    
    cols <- c("black", "red", "blue", "orange", "gray", "purple", "navy", "gold", "cyan")
    if (m > 1) {
        plot(x = xbds, y = x[[2]][xbds + 1], main = main, ylab = ylab, xlab = xlab, ylim = ylim, 
            xlim = xlim, type = type, pch = pch, col = cols[1])
        for (i in 1:(m - 1)) lines(x = xbds, y = x[[i + 2]][xbds + 1], type = type, pch = pch, 
            col = cols[i + 1])
        nncol <- floor(m/3)
        if (nncol == 0) 
            nncol <- 1
        legend(x = "topright", ncol = nncol, legend = paste("Mode", 1:m), pch = pch, col = cols[1:m])
    } else plot(x = xbds, y = x[[2]][xbds + 1], main = main, ylab = ylab, xlab = xlab, ylim = ylim, 
        xlim = xlim, type = type, pch = pch, col = cols[1])
}

tacfplot <- function(fits = list(), modes = "all", xlab = NULL, ylab = NULL, main = NULL, 
    xlim = NULL, ylim = NULL, maxlag = 20, lag0 = FALSE, ...) {
    
    if (length(fits) == 0) 
        stop("need at least one tacvf to compare")
    
    if (length(fits) == 1) {
        stop("for a single fit, only the fit should be provided: recommended to use plot.tacvf with tacf = TRUE")
    }
    
    if (class(fits) == "arfima") {
        if (modes != "all") {
            warning("only one fit provided and modes != 'all':  using all modes")
        }
        warning("tacfplot is intended to compare two fits on the same data: recommended to use plot.tacvf with tacf = TRUE")
        tacvf <- tacvf(fits)
        if (is.null(main)) 
            main <- "The tacf plot of the fit"
        
        plot(tacvf, type = "o", xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim, 
            tacf = TRUE, maxlag = maxlag, lag0 = lag0, ...)
    } else {
        if (modes != "all" && length(modes) != length(fits)) 
            stop("length of 'modes' not equal to length of 'fits', and modes != 'all'")
        numfits <- length(fits)
        if (numfits > 9) 
            stop("comparing more than 9 fits: stopping")
        z1 <- fits[[1]]$z
        num <- length(fits[[1]]$modes)
        for (i in 2:numfits) {
            z2 <- fits[[i]]$z
            if (length(z1) != length(z2) || !all(z1 == z2)) 
                stop("fits are not on the same series")
        }
        lenner <- length(tacvf(fits[[1]])[[1]])
        if (!lag0 && is.null(xlim)) {
            xbds <- if (is.null(maxlag)) 
                1:(lenner - 1) else 1:maxlag
            xlim <- range(xbds)
        } else if (is.null(xlim)) {
            xbds <- if (is.null(maxlag)) 
                0:(lenner - 1) else 0:maxlag
            xlim <- range(xbds)
        } else {
            if (xlim[1] > xlim[2]) {
                warning("upper x limit lower than lower x limit: resetting")
                temp <- xlim[1]
                xlim[1] <- xlim[2]
                xlim[2] <- temp
            }
            
            if (xlim[2] > length(z1) - 1) {
                warning("the bounds of the x limits (lag) end at n - 1; subtracting one from each of the limits")
                xlim <- xlim - 1
            }
            if (xlim[1] < 0) {
                warning("lower x limit (lag) less than 0: resetting to 0")
                xlim[1] <- 0
            }
            xbds <- seq(xlim[1], xlim[2])
        }
        
        tacfs <- vector("list", numfits)
        
        maxnum <- 1
        for (i in 1:numfits) {
            if (modes == "all") {
                maxnum <- max(maxnum, length(fits[[i]]$modes))
                tacfs[[i]] <- tacvf(fits[[i]])
                tacfs[[i]] <- tacfs[[i]][-1]
                for (j in 1:length(tacfs[[i]])) {
                  num <- num + 1
                  tacfs[[i]][[j]] <- tacfs[[i]][[j]]/tacfs[[i]][[j]][1]
                }
            } else {
                tacfs[[i]] <- tacfs[[i]][[modes[i]]]/tacfs[[i]][[modes[i]]][1]
            }
        }
        
        if (is.null(main)) {
            main <- paste("Comparing tacfs of ", numfits, " fits", sep = "")
        }
        namer <- paste("Fit ", 1:numfits)
        if (num > 9) 
            warning("more than 9 tacfs to plot:  plot will be very hard to read")
        
        if (is.null(xlab)) 
            xlab <- "lag"
        if (is.null(ylab)) 
            ylab <- "tacf"
        if (is.null(ylim)) {
            minner <- Inf
            maxxer <- -Inf
            for (i in 1:numfits) {
                if (modes == "all") {
                  for (j in 1:length(tacfs[[i]])) {
                    minner <- min(minner, tacfs[[i]][[j]][xbds + 1])
                    maxxer <- max(maxxer, tacfs[[i]][[j]][xbds + 1])
                  }
                } else {
                  minner <- min(minner, tacfs[[i]][xbds + 1])
                  maxxer <- max(maxxer, tacfs[[i]][xbds + 1])
                }
            }
            ylim <- c(minner, maxxer)
        } else if (ylim[1] > ylim[2]) {
            warning("upper y limit lower than lower y limit: resetting")
            temp <- ylim[1]
            ylim[1] <- ylim[2]
            ylim[2] <- temp
        }
        
        cols <- c("black", "red", "blue", "orange", "gray", "purple", "navy", "gold", "cyan")
        pchlist <- c(16, 17, 18, 0, 1, 2, 3, 4, 5)
        pchs <- pchlist[1:maxnum]
        
        if (modes != "all") {
            plot(x = xbds, y = tacfs[[1]][xbds + 1], main = main, ylab = ylab, xlab = xlab, 
                ylim = ylim, xlim = xlim, type = "o", pch = pchs[modes[1]], col = cols[1])
            for (i in 2:(numfits)) lines(x = xbds, y = tacfs[[i]][xbds + 1], type = "o", 
                pch = pchs[modes[i]], col = cols[i])
            legend(x = "topright", ncol = floor(numfits/3), legend = paste(namer, "mode", 
                modes), pch = pchs[modes], col = cols[1:numfits])
        } else {
            plot(x = xbds, y = tacfs[[1]][[1]][xbds + 1], main = main, ylab = ylab, xlab = xlab, 
                ylim = ylim, xlim = xlim, type = "o", pch = pchs[1], col = cols[1])
            if (length(tacfs[[1]]) > 1) {
                for (j in 2:length(tacfs[[1]])) lines(x = xbds, y = tacfs[[1]][[j]][xbds + 
                  1], type = "o", pch = pchs[j], col = cols[1])
            }
            for (i in 2:numfits) {
                for (j in 1:length(tacfs[[i]])) {
                  lines(x = xbds, y = tacfs[[i]][[j]][xbds + 1], type = "o", pch = pchs[j], 
                    col = cols[i])
                }
            }
            legend(x = "topright", ncol = floor(num/3), legend = c(namer, paste("Mode", 
                1:maxnum)), pch = c(rep(19, numfits), pchs), col = c(cols[1:numfits], rep("black", 
                maxnum)))
        }
    }
    
}

"tacvfARFIMA" <- function(phi = numeric(0), theta = numeric(0), dfrac = numeric(0), phiseas = numeric(0), 
    thetaseas = numeric(0), dfs = numeric(0), H = numeric(0), Hs = numeric(0), alpha = numeric(0), 
    alphas = numeric(0), period = 0, maxlag, useCt = T, sigma2 = 1) {
    if (length(maxlag) == 0) 
        stop("maxlag must be defined")
    if (!IdentInvertQ(phi = phi, theta = theta, dfrac = dfrac, phiseas = phiseas, thetaseas = thetaseas, 
        dfs = dfs, H = H, Hs = Hs, alpha = alpha, alphas = alphas, period = period, ident = F)) {
        return(NULL)
    }
    
    if (length(period) == 0 || is.na(period) || period == 0) {
        if (length(phiseas) > 0 && phiseas != 0) 
            stop("Period = 0, but seasonal phi != 0")
        if (length(thetaseas) > 0 && thetaseas != 0) 
            stop("Period = 0, but seasonal theta != 0")
        if ((length(dfs) > 0 && dfs != 0) || (length(Hs) > 0 && Hs != 0) || (length(alphas) > 
            0 && alphas != 0)) 
            stop("Period = 0, but fractional seasonal parameter != 0")
        return(tacvfFARMA(phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, 
            maxlag = maxlag, useCt = useCt, sigma2 = sigma2))
    }
    if ((period != round(period)) || (period < 2)) 
        stop("Period must be an integer >= 2")
    lagTrunc <- 2 * max(128, nextn(maxlag, factors = 2))
    iseas <- nextn(period, factors = 2)
    model <- tacvfFARMA(phi = phi, theta = theta, dfrac = dfrac, H = H, alpha = alpha, maxlag = (lagTrunc * 
        iseas), nolagtrunc = T, useCt = useCt, sigma2 = 1)
    modelseas <- tacvfFARMA(phi = phiseas, theta = thetaseas, dfrac = dfs, H = Hs, alpha = alphas, 
        maxlag = (lagTrunc * 2), nolagtrunc = T, useCt = useCt, sigma2 = 1)
    
    modelseas <- shift(modelseas, period, useCt = T)
    if (length(modelseas) < length(model)) 
        stop("oops")
    modelseas <- modelseas[1:(lagTrunc * iseas + 1)]
    
    z <- sigma2 * mix(model, modelseas)
    return(z[1:(maxlag + 1)])
}


"tacvfFARMA" <- function(phi = numeric(0), theta = numeric(0), dfrac = numeric(0), H = numeric(0), 
    alpha = numeric(0), maxlag, nolagtrunc = F, useCt = T, sigma2 = 1) {
    if (length(H) + length(dfrac) + length(alpha) > 1) 
        stop("more than one GHD process specified:  stopping")
    if ((length(H) == 0) && (length(dfrac) == 0) && (length(alpha) == 0) && (length(phi) == 
        0) && (length(theta) == 0)) {
        return(c(sigma2, rep(0, maxlag)))
    }
    
    if ((length(dfrac) == 0) && (length(H) == 0) && (length(alpha) == 0)) {
        return(tacvfARMA(phi = phi, theta = theta, maxlag = maxlag, useCt = useCt, sigma2 = sigma2))
    }
    
    onlyHD <- (length(phi) == 0 && length(theta) == 0)
    
    if (nolagtrunc || onlyHD) 
        lagTrunc <- maxlag else {
        lagTrunc <- 2 * max(128, nextn(maxlag, factors = 2))
    }
    
    
    if (length(dfrac) > 0) 
        x <- tacvfFDWN(dfrac = dfrac, maxlag = lagTrunc, useCt = useCt) else if (length(H) > 0) 
        x <- tacvfFGN(H = H, maxlag = lagTrunc, useCt = useCt) else x <- tacvfHD(alpha = alpha, maxlag = lagTrunc, useCt = useCt)
    
    if (onlyHD) 
        return(sigma2 * x[1:(maxlag + 1)])
    
    y <- tacvfARMA(phi = phi, theta = theta, maxlag = lagTrunc, useCt = useCt, sigma2 = 1)
    
    z <- sigma2 * mix(x, y)
    return(z[1:(maxlag + 1)])
}

"symtacvf" <- function(x) {
    c(rev(x[-1])[-1], x)
}

"tacvfARMA" <- function(phi = numeric(0), theta = numeric(0), maxlag = 20, useCt = T, sigma2 = 1) {
    
    if (useCt) {
        dd <- function(x, y, ml) .C("tacvfARMA_C", as.double(x), as.integer(length(x)), 
            as.double(y), as.integer(length(y)), as.integer(ml), res = double(ml + 1))$res
        return(sigma2 * dd(phi, theta, maxlag))
    }
    p <- length(phi)
    q <- length(theta)
    maxlagp1 <- maxlag + 1
    
    if (max(p, q) == 0) {
        return(c(sigma2, numeric(maxlag)))
    }
    r <- max(p, q) + 1
    b <- numeric(r)
    C <- numeric(q + 1)
    
    C[1] <- 1
    theta2 <- c(-1, theta)
    phi2 <- numeric(3 * r)
    phi2[r] <- -1
    if (p > 0) {
        phi2[r + 1:p] <- phi
    }
    if (q > 0) {
        for (k in 1:q) {
            C[k + 1] <- -theta[k]
            if (p > 0) {
                for (i in 1:min(p, k)) {
                  C[k + 1] <- C[k + 1] + phi[i] * C[k + 1 - i]
                }
            }
        }
    }
    
    for (k in 0:q) {
        for (i in k:q) {
            b[k + 1] <- b[k + 1] - theta2[i + 1] * C[i - k + 1]
        }
    }
    
    if (p == 0) {
        g <- c(b, numeric(maxlagp1))[1:maxlagp1]
        
        return(g)
    } else if (p > 0) {
        a <- matrix(numeric(r^2), ncol = r)
        for (i in 1:r) {
            for (j in 1:r) {
                if (j == 1) {
                  a[i, j] <- phi2[r + i - 1]
                } else if (j != 1) {
                  a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 2]
                }
            }
        }
        
        
        g <- solve(a, -b)
        
        if (length(g) <= maxlag) {
            g <- c(g, numeric(maxlagp1 - r))
            
            for (i in (r + 1):maxlagp1) {
                g[i] <- phi %*% g[i - 1:p]
            }
            
            return(sigma2 * g[1:maxlagp1])
        } else if (length(g) >= maxlagp1) {
            
            return(sigma2 * g[1:maxlagp1])
        }
    }
}

"tacvfFDWN" <- function(dfrac, maxlag, useCt = T) {
    
    if (!InvertibleD(dfrac)) {
        warning("Model is non-causal or non-invertible\n")
        return(NULL)
    }
    if (useCt) {
        ta <- function(ds, ma) .C("tacvfFDWN_C", as.double(ds), as.integer(ma), x = double(ma + 
            1))$x
        return(ta(dfrac, maxlag))
    }
    x <- numeric(maxlag + 1)
    x[1] <- gamma(1 - 2 * dfrac)/gamma(1 - dfrac)^2
    for (i in 1:maxlag) {
        x[i + 1] <- ((i - 1 + dfrac)/(i - dfrac)) * x[i]
    }
    return(x)
}

tacvfFGN <- function(H, maxlag, useCt = TRUE) {
    if (!InvertibleH(H)) 
        return(NULL)
    if (useCt) {
        tg <- function(H, ma) .C("tacfFGN_C", as.double(H), as.integer(ma), x = double(ma + 
            1))$x
        return(tg(H, maxlag))
    }
    h2 <- 2 * H
    r <- sapply(0:maxlag, function(k) 0.5 * (abs(k + 1)^h2 - 2 * abs(k)^h2 + abs(k - 1)^h2))
    return(r)
}

## will lgamma make it faster??
Zeta <- function(s, n = 20) {
    
    d <- rep(0, n + 1)
    
    d[1] <- n * (factorial(n - 1))/(factorial(n))
    
    for (i in 1:(n)) {
        d[i + 1] <- d[i] + n * factorial(n + i - 1) * 4^i/(factorial(n - i) * factorial(2 * 
            i))
    }
    
    zeta <- 0
    
    for (i in 0:(n - 1)) {
        zeta <- zeta + (-1)^i * (d[i + 1] - d[n + 1])/(i + 1)^s
    }
    
    zeta <- zeta * (-1)/(d[n + 1] * (1 - 2^(1 - s)))
    zeta
}
## Reference paper. (zetaalgorithm.)

tacvfHD <- function(alpha, maxlag, useCt = TRUE) {
    if (!InvertibleAlpha(alpha)) 
        return(NULL)
    if (useCt) {
        th <- function(alpha, ma) .C("tacfHD_C", as.double(alpha), as.integer(ma), x = double(ma + 
            1))$x
        return(th(alpha, maxlag))
    }
    val <- (-2 * Zeta(alpha))^(-1)
    r <- c(1, val * sapply(1:maxlag, function(k) k^(-alpha)))
    return(r)
} 
