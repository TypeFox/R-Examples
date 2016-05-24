qqstuff <-
function (intensity, binsize = 4, paths = 100, respaths = 1000, 
    plot.q = FALSE, plot.sq = FALSE) 
{
    lambda <- binsize * intensity
    qqinfo <- list()
    vmat <- matrix(0, respaths, length(intensity))
    Av <- vmat
    Bv <- vmat
    bfv <- vmat
    Asqres <- matrix(0, 1, length(intensity))
    Bsqres <- Asqres
    bfsqres <- matrix(0, 1, length(intensity))
    for (i in 1:length(intensity)) {
        vmat[, i] <- rbinom(respaths, binsize, intensity[i])
    }
    for (j in 1:respaths) {
        Av[j, ] <- ansc(vmat[j, ], binsize)
        Bv[j, ] <- free(vmat[j, ], binsize)
        bfv[j, ] <- binhf.wd(vmat[j, ], binsize = binsize)$transformed
    }
    Al <- ansc(lambda, binsize)
    Bl <- free(lambda, binsize)
    bfl <- binhf.wd(lambda, binsize = binsize)$transformed
    vminusl <- vmat - matrix(lambda, nrow = nrow(vmat), ncol = ncol(vmat), 
        byrow = TRUE)
    AvminusAl <- (Av - matrix(Al, nrow = nrow(Av), ncol = ncol(Av), 
        byrow = TRUE))
    BvminusBl <- (Bv - matrix(Bl, nrow = nrow(Bv), ncol = ncol(Bv), 
        byrow = TRUE)) * sqrt(1 * (binsize + 0.5))
    bfvminusbfl <- (bfv - matrix(bfl, nrow = nrow(bfv), ncol = ncol(bfv), 
        byrow = TRUE))
    Asqres <- apply(AvminusAl^2, 2, mean)
    Bsqres <- apply(BvminusBl^2, 2, mean)
    bfsqres <- apply(bfvminusbfl^2, 2, mean)
    vqqnorm <- apply(vminusl[1:paths, ], 1, qqnormy)
    vqqmean <- apply(vqqnorm, 1, mean)
    Aqqnorm <- apply(AvminusAl[1:paths, ], 1, qqnormy)
    Aqqmean <- apply(Aqqnorm, 1, mean)
    Bqqnorm <- apply(BvminusBl[1:paths, ], 1, qqnormy)
    Bqqmean <- apply(Bqqnorm, 1, mean)
    bfqqnorm <- apply(bfvminusbfl[1:paths, ], 1, qqnormy)
    bfqqmean <- apply(bfqqnorm, 1, mean)
    x <- qnorm((1:length(intensity) - 0.5)/length(intensity))
    if (plot.q == TRUE) {
        getOption("device")()
        plot(sort(x), vqqmean, xlab = "Quantiles of Standard normal", 
            ylab = "Mean raw quantiles")
        getOption("device")()
        plot(sort(x), Aqqmean, xlab = "Quantiles of Standard normal", 
            ylab = "Mean Anscombe quantiles")
        getOption("device")()
        plot(sort(x), Bqqmean, xlab = "Quantiles of Standard normal", 
            ylab = "Mean Freeman-Tukey quantiles")
        getOption("device")()
        plot(sort(x), bfqqmean, xlab = "Quantiles of Standard normal", 
            ylab = "Mean Nunes-Nason quantiles")
        getOption("device")()
        plot(sort(x), vqqmean, xlab = "Quantiles of Standard normal", 
            ylab = "Mean raw quantiles")
        points(sort(x), Aqqmean, col = 2)
        points(sort(x), bfqqmean, col = 3)
        points(sort(x), Bqqmean, col = 4)
        abline(b = 1, a = 0)
    }
    if (plot.sq == TRUE) {
        r <- range(c(Asqres, Bsqres, bfsqres))
        getOption("device")()
        plot(Asqres, xlab = "", ylab = "Squared residuals (Anscombe)", 
            type = "l", ylim = r)
        getOption("device")()
        plot(Bsqres, xlab = "", ylab = "Squared residuals (Freeman-Tukey)", 
            type = "l", ylim = r)
        getOption("device")()
        plot(bfsqres, type = "l", xlab = "", ylab = "Squared residuals (Nunes-Nason)", 
            ylim = r)
    }
    qqinfo[[1]] <- vmat
    qqinfo[[2]] <- Av
    qqinfo[[3]] <- Bv
    qqinfo[[4]] <- bfv
    qqinfo[[5]] <- vminusl
    qqinfo[[6]] <- AvminusAl
    qqinfo[[7]] <- BvminusBl
    qqinfo[[8]] <- bfvminusbfl
    qqinfo[[9]] <- Asqres
    qqinfo[[10]] <- Bsqres
    qqinfo[[11]] <- bfsqres
    return(qqinfo)
}

