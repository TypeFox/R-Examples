`pmden` <-
function (x, DISCR = FALSE, verbose = FALSE, bandwidth = -1, 
    extrema.nr = -1, accuracy = mad(x)/1000, extrema.mean = TRUE, 
    maxkuipnr = 19, asympbounds = FALSE, tolerance = 1e-08, localsq = TRUE, 
    locsq.factor = 0.9) 
{
    nsamp <- length(x)
    if (asympbounds || nsamp > max(kuipdiffbounds.x)) 
        currbounds <- kuipdiffbounds[length(kuipdiffbounds.x), 
            ] * sqrt(max(kuipdiffbounds.x))/sqrt(nsamp)
    else {
        currbounds <- double(maxkuipnr)
        for (i in 1:maxkuipnr) currbounds[i] <- approx(kuipdiffbounds.x, 
            kuipdiffbounds[, i], nsamp, rule = 2)$y
    }
    if (maxkuipnr > dim(kuipdiffbounds)[2]) 
        stop("maxkuipnr is too large")
    if (DISCR) {
        datax <- as.double(levels(as.factor(x)))
        N <- length(datax)
        dataemp <- c(0, as.double(summary(as.factor(x), maxsum = N))/nsamp)
        fdist.y <- cumsum(dataemp)
        x <- 0:N
        nsamp <- N + 1
    }
    else {
        dataemp <- c(0, rep(1/(nsamp - 1), nsamp - 1))
        x <- sort(x)
        if (min((x[-1] - x[-nsamp])/(x[nsamp] - x[1])) < 1e-14) {
            newx <- x / accuracy
            j <- 1
            while(j < nsamp)
              {
              ind <- (j:nsamp)[newx[j:nsamp]<floor(newx[j])+1]
              if(length(ind) == 0) # this can happen if newx[j] is extremely large
                k <- j
              else
                k <- max((j:nsamp)[newx[j:nsamp]<floor(newx[j])+1])

              newx[j:k] <- (1:(k-j+1))/(k-j+2)+floor(newx[j])
              j <- k+1
              }
            x <- newx * accuracy
        }
        fdist.y <- c(seq(0, 1, length = nsamp))
    }
    if (bandwidth > 0) 
        eps <- rep(bandwidth, nsamp)
    else {
        currprecision <- 0.5
        eps <- rep(0.5, nsamp)
    }
    eps[1] <- 0
    eps[nsamp] <- 0
    repeat {
        lower <- fdist.y - eps
        upper <- fdist.y + eps
        fts <- tautstring(x, fdist.y, lower, upper, upper[1], 
            lower[length(lower)], extrmean = extrema.mean)
        x.string <- fts$string
        if (sum(x.string[-1] != x.string[-(nsamp - 1)]) > 0) {
            ind1 <- min((1:(nsamp - 2))[x.string[-1] != x.string[-(nsamp - 
                1)]])
            ind2 <- max((1:(nsamp - 2))[x.string[-1] != x.string[-(nsamp - 
                1)]])
            if (x.string[ind1] > x.string[ind1 + 1]) 
                fts$nmax <- fts$nmax + 1
            if (x.string[ind2] < x.string[ind2 + 1]) 
                fts$nmax <- fts$nmax + 1
        }
        lastunif <- approx(fts$knotst, fts$knotsy, x)$y
        if (verbose) {
            par(mfrow = c(2, 1))
            if (DISCR) {
                plot(datax, dataemp[-1], col = "grey")
                lines(datax, x.string, col = "red")
            }
            else {
                hist(x, 40, prob = TRUE)
                lines((rep(x, rep(2, length(x))))[-c(1, 2 * length(x))], 
                  rep(x.string, rep(2, length(x.string))), col = "red")
            }
            plot(x, upper, type = "l")
            lines(fts$knotst, fts$knotsy, col = "red")
            lines(fts$knotst, fdist.y[fts$knotsind], col = "green")
            lines(x, lower)
        }
        if (bandwidth > 0) 
            break
        if (extrema.nr > 0) {
            if (fts$nmax > extrema.nr) 
                eps <- eps + currprecision
            if (currprecision < tolerance) {
                if (fts$nmax <= extrema.nr) 
                  break
            }
            else {
                currprecision <- currprecision/2
                eps[eps > 0] <- eps[eps > 0] - currprecision
            }
        }
        else {
            diff <- cumsum(dataemp) - lastunif
            currkkuip <- kkuip(diff, maxkuipnr)$met
            kuipinds <- c(currkkuip[1], currkkuip[-1] - currkkuip[-maxkuipnr]) > 
                currbounds + 1e-08
            if (sum(kuipinds) != 0) 
                eps[eps > 0] <- (currbounds[kuipinds])[1]/2
            else if (localsq) {
                irmax <- floor(log(nsamp)) + 3
                icomax <- 1
                prp <- exp(-1)
                currsum <- 2 * exp(-1)
                while (log(currsum) < log(0.95)/nsamp) {
                  icomax <- icomax + 1
                  prp <- prp/icomax
                  currsum <- currsum + prp
                }
                kni <- rep(0, nsamp)
                ind1 <- diff(kni, lag = 1) > 0.04/(nsamp^2)
                if (sum(ind1) > 0) 
                  kni[c(FALSE, ind1) | c(ind1, FALSE)] <- 1
                if (nsamp >= 3) {
                  ind2 <- diff(kni, lag = 2) > 0.25/(nsamp^1.5)
                  if (sum(ind2) > 0) 
                    kni[c(FALSE, FALSE, ind2) | c(FALSE, ind2, 
                      FALSE) | c(ind2, FALSE, FALSE)] <- 1
                }
                tmp <- .Fortran("denlocal", as.double(lastunif), 
                  kni = as.integer(kni), as.integer(nsamp), icomax = as.integer(icomax), 
                  irmax = as.integer(irmax), PACKAGE = "ftnonpar")
                if (sum(tmp$kni) == 0) 
                  break
                else eps[tmp$kni == 1] <- eps[tmp$kni == 1] * 
                  locsq.factor
            }
            else break
        }
        if (verbose) {
            print("Press Enter")
            readline()
        }
    }
    list(y = x.string, widthes = upper - fdist.y, nmax = fts$nmax, 
        ind = fts$knotsind, trans = lastunif)
}

