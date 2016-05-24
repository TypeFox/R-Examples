"multiwdwr" <-
function (y, thresh,firstwidth=1,dyad.factor=1.1) 
{
    .C("multiwdwr", y = as.double(y), as.integer(length(y)), 
        as.double(thresh),as.integer(firstwidth),as.double(dyad.factor),
        PACKAGE="ftnonpar")$y
}
"pmreg" <-
function (y, thr.const = 2.3, verbose = FALSE, extrema.nr = -1, 
    bandwidth = -1, sigma = -1, localsqueezing = TRUE, 
    squeezing.factor = 0.5,tolerance=1e-08,extrema.mean = TRUE,DYADIC=TRUE,
    dyad.factor = 1.1,POSTISO=TRUE) 
{
 if (extrema.nr > -1) localsqueezing <- FALSE
    nsamp <- length(y)
    x<-seq(0,1,len=nsamp)
    if (sigma < 0) {
        sigma <- mad((y[-1] - y[-nsamp])/sqrt(2))
        if (verbose) 
            print(c("sigma is ", sigma))
    }
    fdist <- c(0, cumsum(y))/nsamp
    fdistx <- seq(0,1,len=nsamp+1)
    if (bandwidth < 0) 
        d <- max(fdist) - min(fdist)
    else d <- bandwidth
    currprecision <- d
    eps <- rep(d, nsamp + 1)
    lower <- fdist - d
    upper <- fdist + d
    repeat {
        tstring <- tautstring(fdistx, fdist, lower, upper, 0, fdist[nsamp + 1],extrmean=FALSE)
        y.string <- tstring$string
        if((bandwidth<0)&&(extrema.nr<0))
          {
          residuals <- y - y.string
          residuals <- residuals - mean(residuals)
          if(DYADIC) 
            residuals.wr <- multiwdwr(residuals, sqrt(thr.const * log(nsamp)) * sigma,dyad.factor=dyad.factor)
          else
            residuals.wr <- nondymwdr(residuals, sqrt(thr.const * log(nsamp)) * sigma)
          }
        if (verbose) {
            par(mfrow = c(2, 2))
            plot(fdistx, lower, type = "l", ylim = range(c(lower, 
                upper)))
            lines(fdistx, upper)
            if (length(tstring$knotsy) > 0) {
                lines(tstring$knotst, tstring$knotsy, col = "red")
            }
            plot(x, y, col = "grey")
            lines(x, y.string, col = "red")
            if((bandwidth<0)&&(extrema.nr<0))
              plot(x, residuals.wr, type = "l", col = "green")
            print("Press Enter")
            dum <- readline()
        }
        if(bandwidth>0) break
        if(extrema.nr>0)
          {
          if(tstring$nmax>extrema.nr)
            eps<-eps+currprecision 
          if(currprecision<tolerance)
            {
            if(tstring$nmax<=extrema.nr)
            break
            }
          else
            {
            currprecision<-currprecision/2
            eps<-eps-currprecision
            }
          }
        else
          {
          ind <- (abs(residuals.wr) > 1e-10)
          ind2 <- c(FALSE, ind) | c(ind, FALSE)
          if (length(ind[ind == TRUE]) == 0) 
            {
            tstring.bak <- tstring
            while(tstring$nmax == tstring.bak$nmax)
              {
              tstring.bak <- tstring
              eps <- eps * squeezing.factor
              lower <- fdist - eps
              upper <- fdist + eps
              tstring <- tautstring(fdistx, fdist, lower, upper, 0, fdist[nsamp + 1], extrmean = FALSE)
              }
            tstring <- tstring.bak
            y.string <- tstring$string
            break
            }
          if (localsqueezing) 
            eps[ind2] <- eps[ind2] * squeezing.factor
          else
            eps <- eps * squeezing.factor
          }
        lower <- fdist - eps
        upper <- fdist + eps
        }
    if(extrema.mean)
      y.string <- settomean(y.string,y)
    if(POSTISO)
      y.string <- isoanti(y,y.string)

    list(y = y.string, sigma = sigma, widthes = upper - fdist, 
        nmax = tstring$nmax, knotsind = tstring$knotsind, knotsy = tstring$knotsy)
}

"tautstring" <-
function (ttt, fdist, y.low, y.up, y1 = 0.5 * (y.low[1] + y.up[1]), 
    yn = 0.5 * (y.low[length(y.low)] + y.up[length(y.up)]),extrmean=TRUE)  
{
        tmp <- .C("tautstring", as.double(fdist), as.double(ttt), 
            as.double(y.low), as.double(y.up), as.double(y1), 
            as.double(yn), as.integer(length(y.low)), tautstring = double(length(y.low) - 
                1), knotsind = integer(length(y.low)), knotst = double(length(y.low)), 
            knotsy = double(length(y.low)), nknots = integer(1), 
            nmax=integer(1),extrmean=as.integer(extrmean),PACKAGE="ftnonpar")
        list(string = tmp$tautstring, knotsind = tmp$knotsind[1:tmp$nknots], 
            knotst = tmp$knotst[1:tmp$nknots], knotsy = tmp$knotsy[1:tmp$nknots], 
            nknots = tmp$nknots,nmax=tmp$nmax)
}

"pmlogreg" <-
function (y, thr.const = 2.5, verbose = FALSE, extrema.nr = -1, bandwidth = -1, 
    localsqueezing = TRUE, squeezing.factor = 0.5, tolerance = 0.001,extrema.mean=TRUE) 
{
    if (extrema.nr > -1) localsqueezing <- FALSE
    nsamp <- length(y)
    x <- seq(0, 1, len = nsamp)
    fdist <- c(0, cumsum(y))/nsamp
    fdistx <- seq(0, 1, len = nsamp + 1)
    if (bandwidth < 0) 
        d <- 0.5 * (max(fdist) - min(fdist))
    else d <- bandwidth
    currprecision <- d
    eps <- rep(d, nsamp + 1)
    lower <- fdist - d
    upper <- fdist + d
    repeat {
        tstring <- tautstring(fdistx, fdist, lower, upper, 0, 
            fdist[nsamp + 1],extrmean=extrema.mean)
        y.string <- tstring$string
        if ((bandwidth < 0) && (extrema.nr < 0)) {
            residuals <- (y - y.string)/(sqrt(max(0.000001,y.string*(1-y.string))))
            residuals.wr <- multiwdwr(residuals, sqrt(thr.const * 
                log(nsamp)) )
        }
        if (verbose) {
            par(mfrow = c(2, 2))
            plot(fdistx, lower, type = "l", ylim = range(c(lower, 
                upper)))
            lines(fdistx, upper)
            if (length(tstring$knotsy) > 0) {
                lines(tstring$knotst, tstring$knotsy, col = "red")
            }
            plot(x, y, col = "grey")
            lines(x, y.string, col = "red")
            if ((bandwidth < 0) && (extrema.nr < 0)) 
                plot(x, residuals.wr, type = "l", col = "green")
            print("Press Enter")
            dum <- readline()
        }
        if (bandwidth > 0) 
            break
        if (extrema.nr > 0) {
            if (tstring$nmax > extrema.nr) 
                eps <- eps + currprecision
            if (currprecision < tolerance) {
                if (tstring$nmax <= extrema.nr) 
                  break
            }
            else {
                currprecision <- currprecision/2
                eps <- eps - currprecision
            }
        }
        else {
            ind <- (abs(residuals.wr) > 1e-10)
            ind2 <- c(FALSE, ind) | c(ind, FALSE)
            if (length(ind[ind == TRUE]) == 0) 
                break
            if (localsqueezing) 
                eps[ind2] <- eps[ind2] * squeezing.factor
            else eps <- eps * squeezing.factor
        }
        lower <- fdist - eps
        upper <- fdist + eps
    }
    list(y = y.string, widthes = upper - fdist, 
        nmax = tstring$nmax, knotsind = tstring$knotsind, knotsy = tstring$knotsy)
}
"isoanti" <-
function(x,fhat){

        n <- length(x)
        kni <- c(1,(2:n)[fhat[-1]!=fhat[-n]],n+1)     
        nkn <- length(kni) - 1
        if (nkn >= 3) {
            kny <- fhat[kni[1:nkn]]
            currkn <- 1
            while (currkn < nkn) {
                currsign <- sign(kny[currkn + 1] - kny[currkn])
                nextkn <- currkn + 1
                while ((nextkn < nkn) && ((currsign * (kny[nextkn +
                  1] - kny[nextkn]) > 0))) nextkn <- nextkn +
                  1
                if (nextkn - currkn > 1) {
                  if (currsign == 1)
                    tmp <- isoreg(x[kni[currkn]:(kni[nextkn +
                      1] - 1)])$yf
                  else tmp <- -isoreg(-x[kni[currkn]:(kni[nextkn +
                    1] - 1)])$yf
                  maxy <- max(fhat[kni[currkn]:(kni[nextkn +
                    1] - 1)])
                  miny <- min(fhat[kni[currkn]:(kni[nextkn +
                    1] - 1)])
                  if (max(tmp) > maxy)
                    tmp[tmp > maxy] <- maxy
                  if (min(tmp) < miny)
                    tmp[tmp < miny] <- miny
                  fhat[kni[currkn]:(kni[nextkn + 1] - 1)] <- tmp
                }
                currkn <- nextkn
            }
        }
fhat
}
.First.lib <- function(lib, pkg) {
  if(version$major==0)
    stop("This version for R 1.00 or later")
  library.dynam("ftnonpar", pkg, lib)
}
