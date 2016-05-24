
"grpSpeedFilter" <- function(x, speed.thr, window=5, ...)
{
    ## Value: Do stage one on matrix x (assuming it's a single unit),
    ## return a logical; whether each observation in x passed the test
    ## --------------------------------------------------------------------
    ## Arguments: x=matrix with cols: POSIXct, lon, lat; speed.thr=speed
    ## threshold (m/s), window=size of window to test; ...=arguments passed
    ## to distSpeed(), namely only 'method' for now.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (!window %% 2) stop ("window size must be an odd number")
    if (nrow(x) < window) stop ("there are fewer rows than window size")
    tpos <- window %/% 2                # test subscript - 1
    testfun <- function(k) {            # k=matrix with group to test
        mid <- tpos + 1                 # subscript of pt to test
        ref <- c(-seq(tpos), seq(tpos)) # subscripts of pts to test against
        mid.rep <- rep(mid, length(ref))
        speeds <- distSpeed(k[mid.rep, ], k[mid + ref, ], ...)[, 3]
        all(speeds > speed.thr, na.rm=TRUE) # TRUE if all speeds > thr
    }
    pass <- !logical(nrow(x))         # all pass at start up
    ## define all test rows and subscript for forward movement of window
    testrows <- seq(1 + tpos, nrow(x) - tpos); i <- 1
    for (j in testrows) {
        test <- testfun(x[c(i:(i + tpos - 1), j:(j + tpos)), ])
        if (test) {
            pass[j] <- FALSE
        } else {
            i <- i + 1
        }
    }
    pass
}


"rmsDistFilter" <- function(x, speed.thr, window=5, dist.thr, ...)
{
    ## Value: Run McConnell et al's filter and Austin et al's last stage,
    ## return 2-col matrix of logicals; whether each observation passed
    ## each test.
    ## --------------------------------------------------------------------
    ## Arguments: x=matrix with cols: POSIXct, lon, lat; speed.thr=speed
    ## threshold (m/s), window=size of window to test; dist.thr=distance
    ## threshold (km); ...=arguments passed to distSpeed(), namely only
    ## 'method' for now.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (!window %% 2) stop ("window size must be an odd number")
    wdw.errmess <- "there are fewer rows than window size"
    if (nrow(x) < window) stop (wdw.errmess)
    tpos <- window %/% 2                      # test subscript - 1
    ref <- c(-seq(tpos), seq(tpos))           # reference points for test
    travel.fun <- function(k) {           # k=subscripts
        xmid <- k[1]                      # 1st is the middle
        xmid.rep <- rep(xmid, length(k) - 1)
        others <- k[-1]                 # 1st against all other positions
        tr <- distSpeed(x[xmid.rep, ], x[others, ], ...)
        tr[, c(1, 3)]
    }
    testrows <- seq(tpos + 1, nrow(x) - tpos) # subscripts of locs to test
    rmsSwitch <- distPass <- !logical(length(testrows)) # to begin looping

    while (any(rmsSwitch)) {           # stop when switch is all FALSE
        switchidx <- which(rmsSwitch)
        if ((length(switchidx) + (tpos * 2)) < window) stop (wdw.errmess)
        testrows.new <- testrows[rmsSwitch]
        idx <- seq_along(testrows.new)  # index the above
        testidx.mtx <- testidx <- c(idx, sapply(ref, "+", idx))
        testidx.inner <- testidx >= 1 & testidx <= length(idx)
        testidx.mtx[testidx.inner] <- testrows.new[testidx[testidx.inner]]
        testidx.low <- testidx < 1
        testidx.mtx[testidx.low] <- testrows[1] - (1 - testidx[testidx.low])
        testidx.high <- testidx > length(idx)
        testidx.mtx[testidx.high] <- testrows[length(testrows)] +
            (testidx[testidx.high] - length(idx))
        testidx.mtx <- matrix(testidx.mtx, nrow=length(idx))
        travel <- apply(testidx.mtx, 1, travel.fun)
        if (dim(travel)[1] > 2) {
            dist.refs <- seq(length(ref))
            speed.refs <- seq(length(ref) + 1, nrow(travel))
        } else {
            dist.refs <- 1
            speed.refs <- 2
        }
        dists <- as.matrix(travel[dist.refs, ])
        speeds <- as.matrix(travel[speed.refs, ])
        ## root mean square value
        rms <- apply(speeds, 2, function(k) { # do this for every test group
            sqrt(sum(k ^ 2, na.rm=TRUE) / length(k))
        })
        ## prelim passing locs in current set
        rmsPass <- rms <= speed.thr
        ## find series of adj prelim failing locs in current set
        rmsPass.rle <- rle(rmsPass)
        bad <- which(!rmsPass.rle$values & rmsPass.rle$lengths > 1)
        if (length(bad) < 1) break  # stop looping if no adjc failing locs
        beg <- rev(length(rmsPass) + 1 - cumsum(rev(rmsPass.rle$lengths)))
        end <- cumsum(rmsPass.rle$lengths)
        for (j in bad) {
            bads <- seq(beg[j], end[j])
            Vi <- rms[bads]
            maxVi.idx <- which(Vi == max(Vi))
            peaks <- bads[maxVi.idx]
            ok <- bads[!bads %in% peaks]
            rmsPass[ok] <- TRUE
        }
        rmsSwitch[switchidx[!rmsPass]] <- FALSE # set failing locs in orig set
        ## Distance filter
        distt <- apply(dists, 2, function(k) sum(k, na.rm=TRUE) / length(k))
        distFail <- distt > dist.thr
        distPass[switchidx[distFail]] <- FALSE
    }

    rmsPass <- rmsSwitch
    ## Top and bottom ends pass as these can't be tested
    untested <- matrix(!logical(tpos * 2), ncol=2)
    rbind(untested, cbind(rmsPass, distPass), untested)
}


"austFilter" <- function(time, lon, lat, id=gl(1, 1, length(time)),
                         speed.thr, dist.thr, window=5, ...)
{
    ## Value: A matrix with logicals indicating whether each reading failed
    ## each filter.  This runs the filters in Austin et al. (2003).
    ## Results are presented from each filter, independently of the others
    ## for flexibility.
    ## --------------------------------------------------------------------
    ## Arguments: lat and lon=latitude and longitude vectors in degrees;
    ## time=POSIXct object with times for each point; id=factor identifying
    ## sections of the data to be treated separately; speed.thr=speed
    ## threshold (m/s); dist.thr=distance threshold (km); window=size of
    ## window to test; ...=arguments passed to grpSpeedFilter() and
    ## rmsDistFilter(), namely only 'method' for now.
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    ## FIRST STAGE ********************************************************
    locs <- data.frame(time, lon, lat)

    ## Do first stage over each seal's data, returns vector as long as locs
    first <- unlist(by(locs, id, grpSpeedFilter, speed.thr, window, ...),
                    use.names=FALSE)

    ## SECOND AND THIRD STAGES ********************************************
    good <- which(first)               # native subscripts that passed
    last <- do.call(rbind, by(locs[good, ], id[good], rmsDistFilter,
                              speed.thr, window, dist.thr, ...))
    filter123 <- cbind(firstPass=first,
                       secondPass=first, # 2nd and 3rd start the same as 1st
                       thirdPass=first)
    filter123[good, 2:3] <- last
    filter123
}

## TEST ZONE --------------------------------------------------------------
