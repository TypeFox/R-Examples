`cddews` <-
function (data, filter.number = 1, family = "DaubExPhase", switch = "direction", 
    correct = TRUE, verbose = FALSE, smooth = TRUE, sm.filter.number = 4, 
    sm.family = "DaubExPhase", levels = 3:6, type = "hard", policy = "LSuniversal", 
    by.level = FALSE, value = 0, dev = var) 
{
    now <- proc.time()[1:2]
    if (nrow(data) != ncol(data)) 
        stop(paste("Sorry, but imwd has only been coded for square images!"))
    data.wd <- imwd(data, filter.number = filter.number, family = family, 
        type = "station")
    RawPer <- getdata(data.wd, switch = switch)
    if (smooth == TRUE) {
        cat("Now starting to smooth\n")
        test <- RawPer
        NS <- dim(test)[1]
        for (i in 1:NS) {
            tmp <- test[i, , ]
            tmp.imwd <- imwd(tmp, filter.number = sm.filter.number, 
                family = sm.family)
            tmp.imwdTH <- threshold.imwd(tmp.imwd, levels = levels, 
                type = type, policy = policy, value = value, 
                by.level = by.level, dev = dev, compression = FALSE)
            tmp.imwr <- imwr(tmp.imwdTH)
            test[i, , ] <- tmp.imwr
        }
        RawPer <- test
        if (correct == FALSE) {
            cat("OK, so you've chosen to use the raw (uncorrected) periodogram!\n")
            l <- list(S = RawPer, datadim = dim(data), filter.number = filter.number, 
                family = "DaubExPhase", structure = switch, nlevels = data.wd$nlevels, 
                correct = correct, smooth = smooth, sm.filter.number = sm.filter.number, 
                sm.family = sm.family, levels = levels, type = type, 
                policy = policy, date = date())
        }
        if (correct == TRUE) {
            A <- D2Amat(-data.wd$nlevels, filter.number = data.wd$filter$filter.number, 
                family = data.wd$filter$family, switch = switch, 
                verbose = verbose)
            Ainv <- solve(A)
            first.last.d <- data.wd$fl.dbase$first.last.d
            first.last.c <- data.wd$fl.dbase$first.last.c
            firstD <- first.last.d[data.wd$nlevels, 1]
            lastD <- first.last.d[data.wd$nlevels, 2]
            LengthD <- lastD - firstD + 1
            LEVELS <- data.wd$nlevels
            TMP <- matrix(aperm(RawPer), nrow = 3 * LEVELS, ncol = LengthD^2,byrow = TRUE)
#            TMP <- matrix(c(t(RawPer)), nrow = 3 * LEVELS, ncol = LengthD^2,byrow = TRUE)
            TMP2 <- Ainv %*% TMP
            data2 <- array(0, dim(RawPer))
            for (i in (1:(3 * LEVELS))) {
                data2[i, , ] <- matrix(TMP2[i, ], nrow = LengthD, 
                  ncol = LengthD, byrow = TRUE)
            }
            speed <- proc.time()[1:2] - now
            cat("Took ", sum(speed), "seconds \n")
            l <- list(S = data2, datadim = dim(data), filter.number = filter.number, 
                family = "DaubExPhase", structure = switch, nlevels = data.wd$nlevels, 
                correct = correct, smooth = smooth, sm.filter.number = sm.filter.number, 
                sm.family = sm.family, levels = levels, type = type, 
                policy = policy, date = date())
        }
        class(l) <- "cddews"
        return(l)
    }
    if (smooth == FALSE) {
        if (correct == FALSE) {
            cat("OK, so you've chosen to use the raw (uncorrected periodogram!\n")
            l <- list(S = RawPer, datadim = dim(data), filter.number = filter.number, 
                family = "DaubExPhase", structure = switch, nlevels = data.wd$nlevels, 
                correct = correct, smooth = smooth, date = date())
        }
        if (correct == TRUE) {
            A <- D2Amat(-data.wd$nlevels, filter.number = data.wd$filter$filter.number, 
                family = data.wd$filter$family, switch = switch, 
                verbose = verbose)
            Ainv <- solve(A)
            first.last.d <- data.wd$fl.dbase$first.last.d
            first.last.c <- data.wd$fl.dbase$first.last.c
            firstD <- first.last.d[data.wd$nlevels, 1]
            lastD <- first.last.d[data.wd$nlevels, 2]
            LengthD <- lastD - firstD + 1
            LEVELS <- data.wd$nlevels
            TMP <- matrix(aperm(RawPer), nrow = 3 * LEVELS, ncol = LengthD^2, byrow = TRUE)
#            TMP <- matrix(c(t(RawPer)), nrow = 3 * LEVELS, ncol = LengthD^2, byrow = TRUE)
            TMP2 <- Ainv %*% TMP
            data2 <- array(0, dim(RawPer))
            for (i in (1:(3 * LEVELS))) {
                data2[i, , ] <- matrix(TMP2[i, ], nrow = LengthD, 
                  ncol = LengthD, byrow = TRUE)
            }
            speed <- proc.time()[1:2] - now
            cat("Took ", sum(speed), "seconds \n")
            l <- list(S = data2, datadim = dim(data), filter.number = filter.number, 
                family = "DaubExPhase", structure = switch, nlevels = data.wd$nlevels, 
                correct = correct, smooth = smooth, date = date())
        }
        class(l) <- "cddews"
        return(l)
    }
}

