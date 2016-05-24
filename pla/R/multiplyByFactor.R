.multiplyByFactor <- function(KP, factor, nms,
                              sampleLabels, testSampleLabels, nTestSamples) {
    if ((nTestSamples) > 1) {
        KP <- matrix(unlist(KP), ncol = 5)
        dimnames(KP)[[1]] <- sampleLabels[-1]
        dimnames(KP)[[2]] <- nms
        if (length(factor) == 1)
            factor <- ifelse(is.na(factor) | is.nan(factor),
                             1/ KP[, 3], factor)
        if (length(factor) == 1) {
            ## cat(" *** PhEur325: A f == 1 T > 1 *** \n")
            KP <- KP * factor
        } else {
            ## cat(" *** PhEur325: B f > 1  T > 1 *** \n")
            if (length(factor) == 3 * nTestSamples) {
                A <- matrix(rep(c(KP),
                                rep(length(factor) / nTestSamples, length(c(KP)))),
                            ncol = 5)
                replace <- is.na(factor) | is.nan(factor)
                factor[replace] <- (100 / A[, 3])[replace]
                B <- matrix(factor, nrow = length(factor), ncol = 5)
                KP <- matrix(A * B, ncol = 5)
                nm1 <- rep(testSampleLabels,
                           rep(length(factor) / nTestSamples, nTestSamples))
                nm2 <- names(factor)
                dimnames(KP)[[1]] <- paste(nm1, nm2, sep = "")
                dimnames(KP)[[2]] <- nms
            } else {
                sameLength <- length(factor) == length(sampleLabels[-1])
                if (sameLength)
                    sameLength <- all(names(factor) == sampleLabels[-1])
                if (sameLength) {
                    A <- matrix(rep(c(KP), rep(1, length(c(KP)))), ncol = 5)
                    B <- matrix(factor, nrow = length(c(KP)) * 1 / 5,
                                ncol = 5)
                    replace <- is.na(B[, 3]) | is.nan(B[, 3])
                    B[replace, 2] <- (100 / A[replace, 3])
                    B[replace, 3] <- (100 / A[replace, 3])
                    B[replace, 4] <- (100 / A[replace, 3])
                    KP <- matrix(A * B, ncol = 5)
                    dimnames(KP)[[1]] <- sampleLabels[-1]
                    dimnames(KP)[[2]] <- nms
                } else {
                    A <- matrix(rep(c(KP),
                                    rep(length(factor), length(c(KP)))),
                                ncol = 5)
                    B <- matrix(factor,
                                nrow = length(c(KP)) * length(factor) / 5,
                                ncol = 5)
                    replace <- is.na(B[, 3]) | is.nan(B[, 3])
                    B[replace, 2] <- (100 / A[replace, 3])
                    B[replace, 3] <- (100 / A[replace, 3])
                    B[replace, 4] <- (100 / A[replace, 3])
                    KP <- matrix(A * B, ncol = 5)
                    nm1 <- rep(testSampleLabels, rep(length(factor), nTestSamples))
                    nm2 <- rep(names(factor), nTestSamples)
                    dimnames(KP)[[1]] <- paste(nm1, nm2, sep = "")
                    dimnames(KP)[[2]] <- nms
                }
            }
        }
    } else {
        if (length(factor) == 1) {
            ## cat(" *** PhEur325: C f == 1 T==1 *** \n")
            factor <- ifelse(is.na(factor) | is.nan(factor), 1/ KP[3], factor)
            KP <- KP * factor
            ## nms <- c("exp(Width", "Lower", "Potency", "Upper", "exp(CM")
            ## names(KP) <- paste(nms, rep(testSampleLabels, 5),
            ##                    c(")", "", "", "", ")"),
            ##                    sep = "")
            names(KP) <- nms
        } else {
            ## cat(" *** PhEur325: D f > 1 T == 1 *** \n")
            A <- matrix(rep(c(KP), rep(length(factor), length(c(KP)))),
                        ncol = 5)
            replace <- is.na(factor) | is.nan(factor)
            factor[replace] <- (100 / A[, 3])[replace]
            B <- matrix(rep(factor), length(c(KP)) * length(factor) / 5,
                        ncol = 5)
            KP <- A * B
            dimnames(KP)[[2]] <- nms
            dimnames(KP)[[1]] <- names(factor)
            ## dimnames(KP)[[2]] <- nms
        }
    }
    return(KP)
}

