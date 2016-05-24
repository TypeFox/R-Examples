pheur325 <- function (data,
                      lmInteraction = NULL,
                      lmRestricted = NULL,
                      dr = 2,
                      response = "Response",
                      sampleLabels = levels(unlist(data["Sample"])),
                      indexOfReference = 1,
                      StdName = sampleLabels[indexOfReference],
                      sampleStepName = "SampleStep",
                      dfAdj = 0,
                      factor = 1,
                      alpha = 0.05) {
    nTestSamples <- length(sampleLabels) - 1
    testSampleLabels <- sampleLabels[-indexOfReference]
    testSampleLabelsAlt <- paste0("", testSampleLabels, " - ")
    Means <- sapply(split(data[response], data[sampleStepName], drop = TRUE), 
                    FUN = function(x) mean(t(x), na.rm = TRUE))
    nms <- names(Means)
    splitNms <- strsplit(nms, ":")
    SAMPLE <- unlist(lapply(splitNms,
                            FUN = function (x)
                            ifelse(length(x) > 1, x[1],
                                   substr(x, 1, 1)
                                   )))
    DILUTION <- as.numeric(unlist(lapply(splitNms,
                                         FUN = function (x)
                                         ifelse(length(x) > 1, x[2],
                                                substr(x, 2, nchar(x))
                                                ))))
    ## d <- max(DILUTION)
    d <- length(unique((DILUTION)))
    names(d) <- "#Doses"
    n <- dim(data)[1] / length(DILUTION)
    names(n) <- "#Replicates"
    h <- length(DILUTION) / d
    names(h) <- "#Samples"
    useLm <- FALSE
    if (any(is.na(data[response]))
        | length(which(!is.na(data[response]))) < h * d * n) {
        useLm <- TRUE
        namesCoefs <- names(coef(lmRestricted))
        i <- (substr(namesCoefs, 1, 14) == "factor(Sample)") &
            unlist(lapply(strsplit(namesCoefs, ":"), length)) < 2
        gammas <- coef(lmRestricted)[i]
        beta <- coef(lmRestricted)["Z"]
        MTcomplete <- gammas / beta
        names(MTcomplete) <- paste(rep("Log(Potency)", nTestSamples),
                                   testSampleLabels, sep = "")
        MT <- MTcomplete
        b  <- beta
    } else {
        X <- matrix(Means, nrow = d,
                    dimnames = list(paste("Dose", 1:d), sampleLabels))
        P <- apply(X, 2, sum)
        D <- matrix(rep(1:d, length(sampleLabels)), nrow = d,
                    dimnames = list(paste("Dose", 1:d), sampleLabels))
        L <- apply(X * D, 2, sum) - (d + 1)/2 * P
        H <- c(n/d, 12 * n/(d^3 - d))
        names(H) <- c("P", "L")
        b <- H["L"] * (sum(L))/(log(dr) * n * h)
        names(b) <- NULL
        MT <- (P[names(P) != StdName] - P[names(P) == StdName])/(d * b)
        names(MT) <- paste(rep("Log(Potency)", nTestSamples), testSampleLabels,
                           sep = "")
    }
    Anova <- anova(lmInteraction)
    SStot <- sum(Anova[, "Sum Sq"])
    SSsmp <- Anova["factor(Sample)", "Sum Sq"]
    SSreg <- Anova["Dilution", "Sum Sq"]
    SSres <- Anova["Residual", "Sum Sq"]
    DFres <- Anova["Residual", "Df"] - dfAdj
    DFsmp <- Anova["factor(Sample)", "Df"]
    DFreg <- Anova["Dilution", "Df"]
    DFtot <- sum(Anova[, "Df"])
    R     <- sqrt((SSsmp + SSreg) / SStot)
    Radj  <- sqrt(1 - DFtot * (SStot-SSsmp-SSreg) /
                  ((DFtot-DFsmp-DFreg) * SStot))
    s2    <- SSres/DFres
    qt    <- qt(1-alpha/2, DFres)
    C     <- SSreg/(SSreg - s2 * qt^2)
    names(C) <- NULL
    CM <- C * MT
    names(CM) <- paste(rep("CM", nTestSamples), testSampleLabels)
    if (useLm) {
        ## if (nTestSamples > 1)
        ##     for (i in 1:nTestSamples)
        ##         F <- Fiellers(lmRestricted, Which = i, alpha = alpha)
        F <- Fiellers(lmRestricted, alpha = alpha)
        if (nTestSamples == 1) {
            W <- (F[3] - F[1]/2)
            K <- c(W, F[1], MT, F[3], CM)
        } else {
            W <- c(F[,3] - F[,1])/2
            K <- cbind(W, F[,1], MT, F[,3], CM)
        }
        V0 <- (W^2/(C - 1) - CM^2)/2
        nxd <- SSreg / (V0 * b^2)
        V <- NA
    } else {
        V <- SSreg/(b^2 * d * n)
        W <- sqrt((C - 1) * (CM^2 + 2 * V))
        K <- c(W, CM - W, MT, CM + W, CM)
    }
    names(V) <- NULL
    names(W) <- paste(rep("Width, Conf.int., Log", nTestSamples), testSampleLabels)
    nms <- c("Width", "Log(Lower)", "Log(Potency)", "Log(Upper)", "CM")
    names(K) <- paste(rep(nms, rep(nTestSamples, 5)), rep(testSampleLabelsAlt, 5),
                      sep = "")
    if ((nTestSamples) > 1) {
        K <- matrix(unlist(K), ncol = 5)
        dimnames(K)[[1]] <- testSampleLabels
        dimnames(K)[[2]] <- nms
    }
    nms <- c("exp(Width)", "Lower", "Potency", "Upper", "exp(CM)")
    KP <- exp(K)
    if ((nTestSamples) > 1) {
        KP <- matrix(unlist(KP), ncol = 5)
        dimnames(KP)[[1]] <- testSampleLabels
        dimnames(KP)[[2]] <- nms
    } else
        names(KP) <- nms
    KP <- .multiplyByFactor(KP, factor, nms,
                            sampleLabels, testSampleLabelsAlt, nTestSamples)
    if (useLm)
        inputList <- list(Means = Means, SAMPLE = SAMPLE, DILUTION = DILUTION,
                          d = d, n = n, h = h, nd = nxd, V = V0)
    else
        inputList <- list(Means = Means, SAMPLE = SAMPLE, DILUTION = DILUTION,
                          d = d, n = n, h = h, MeansX = X, 
                          P = P, D = D, L = L, H = H, PL = rbind(P, L))
    SS <- list(SStot = SStot, SSreg = SSreg,  SSres = SSres, SSsmp = SSsmp,
               R = R, Radj = Radj)
    latexSS <- c("SS\\textsubscript{total}",
                 "SS\\textsubscript{regression}",
                 "SS\\textsubscript{residuals}",
                 "SS\\textsubscript{preparation}", 
                 "R", "R\\textsubscript{adj}")
    mostattributes(SS) <- list(names = names(SS), namesLaTeX = latexSS)
    reg <- list(DFres = DFres, b = b, s2 = s2, qt = qt, C = C, V = V)
    latexReg <- c("DF\\textsubscript{residuals}", "b ($\\hat{\\beta}$)",
                  "s\\textsuperscript{2}", "t(DF, 5 \\%)", "C", "V")
    mostattributes(reg) <- list(names = names(reg), namesLaTeX = latexReg)
    return(list(input = inputList, SS = SS, reg = reg, K = K, KP = KP))
}
