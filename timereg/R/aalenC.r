aalenBaseC <- function(times, fdata, designX, status, id, clusters, robust = 0, 
    sim = 0, retur = 0, antsim = 1000, weighted.test = 1, covariance = 0, 
    resample.iid = 0, namesX = NULL, silent = 0, scale = 1) 
{ ## {{{
    Ntimes <- length(times)
    designX <- as.matrix(designX)
    if (is.matrix(designX) == TRUE) 
        p <- as.integer(dim(designX)[2])
    if (is.matrix(designX) == TRUE) 
        nx <- as.integer(dim(designX)[1])
    if (robust == 0 & sim >= 1) 
        robust <- 1
    devi <- rep(0, 1)
    cumint <- matrix(0, Ntimes, p + 1)
    Vcumint <- cumint
    if (retur == 1) 
        cumAi <- matrix(0, Ntimes, fdata$antpers)
    else cumAi <- 0
    test <- matrix(0, antsim, 3 * p)
    testOBS <- rep(0, 3 * p)
    testval <- c()
    unifCI <- c()
    rani <- -round(runif(1) * 10000)
    if (sim >= 1) 
        simUt <- matrix(0, Ntimes, 50 * p)
    else simUt <- NULL
    Ut <- matrix(0, Ntimes, p + 1)
    if (covariance == 1) 
        covs <- matrix(0, Ntimes, p * p)
    else covs <- 0
    if (resample.iid == 1) {
        B.iid <- matrix(0, Ntimes, fdata$antclust * p)
    }
    else B.iid <- NULL
    if (robust == 2) {
        aalenout <- .C("aalen", as.double(times), as.integer(Ntimes), 
            as.double(designX), as.integer(nx), as.integer(p), 
            as.integer(fdata$antpers), as.double(fdata$start), 
            as.double(fdata$stop), as.double(cumint), as.double(Vcumint), 
            as.integer(status), PACKAGE = "timereg")
        robV <- NULL
        cumAI <- NULL
        test <- NULL
    }
    else {
        robVar <- Vcumint
        aalenout <- .C("robaalenC", as.double(times), as.integer(Ntimes), 
            as.double(designX), as.integer(nx), as.integer(p), 
            as.integer(fdata$antpers), as.double(fdata$start), 
            as.double(fdata$stop), as.double(cumint), as.double(Vcumint), 
            as.double(robVar), as.integer(sim), as.integer(antsim), 
            as.integer(retur), as.double(cumAi), as.double(test), 
            as.integer(rani), as.double(testOBS), as.integer(status), 
            as.double(Ut), as.double(simUt), as.integer(id), 
            as.integer(weighted.test), as.integer(robust), as.integer(covariance), 
            as.double(covs), as.integer(resample.iid), as.double(B.iid), 
            as.integer(clusters), as.integer(fdata$antclust), 
            as.double(devi), as.integer(silent), PACKAGE = "timereg")
        if (covariance == 1) {
            covit <- matrix(aalenout[[26]], Ntimes, p * p)
            cov.list <- list()
            for (i in 1:Ntimes) cov.list[[i]] <- matrix(covit[i, 
                ], p, p)
        }
        else cov.list <- NULL
        if (resample.iid == 1) {
            covit <- matrix(aalenout[[28]], Ntimes, fdata$antclust * 
                p)
            B.iid <- list()
            for (i in (0:(fdata$antclust - 1)) * p) {
                B.iid[[i/p + 1]] <- as.matrix(covit[, i + (1:p)])
                colnames(B.iid[[i/p + 1]]) <- namesX
            }
        }
        robV <- matrix(aalenout[[11]], Ntimes, p + 1)
        if (retur == 1) {
            cumAi <- matrix(aalenout[[15]], Ntimes, fdata$antpers * 
                1)
            cumAi <- list(time = times, dM = cumAi, dM.iid = cumAi)
        }
        else cumAi <- NULL
    }
    if (sim >= 1) {
        Uit <- matrix(aalenout[[21]], Ntimes, 50 * p)
        UIt <- list()
        for (i in (0:49) * p) UIt[[i/p + 1]] <- as.matrix(Uit[, 
            i + (1:p)])
        Ut <- matrix(aalenout[[20]], Ntimes, (p + 1))
        test <- matrix(aalenout[[16]], antsim, 3 * p)
        testOBS <- aalenout[[18]]
        for (i in 1:(3 * p)) testval <- c(testval, pval(test[, 
            i], testOBS[i]))
        for (i in 1:p) unifCI <- as.vector(c(unifCI, percen(test[, 
            i], 0.95)))
        pval.testBeq0 <- as.vector(testval[1:p])
        pval.testBeqC <- as.vector(testval[(p + 1):(2 * p)])
        pval.testBeqC.is <- as.vector(testval[(2 * p + 1):(3 * 
            p)])
        obs.testBeq0 <- as.vector(testOBS[1:p])
        obs.testBeqC <- as.vector(testOBS[(p + 1):(2 * p)])
        obs.testBeqC.is <- as.vector(testOBS[(2 * p + 1):(3 * 
            p)])
        sim.testBeq0 <- as.matrix(test[, 1:p])
        sim.testBeqC <- as.matrix(test[, (p + 1):(2 * p)])
        sim.testBeqC.is <- as.matrix(test[, (2 * p + 1):(3 * 
            p)])
    }
    else {
        test <- NULL
        unifCI <- NULL
        Ut <- NULL
        UIt <- NULL
        pval.testBeq0 <- NULL
        pval.testBeqC <- NULL
        obs.testBeq0 <- NULL
        obs.testBeqC <- NULL
        sim.testBeq0 <- NULL
        sim.testBeqC <- NULL
        sim.testBeqC.is <- NULL
        pval.testBeqC.is <- NULL
        obs.testBeqC.is <- NULL
    }
    cumint <- matrix(aalenout[[9]], Ntimes, p + 1)
    Vcumint <- matrix(aalenout[[10]], Ntimes, p + 1)
    devi <- aalenout[[31]]
    list(cum = cumint, var.cum = Vcumint, robvar.cum = robV, 
        residuals = cumAi, pval.testBeq0 = pval.testBeq0, obs.testBeq0 = obs.testBeq0, 
        pval.testBeqC = pval.testBeqC, pval.testBeqC.is = pval.testBeqC.is, 
        obs.testBeqC = obs.testBeqC, obs.testBeqC.is = obs.testBeqC.is, 
        sim.testBeq0 = sim.testBeq0, sim.testBeqC = sim.testBeqC, 
        sim.testBeqC.is = sim.testBeqC.is, conf.band = unifCI, 
        test.procBeqC = Ut, sim.test.procBeqC = UIt, covariance = cov.list, 
        B.iid = B.iid, deviance = devi)
} ## }}}

