################################################################################
## package 'secr'
## pmixProfileLL.R
## 2015-11-17
## 2015-12-30 NOTE: better to generalize this to any real parameter
################################################################################

## fit CL model with range of fixed beta values for mixing proportion
## ... arguments passed to secr.fit
pmixProfileLL <- function (CH, model = list(g0~h2, sigma~h2), CL = TRUE,
                           pmvals = seq(0.01, 0.99, 0.01), pmi = 5, ncores = 1, ...) {
    oneLL <- function (pm, args, pmi) {
        ## CL = TRUE, model = list(g0~h2, sigma~h2)
        args$details$fixedbeta[pmi] <- logit(pm)
        do.call(secr.fit, args)
    }

    args <- list(...)
    args$capthist <- CH
    if (is.null(args$details))
        args$details <- vector('list')
    if (!is.null(args$details$fixedbeta))
        warning("overriding fixedbeta in input")

# Alternatively, get pmi automatically NOT READY
#     if ('formula' %in% class(model)) model <- list(model)
#     model <- stdform (model)  ## named, no LHS
#     pnames <- pnames[!(pnames %in% fnames)]   ## drop fixed real parameters
#     model <- defaultmodel[pnames]             ## select real parameters
#     valid.model(model, CL, detectfn, hcov, details$userdist, names(sessioncov))
#
#     design <- secr.design.MS(CH, models = stdform(model), timecov = args$timecov, sessioncov =
#                                  args$sessioncov, groups = args$groups, hcov = args$hcov)
#     np <- sapply(design$designMatrices, ncol)
#     NP <- sum(np)
#     parindx <- split(1:NP, rep(1:length(np), np))
#     names(parindx) <- names(np)[np>0]
#     pmi <- parindx$pmix[1]

    ## pmi = index of pmix in beta vector
    # args$CL <- TRUE
    # model$D <- NULL
    # model$noneuc <- NULL

    args$model <- model
    args$CL <- CL
    args$details$fixedbeta <- rep(NA, pmi)  #! assume last

    npm <- length(pmvals)
    if (ncores == 1)
        outCL <- lapply(pmvals, oneLL, args, pmi)
    else {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        outCL <- parLapply(clust, pmvals, oneLL, args, pmi)
        stopCluster(clust)
    }
    sapply(outCL, logLik)
}

