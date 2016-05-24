drcfit <- function(data, chooseone=TRUE,
        probit = TRUE, logit = FALSE, weibull = FALSE,
        linlogit = FALSE, level = 0.95,
        showED50 = FALSE,
        EDx = NULL)
{
    if(!is.null(data$ok)) data <- subset(data,ok!="no fit") # Don't use data
                                                            # with ok set to
                                                            # "no fit"
    substances <- levels(data$substance)

    ri <- rix <- 0                  # ri is the index over the result rows
                                    # rix is used later to check if any
                                    # model result was appended
    rsubstance <- array()           # the substance names in the results
    rndl <- vector()                # number of dose levels
    rn <- vector()                  # mean number of replicates 
                                    # in each dose level
    runit <- vector()               # vector of units for each result row
    rlhd <- rlld <- vector()        # highest and lowest doses tested
    mtype <- array()                # the modeltypes
    sigma <- array()                # the standard deviation of the residuals
    logED50 <- vector()
    logED50low <- logED50high <- vector()
    a <- b <- c <- vector()

    models <- list()                  # a list containing the dose-response models

    splitted <- split(data,data$substance)
    for (i in substances) {
        tmp <- splitted[[i]]
        fit <- FALSE
        if (length(tmp) != 0) {
            unit <- levels(as.factor(as.vector(tmp$unit)))
            message("\n",i,": Fitting data...\n")
        } else {
            unit <- ""
            message("\n",i,": No data\n")
        }
        if (length(unit) == 0) {
            unit <- NA
        }
        if (length(unit) > 1) {
            message("More than one unit for substance ",i,", halting\n\n")
            break
        }
        if (length(tmp$response) == 0) {
            nodata = TRUE
        } else {
            nodata = FALSE
        }
        rix <- ri
        if (nodata) {
            n <- ndl <- 0
        } else {
            ndl <- length(levels(factor(tmp$dose)))
            n <- length(tmp$response)
            highestdose <- max(tmp$dose)
            lowestdose <- min(tmp$dose)
            lhd <- log10(highestdose)
            lld <- log10(lowestdose)
            responseathighestdose <- mean(subset(tmp,dose==highestdose)$response)
            responseatlowestdose <- mean(subset(tmp,dose==lowestdose)$response)
            if (responseathighestdose < 0.5) {
                inactive <- FALSE
                if (responseatlowestdose < 0.5) {
                    active <- TRUE
                } else {
                    active <- FALSE
                    if (linlogit)
                    {
                        m <- try(drm(response ~ dose, data = tmp, fct = BC.4(fixed = c(NA, 1, NA, NA))))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                fit <- TRUE
                                ri <- ri + 1
                                mtype[[ri]] <- "linlogit"
                                models[[ri]] <- m
                                s <- summary(m)
                                sigma[[ri]] <- s$rseMat[1, "rse"]
                                rsubstance[[ri]] <- i
                                rndl[[ri]] <- ndl
                                rn[[ri]] <- n
                                runit[[ri]] <- unit
                                rlld[[ri]] <- log10(lowestdose)
                                rlhd[[ri]] <- log10(highestdose)
                                logED50[[ri]] <- NA
                                logED50low[[ri]] <- NA
                                logED50high[[ri]] <- NA
                                a[[ri]] <- coef(m)[[2]]
                                b[[ri]] <- coef(m)[[1]]
                                c[[ri]] <- coef(m)[[3]]
                                ED50 <- try(ED(m, 50, interval = "delta", 
                                               lower = lowestdose / 10,
                                               upper = highestdose * 10,
                                               display = FALSE))
                                if (!inherits(ED50, "try-error")) {
                                    logED50[[ri]] <- log10(ED50["1:50", "Estimate"])
                                    logED50low[[ri]] <- log10(ED50["1:50", "Lower"])
                                    logED50high[[ri]] <- log10(ED50["1:50", "Upper"])
                                    if (logED50[[ri]] > rlhd[[ri]]) {
                                        mtype[[ri]] <- "no fit"
                                    }
                                } 
                            }
                        }
                    }
                    if (probit)
                    {
                        m <- try(drm(response ~ dose, data = tmp, 
                                     fct = LN.2()))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                fit <- TRUE
                                ri <- ri + 1
                                models[[ri]] <- m
                                s <- summary(m)
                                sigma[[ri]] <- s$rseMat[1, "rse"]
                                rsubstance[[ri]] <- i
                                rndl[[ri]] <- ndl
                                rn[[ri]] <- n
                                runit[[ri]] <- unit
                                rlld[[ri]] <- log10(lowestdose)
                                rlhd[[ri]] <- log10(highestdose)
                                logED50[[ri]] <- log10(coef(m)[[2]])
                                a[[ri]] <- coef(m)[[2]]
                                b[[ri]] <- coef(m)[[1]]
                                c[[ri]] <- NA
                                if (logED50[[ri]] > rlhd[[ri]]) {
                                    mtype[[ri]] <- "no fit"
                                    logED50[[ri]] <- NA
                                    logED50low[[ri]] <- NA
                                    logED50high[[ri]] <- NA
                                } else {
                                    mtype[[ri]] <- "probit"
                                    ED50 <- ED(m, 50, interval = "delta", display = FALSE)
                                    logED50low[[ri]] <- log10(ED50["1:50", "Lower"])
                                    logED50high[[ri]] <- log10(ED50["1:50", "Upper"])
                                }
                            }
                        }
                    }
                    if (logit)
                    {
                        m <- try(drm(response ~ dose, data = tmp, fct = LL.2()))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                fit <- TRUE
                                ri <- ri + 1
                                models[[ri]] <- m
                                s <- summary(m)
                                sigma[[ri]] <- s$rseMat[1, "rse"]
                                rsubstance[[ri]] <- i
                                rndl[[ri]] <- ndl
                                rn[[ri]] <- n
                                runit[[ri]] <- unit
                                rlld[[ri]] <- log10(lowestdose)
                                rlhd[[ri]] <- log10(highestdose)
                                logED50[[ri]] <- log10(coef(m)[[2]])
                                a[[ri]] <- coef(m)[[2]]
                                b[[ri]] <- coef(m)[[1]]
                                c[[ri]] <- NA
                                if (logED50[[ri]] > rlhd[[ri]]) {
                                    mtype[[ri]] <- "no fit"
                                    logED50[[ri]] <- NA
                                    logED50low[[ri]] <- NA
                                    logED50high[[ri]] <- NA
                                    a[[ri]] <- NA
                                    b[[ri]] <- NA
                                } else {
                                    mtype[[ri]] <- "logit"
                                    ED50 <- ED(m, 50, interval = "delta", display = FALSE)
                                    logED50low[[ri]] <- log10(ED50["1:50", "Lower"])
                                    logED50high[[ri]] <- log10(ED50["1:50", "Upper"])
                                }
                            }
                        }

                    }
                    if (weibull)
                    {
                        m <- try(drm(response ~ dose, data = tmp, fct = W1.2()))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                fit <- TRUE
                                ri <- ri + 1
                                models[[ri]] <- m
                                s <- summary(m)
                                sigma[[ri]] <- s$rseMat[1, "rse"]
                                rsubstance[[ri]] <- i
                                rndl[[ri]] <- ndl
                                rn[[ri]] <- n
                                runit[[ri]] <- unit
                                rlld[[ri]] <- log10(lowestdose)
                                rlhd[[ri]] <- log10(highestdose)
                                c[[ri]] <- NA
                                ED50 <- ED(m, 50, interval = "delta", display = FALSE)
                                logED50[[ri]] <- log10(ED50["1:50", "Estimate"])
                                if (logED50[[ri]] > rlhd[[ri]]) {
                                    mtype[[ri]] <- "no fit"
                                    logED50[[ri]] <- NA
                                    logED50low[[ri]] <- NA
                                    logED50high[[ri]] <- NA
                                    a[[ri]] <- NA
                                    b[[ri]] <- NA
                                } else {
                                    mtype[[ri]] <- "weibull"
                                    logED50low[[ri]] <- log10(ED50["1:50", "Lower"])
                                    logED50high[[ri]] <- log10(ED50["1:50", "Upper"])
                                    a[[ri]] <- logED50[[ri]]
                                    b[[ri]] <- coef(m)[[1]]
                                }
                            }
                        }

                    }
                }
            } else {
                inactive <- TRUE
            }
        }
        if (ri == rix) {          # if no entry was appended for this substance
            ri <- ri + 1
            rsubstance[[ri]] <- i
            rndl[[ri]] <- ndl
            rn[[ri]] <- n
            if (nodata) {
                rlld[[ri]] <- rlhd[[i]] <- NA
                mtype[[ri]] <- "no data"
                runit[[ri]] <- NA
            } else {
                rlld[[ri]] <- log10(lowestdose)
                rlhd[[i]] <- log10(highestdose)
                runit[[ri]] <- unit
                if (inactive) {
                    mtype[[ri]] <- "inactive"
                } else {
                    if (active) {
                        mtype[[ri]] <- "active"
                    } else {
                        mtype[[ri]] <- "no fit"
                    }
                }
            }
            sigma[[ri]] <- NA
            logED50[[ri]] <- NA
            logED50low[[ri]] <- NA
            logED50high[[ri]] <- NA
            a[[ri]] <- NA
            b[[ri]] <- NA
            c[[ri]] <- NA
        }
    }

    results <- data.frame(rsubstance, rndl, rn, rlld, rlhd, mtype, 
        logED50, logED50low, logED50high, runit, sigma, a, b)
    lower_level_percent = paste(100 * (1 - level)/2, "%", sep = "")
    upper_level_percent = paste(100 * (1 + level)/2, "%", sep = "")
    names(results) <- c("Substance","ndl","n","lld","lhd","mtype","logED50",
        lower_level_percent, upper_level_percent,
        "unit","sigma","a","b")

    if (linlogit) {
        results$c <- c
    }

    if (showED50) {
        results[c("ED50", paste("ED50", c(lower_level_percent, upper_level_percent)))] <-
          10^results[7:9]
    }

    if (!is.null(EDx)) {
        for (row.i in 1:ri) {
            m <- models[[row.i]]
            mtype <- as.character(results[row.i, "mtype"])
            if (mtype %in% c("probit", "logit", "weibull", "linlogit")) {
                for (EDi in EDx) {
                    EDx.drc = try(ED(m, EDi, interval = "delta", display = FALSE, level = level))
                    if (!inherits(EDx.drc, "try-error")) {
                        results[row.i, paste0("EDx", EDi)] <- EDx.drc[paste0("1:", EDi), "Estimate"]
                        results[row.i, paste0("EDx", EDi, " ", lower_level_percent)] <- EDx.drc[paste0("1:", EDi), 
                                                                                                "Lower"]
                        results[row.i, paste0("EDx", EDi, " ", upper_level_percent)] <- EDx.drc[paste0("1:", EDi), 
                                                                                                "Upper"]
                    }
                }
            }
        }
    }

    attr(results, "models") <- models
    return(results)
}
# vim: set ts=4 sw=4 expandtab:
