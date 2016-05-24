if(getRversion() >= '2.15.1') utils::globalVariables(c("ok", "dose"))
drfit <- function(data, startlogED50 = NA, chooseone=TRUE,
        probit = TRUE, logit = FALSE, weibull = FALSE,
        linlogit = FALSE, level = 0.95,
        linlogitWrong = NA, allWrong = NA,
        ps0 = 1, ls0 = 0.5, ws0 = 0.5,
        b0 = 2, f0 = 0,
        showED50 = FALSE,
        EDx = NULL, EDx.tolerance = 1e-4)
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
            if (is.na(startlogED50[i])){
                w <- 1/abs(tmp$response - 0.3)
                startlogED50[[i]] <- sum(w * log10(tmp$dose))/sum(w)
            }
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
                    if (linlogit && 
                        length(subset(linlogitWrong,linlogitWrong == i))==0 &&
                        length(subset(allWrong,allWrong == i))==0) {
                        m <- try(nls(response ~ linlogitf(dose,1,f,logED50,b),
                                data=tmp, algorithm="port",
                                start=list(f=f0,logED50=startlogED50[[i]],b=b0)))
                        if (!inherits(m, "try-error")) {
                            fit <- TRUE
                            ri <- ri + 1
                            models[[ri]] <- m
                            s <- summary(m)
                            sigma[[ri]] <- s$sigma
                            rsubstance[[ri]] <- i
                            rndl[[ri]] <- ndl
                            rn[[ri]] <- n
                            runit[[ri]] <- unit
                            rlld[[ri]] <- log10(lowestdose)
                            rlhd[[ri]] <- log10(highestdose)
                            logED50[[ri]] <- coef(m)[["logED50"]]
                            if (logED50[[ri]] > rlhd[[ri]]) {
                                mtype[[ri]] <- "no fit"
                                logED50[[ri]] <- NA
                                logED50low[[ri]] <- NA
                                logED50high[[ri]] <- NA
                                a[[ri]] <- NA
                                b[[ri]] <- NA
                                c[[ri]] <- NA
                            } else {
                                mtype[[ri]] <- "linlogit"
                                logED50conf <- try(confint(m,"logED50",level=level))
                                if (!inherits(logED50conf, "try-error")) {
                                    logED50low[[ri]] <- logED50conf[[1]]
                                    logED50high[[ri]] <- logED50conf[[2]]
                                } else {
                                    logED50low[[ri]] <- NA
                                    logED50high[[ri]] <- NA
                                }
                                a[[ri]] <- coef(m)[["logED50"]]
                                b[[ri]] <- coef(m)[["b"]]
                                c[[ri]] <- coef(m)[["f"]]
                            }
                        }
                    }

                    if (probit &&
                        length(subset(allWrong,allWrong == i))==0) {
                        m <- try(nls(response ~ pnorm(-log10(dose),-logED50,scale),
                                    data=tmp, algorithm="port",
                                    start=list(logED50=startlogED50[[i]],scale=ps0)))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                fit <- TRUE
                                ri <- ri + 1
                                models[[ri]] <- m
                                s <- summary(m)
                                sigma[[ri]] <- s$sigma
                                rsubstance[[ri]] <- i
                                rndl[[ri]] <- ndl
                                rn[[ri]] <- n
                                runit[[ri]] <- unit
                                rlld[[ri]] <- log10(lowestdose)
                                rlhd[[ri]] <- log10(highestdose)
                                logED50[[ri]] <- coef(m)[["logED50"]]
                                c[[ri]] <- NA
                                if (logED50[[ri]] > rlhd[[ri]]) {
                                    mtype[[ri]] <- "no fit"
                                    logED50[[ri]] <- NA
                                    logED50low[[ri]] <- NA
                                    logED50high[[ri]] <- NA
                                    a[[ri]] <- NA
                                    b[[ri]] <- NA
                                } else {
                                    mtype[[ri]] <- "probit"
                                    logED50conf <- try(confint(m,"logED50",level=level))
                                    if (!inherits(logED50conf, "try-error")) {
                                        logED50low[[ri]] <- logED50conf[[1]]
                                        logED50high[[ri]] <- logED50conf[[2]]
                                    } else {
                                        logED50low[[ri]] <- NA
                                        logED50high[[ri]] <- NA
                                    }
                                    a[[ri]] <- coef(m)[["logED50"]]
                                    b[[ri]] <- coef(m)[["scale"]]
                                }
                            }
                        }
                    }

                    if (logit &&
                        length(subset(allWrong,allWrong == i))==0) {
                        m <- try(nls(response ~ plogis(-log10(dose),-logED50,scale),
                                data=tmp, algorithm="port",
                                start=list(logED50=startlogED50[[i]],scale=ls0)))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                fit <- TRUE
                                ri <- ri + 1
                                models[[ri]] <- m
                                s <- summary(m)
                                sigma[[ri]] <- s$sigma
                                rsubstance[[ri]] <- i
                                rndl[[ri]] <- ndl
                                rn[[ri]] <- n
                                runit[[ri]] <- unit
                                rlld[[ri]] <- log10(lowestdose)
                                rlhd[[ri]] <- log10(highestdose)
                                logED50[[ri]] <- a[[ri]] <- coef(m)[["logED50"]]
                                b[[ri]] <- coef(m)[["scale"]]
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
                                    logED50conf <- try(confint(m,"logED50",level=level))
                                    if (!inherits(logED50conf, "try-error")) {
                                        logED50low[[ri]] <- logED50conf[[1]]
                                        logED50high[[ri]] <- logED50conf[[2]]
                                    } else {
                                        logED50low[[ri]] <- NA
                                        logED50high[[ri]] <- NA
                                    }
                                }
                            }
                        }
                    }

                    if (weibull &&
                        length(subset(allWrong,allWrong == i))==0) {
                        m <- try(nls(response ~ pweibull(-log10(dose)+location,shape),
                                data=tmp, algorithm="port",
                                start=list(location=startlogED50[[i]],shape=ws0)))
                        if (chooseone==FALSE || fit==FALSE) {
                            if (!inherits(m, "try-error")) {
                                ri <- ri + 1
                                models[[ri]] <- m
                                a[[ri]] <- coef(m)[["location"]]
                                b[[ri]] <- coef(m)[["shape"]]
                                sqrdev <- function(logdose) {
                                    (0.5 - pweibull( - logdose + a[[ri]], b[[ri]]))^2 
                                }
                                logED50[[ri]] <- nlm(sqrdev,startlogED50[[i]])$estimate
                                if (sqrdev(logED50[[ri]]) > 0.1) {
                                    message("\nCan't find ED50 for fitted weibull model of ",i,
                                        " data\nwith startlogED50 ", startlogED50[[i]],"\n")
                                    ri <- ri - 1
                                    length(a) <- length(b) <- ri
                                    length(logED50) <- ri
                                } else {
                                    c[[ri]] <- NA
                                    fit <- TRUE
                                    s <- summary(m)
                                    sigma[[ri]] <- s$sigma
                                    rsubstance[[ri]] <- i
                                    rndl[[ri]] <- ndl
                                    rn[[ri]] <- n
                                    runit[[ri]] <- unit
                                    rlld[[ri]] <- log10(lowestdose)
                                    rlhd[[ri]] <- log10(highestdose)
                                    logED50low[[ri]] <- NA
                                    logED50high[[ri]] <- NA
                                    if (logED50[[ri]] > rlhd[[ri]]) {
                                        mtype[[ri]] <- "no fit"
                                        logED50[[ri]] <- NA
                                        a[[ri]] <- NA
                                        b[[ri]] <- NA
                                    } else {
                                        mtype[[ri]] <- "weibull"
                                    }
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
            if (mtype[[row.i]] %in% c("probit", "logit", "weibull", "linlogit")) {
                for (ED in EDx) {
                    of <- function(x) {
                        abs(predict(models[[row.i]], data.frame(dose = 10^x)) - 
                            (1 - (ED/100)))
                    } 
                    # Search over interval starting an order of magnitude below 
                    # the lowest dose up to one order of magnitude above the
                    # highest dose
                    o = optimize(of, 
                                 results[row.i, c("lld", "lhd")] + c(-1, 1))
                    # Only keep results within the tolerance
                    if ((o$objective) < EDx.tolerance) {
                        logdose.ED = o$minimum
                        results[row.i, paste0("EDx", ED)] <- 10^logdose.ED 
                    }
                }
            }
        }
    }

    rownames(results) <- 1:ri
    attr(results, "models") <- models
    return(results)
}
# vim: set ts=4 sw=4 expandtab:
