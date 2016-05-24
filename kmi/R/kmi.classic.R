### function to get the censoring times and distribution
### from competing risks right-censored data

kmi.classic <- function(y, x, etype, failcode, epsilon,
                        bootstrap, nboot) {
    if (!is.Surv(y)) stop("y must be a Surv object")
    if (attr(y, "type") != "right") stop("Can only handle right censored data")
    if (is.null(etype)) stop("Argument 'etype' is missing with no default")
    ## Depending on how the model is specified (see example(survfit)),
    ## etype might not be 0 when the observation is censored
    ## etype[y[, 2] == 0] <- 0
    cens.times <- sort(unique(y[, 1][y[, 2] == 0]))
    ind <- which(y[, 2] == 0 | etype == failcode)
    ## itimes are the time that need imputation
    ## otimes don't need imputation
    itimes <- y[-ind, 1]
    otimes <- y[ind, 1]
    cn <- colnames(x)

    xx <- x[-ind, , drop = FALSE]
    ## let's deal with missing values in a really ugly way,
    ## i.e., mean imputation
    if (!is.null(cn)) {
        if (!all(!is.na(xx))) {
            warning("Missing values in the variable(s) used for modelling the censoring distribution.\nMean imputation used")
            mm <- apply(x, 2, mean, na.rm = TRUE)
            for (i in seq_len(ncol(xx))) {
                xx[which(is.na(xx[, i])), i] <- mm[i]
            }
        }
    }
    
    if (bootstrap) { # simple bootstrap with remplacement here
        index <- lapply(seq_len(nboot), function(k) {
            sample(seq_len(nrow(y)), nrow(y),
                   replace = TRUE)
        })## might save some time to compute the index within the loop

        ff <- formula(Surv(y[index[[l]], 1],
                           y[index[[l]], 2] == 0) ~ 1)
        
        if (!is.null(cn)) {
            g <- array(0, dim = c(length(cens.times), length(itimes), nboot))
            ff <- update.formula(ff, paste(". ~", paste(cn, collapse = "+")))
            for (l in seq_len(nboot)) {
                temp <- coxph(ff, as.data.frame(x))
                tmp <- summary(survfit(temp, as.data.frame(xx)))
                ordre <- findInterval(cens.times, tmp$time)
                ordre[ordre == 0] <- NA
                g[,, l] <- tmp$surv[ordre, ]
                g[,, l][is.na(g[, , l])] <- 1
            }
            g <- apply(g, c(1, 2), mean)
            gg <- rbind(1, g)
        } else {
            g <-  matrix(0, nrow = nboot, ncol = length(cens.times))
            for (l in seq_len(nboot)) {
                tmp <- summary(survfit(Surv(y[index[[l]], 1],
                                            y[index[[l]], 2] == 0) ~ 1))
                ordre <- findInterval(cens.times, tmp$time)
                ordre[ordre == 0] <- NA
                g[l, ] <- tmp$surv[ordre]
                g[l, ][is.na(g[l, ])] <- 1
            }
            g <- apply(g, 2, mean)
            gg <- matrix(rep(c(1, g), length(itimes)),
                         nrow = length(g) + 1)
        }
        
    } else {
            
        ff <- formula(Surv(y[, 1], y[, 2] == 0) ~ 1)
        if (!is.null(cn)) {

            ff <- update.formula(ff, paste(". ~", paste(cn, collapse = "+")))
            temp <- coxph(ff, as.data.frame(x))
            g <- summary(survfit(temp, as.data.frame(xx)),
                         times = cens.times, extend = TRUE)$surv
            gg <- rbind(1, g)
        } else {
            g <- summary(survfit(Surv(y[, 1], y[, 2] == 0) ~ 1))$surv
            gg <- matrix(rep(c(1, g), length(itimes)),
                         nrow = length(g) + 1)
        }
    }
    
    a <- FALSE
    if (y[, 2][which.max(y[, 1])] != 0) { # will be true if the last time is an event
        cens.times <- c(cens.times, max(y[, 1]) + epsilon)
        a <- TRUE
    }
    
    list(gg = gg, cens.times = cens.times, itimes = itimes,
         otimes = otimes, place = ind, a = a)
}
