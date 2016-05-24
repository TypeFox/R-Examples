cpf <- function(formula, data, subset, na.action, conf.int = 0.95, failcode) {
    
    call <- match.call()
    if ((mode(call[[2]]) == "call" && call[[2]][[1]] == as.name("Hist")) ||
        inherits(formula, "Hist")) {
        formula <- eval(parse(text = paste(deparse(call[[2]]), 1, sep="~")))
        environment(formula) <- parent.frame()
    }

    m <- match.call(expand.dots = FALSE)
    m$formula <- formula
    m$failcode <- NULL
    m[[1]] <- as.name("prodlim")
    fit <- eval(m, parent.frame())
    if (attr(fit$model.response, "model") != "competing.risks")
        stop("This is not competing risks data")
    if (attr(fit$model.response, "cens.type") != "rightCensored")
        stop("This function gives conditional probability estimates \nonly for right censored data")
    n <- sum(fit$n.risk[c(1, fit$size.strata[-length(fit$size.strata)] + 1)])
    
### Calcul de CP
    if (missing(failcode)) {
        failcode <- ind <- as.numeric(names(fit$n.event)[1])
    }
    else {
        if (!(failcode %in% names(fit$cuminc))) {
            stop("'failcode' is not amongst the event types observed in the data")
        }
        ind <- which(failcode == names(fit$cuminc))
    }
    if (length(fit$cuminc) > 2) {
        cif.other <- do.call(fit$cuminc[-ind], "+")
        n.event.other <- do.call(fit$n.event[-ind], "+")
##         cif.other <- rowSums(fit$cuminc[, -ind])
##         n.event.other <- rowSums(fit$n.event[, -ind])
    }
    else {
        cif.other <- fit$cuminc[[-ind]]
        n.event.other <- fit$n.event[[-ind]]
    }
    n.event <- data.frame(fit$n.event[[ind]], n.event.other)
    names(n.event) <- c(as.name(failcode), "other")
    cp <- fit$cuminc[[ind]] / (1 - cif.other)
    
### Variance computation
    ## getting S(t-)
    size <- cumsum(fit$size.strata) + 1 
    sminus <- c(1, fit$surv[-length(fit$surv)])
    sminus[size[-length(size)]] <- 1
    if (length(fit$size.strata) == 1) {
        var.cp <- (sminus^2 / (1 - cif.other)^4) *
            cumsum(((1 - cif.other)^2 * fit$n.event[[ind]] + fit$cuminc[[ind]]^2 * n.event.other) /
                   (fit$n.risk * (fit$n.risk - 1)))
    }
    else {
        ns <- c(0, cumsum(fit$size.strata))
        var.cp <- lapply(1:length(fit$size.strata), function(i) {
            temp <- (sminus[(ns[i] + 1):ns[i + 1]]^2 / (1 - cif.other[(ns[i]+1):ns[i+1]])^4) *
                cumsum(((1 - cif.other[(ns[i]+1):ns[i+1]])^2 * fit$n.event[[ind]][(ns[i]+1):ns[i+1]] +
                        fit$cuminc[[ind]][(ns[i]+1):ns[i+1]]^2 * n.event.other[(ns[i]+1):ns[i+1]]) /
                       (fit$n.risk[(ns[i]+1):ns[i+1]] * (fit$n.risk[(ns[i]+1):ns[i+1]] - 1)))
            temp
        })
        var.cp <- do.call(c, var.cp)
    }

    level <- 1 - conf.int
    upper <- cp + qnorm(1 - level/2) * sqrt(var.cp)
    lower <- cp - qnorm(1 - level/2) * sqrt(var.cp)

    if (length(fit$size.strata) == 2) {
        max.times <- min(fit$time[cumsum(fit$size.strata)])
        times <- sort(unique(fit$time[fit$time <= max.times]))
        ind.G1 <- 1:fit$size.strata[1]
        ind.G2 <- (fit$size.strata[1] + 1):sum(fit$size.strata)
        ind1 <- 1:length(times)
        ind2 <- (length(times) + 1):(2 * length(times))
        
        ## censoring distribution
        n.group <- fit$n.risk[c(1, fit$size.strata[1] + 1)]
        cens.distr <- c(cumprod((fit$n.risk[ind.G1] - fit$n.lost[ind.G1]) / fit$n.risk[ind.G1]),
                        cumprod((fit$n.risk[ind.G2] - fit$n.lost[ind.G2]) / fit$n.risk[ind.G2]))
        cens.distr <- c(1, cens.distr[-length(cens.distr)])
        cens.distr[size[-length(size)]] <- 1
        ind.time1 <- findInterval(times, fit$time[ind.G1])
        ind.time2 <- findInterval(times, fit$time[ind.G2])
        t.cp <- c(ifelse(ind.time1 == 0, 0, cp[ind.G1][ind.time1]),
                  ifelse(ind.time2 == 0, 0, cp[ind.G2][ind.time2]))
        t.cens.distr <- c(ifelse(ind.time1 == 0, 0, cens.distr[ind.G1][ind.time1]),
                          ifelse(ind.time2 == 0, 0, cens.distr[ind.G2][ind.time2]))
        
        ## test
        weight <- t.cens.distr[ind1] * t.cens.distr[ind2] /
            ((n.group[1] / n) * t.cens.distr[ind1] + (n.group[2] / n) + t.cens.distr[ind2])
        wcp <- sqrt(n.group[1] * n.group[2] / (n.group[1] + n.group[2])) *
            sum(weight * (t.cp[ind1] - t.cp[ind2]))
        
        ## test's variance
        t.surv <- c(ifelse(ind.time1 == 0, 1, fit$surv[ind.G1][ind.time1]),
                    ifelse(ind.time2 == 0, 1, fit$surv[ind.G2][ind.time2]))
        t.cif <- c(ifelse(ind.time1 == 0, 0, fit$cuminc[[ind]][ind.G1][ind.time1]),
                   ifelse(ind.time2 == 0, 0, fit$cuminc[[ind]][ind.G2][ind.time2]))
        t.cif.other <- c(ifelse(ind.time1 == 0, 0, cif.other[ind.G1][ind.time1]),
                         ifelse(ind.time2 == 0, 0, cif.other[ind.G2][ind.time2]))
        t.n.event <- c(ifelse(ind.time1 == 0, 0, fit$n.event[[ind]][ind.G1][ind.time1]),
                       ifelse(ind.time2 == 0, 0, fit$n.event[[ind]][ind.G2][ind.time2]))
        t.n.event.other <- c(ifelse(ind.time1 == 0, 0, n.event.other[ind.G1][ind.time1]),
                             ifelse(ind.time2 == 0, 0, n.event.other[ind.G2][ind.time2]))
        t.n.risk <- c(ifelse(ind.time1 == 0, n.group[1], fit$n.risk[ind.G1][ind.time1]),
                      ifelse(ind.time2 == 0, n.group[2], fit$n.risk[ind.G2][ind.time2]))

        int1 <- intt(weight, t.surv[ind1], t.cif.other[ind1], times)
        int2 <- intt(weight, t.surv[ind2], t.cif.other[ind2], times)
        
        sigma.wcp <- c(sum(int1^2 * ((1 - t.cif.other[ind1])^2 * t.n.event[ind1] +
                                 t.cif[ind1]^2 * t.n.event.other[ind1]) /
                       (t.n.risk[ind1] * (t.n.risk[ind1] - 1)), na.rm = TRUE),
                       sum(int2^2 * ((1 - t.cif.other[ind2])^2 * t.n.event[ind2] +
                                 t.cif[ind2]^2 * t.n.event.other[ind2]) /
                       (t.n.risk[ind2] * (t.n.risk[ind2] - 1)), na.rm = TRUE))
        var.wcp <- (n.group[2] / n) * n*sigma.wcp[1] + (n.group[1] / n) * n*sigma.wcp[2]
        
        z <- wcp / sqrt(var.wcp)
        pvalue <- 2 * pnorm(-abs(z))

        strata <- paste(colnames(fit$X), fit$X[, ], sep = "=")

        zzz <- list(cp = cp, var = var.cp, time = fit$time, lower = lower,
                    upper = upper, n.risk = fit$n.risk,
                    n.event = n.event, n.lost = fit$n.lost, size.strata = fit$size.strata,
                    X = fit$X, strata = strata, call = call, z = z, p = pvalue, failcode = failcode)
        class(zzz) <- "cpf"
    }
    
    else {
        if (length(fit$size.strata) > 2) {
            warning("The test is only available for comparing 2 samples")
        }

        strata <- paste(colnames(fit$X), fit$X[, ], sep = "=")
        
        zzz <- list(cp = cp, var = var.cp, time = fit$time, lower = lower,
                    upper = upper, n.risk = fit$n.risk,
                    n.event = n.event, n.lost = fit$n.lost, size.strata = fit$size.strata,
                    X = fit$X, strata = strata, call = call, z = NULL, p = NULL, failcode = failcode)
        class(zzz) <- "cpf"
    }
    
    return(zzz)
    
}

intt <- function(weight, surv, fOther, times) {
    lt <- length(times)
    res <- sapply(1:lt, function(i) {
        sum((weight[i:lt] * surv[i:lt]) / (1 - fOther[i:lt])^2, na.rm = TRUE)
     })
    res
}


### Extract method for cpf objects
"[.cpf" <- function(x, ..., drop = FALSE) {

    if (length(x$strata) == 0) {
        warning("'cpf' has only 1 conditional probability curve")
        return(x)
    }
    
    if (missing(..1)) i <- NULL else i <- sort(..1)
    
    if (is.null(i)) ind <- rep(TRUE, length(x$time))
    else {
        if (is.character(i)) {
            where <- match(i, x$strata)
            if (any(is.na(where)))
                stop(paste("subscript(s)",
                           paste(i[is.na(where)], sep = " "),
                           "not matched"))
        }
        else {
            where <- i
            if (max(where) > dim(x$X)[1]) {
                stop(paste("'cpf' object has only",
                           length(x$size.strata),
                           "conditional probability curves", sep = " "))
            }
        }
        
        lstrat <- rep(1:dim(x$X)[1], x$size.strata)
        ind <- lstrat %in% where
        
        if (length(where) <= 1) {
            x$X <- NULL
            x$z <- x$p <- NULL
        }
        else {
            x$X <- x$X[where, , drop = drop]
        }
        
        x$cp <- x$cp[ind]
        x$var <- x$var[ind]
        x$time <- x$time[ind]
        x$lower <- x$lower[ind]
        x$upper <- x$upper[ind]
        x$n.risk <- x$n.risk[ind]
        x$n.event <- x$n.event[ind, , drop = drop]
        x$n.lost <- x$n.lost[ind]
        x$size.strata <- x$size.strata[where]
        x$strata <- x$strata[where]
    }
    
    x
}
