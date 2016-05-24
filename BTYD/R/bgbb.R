################################################################################ Beta-Geometric Beta-Binomial Functions

library(hypergeo)

bgbb.rf.matrix.LL <- function(params, rf.matrix) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.rf.matrix.LL")
    
    tryCatch(x.vec <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x.vec <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.periods.vec <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(n.custs.vec <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.rf.matrix.LL: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    LL.sum <- sum(n.custs.vec * bgbb.LL(params, x.vec, t.x.vec, n.periods.vec))
    
    return(LL.sum)
}

bgbb.LL <- function(params, x, t.x, n.cal) {
    
    max.length <- max(length(x), length(t.x), length(n.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.LL")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    denom.ab <- lbeta(alpha, beta)
    denom.gd <- lbeta(gamma, delta)
    
    indiv.LL.sum <- lbeta(alpha + x, beta + n.cal - x) - denom.ab + lbeta(gamma, 
        delta + n.cal) - denom.gd
    
    check <- n.cal - t.x - 1
    addition <- function(alpha, beta, gamma, delta, denom.ab, denom.gd, x, t.x, check) {
        ii <- 0:check
        log(sum(exp(lbeta(alpha + x, beta + t.x - x + ii) - denom.ab + lbeta(gamma + 
            1, delta + t.x + ii) - denom.gd)))
    }
    # for every element of vectors for which t.x<n.cal, add the result of 'addition'
    # in logspace.  addLogs defined in dc.R. addLogs(a,b) = log(exp(a) + exp(b))
    for (i in 1:max.length) {
        if (check[i] >= 0) 
            indiv.LL.sum[i] <- addLogs(indiv.LL.sum[i], addition(alpha, beta, gamma, 
                delta, denom.ab, denom.gd, x[i], t.x[i], check[i]))
    }
    return(indiv.LL.sum)
}


bgbb.EstimateParameters <- function(rf.matrix, par.start = c(1, 1, 1, 1), max.param.value = 1000) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), par.start, "bgbb.EstimateParameters")
    
    bgbb.eLL <- function(params, rf.matrix, max.param.value) {
        params <- exp(params)
        params[params > max.param.value] <- max.param.value
        return(-1 * bgbb.rf.matrix.LL(params, rf.matrix))
    }
    logparams <- log(par.start)
    results <- optim(logparams, bgbb.eLL, rf.matrix = rf.matrix, max.param.value = max.param.value, 
        method = "L-BFGS-B")
    estimated.params <- exp(results$par)
    estimated.params[estimated.params > max.param.value] <- max.param.value
    return(estimated.params)
}

bgbb.pmf <- function(params, n, x) {
    
    max.length <- max(length(n), length(x))
    if (max.length%%length(n)) 
        warning("Maximum vector length not a multiple of the length of n")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.pmf")
    
    if (any(n < 0) || !is.numeric(n)) 
        stop("n must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    n <- rep(n, length.out = max.length)
    x <- rep(x, length.out = max.length)
    
    if (any(x > n)) {
        stop("bgbb.pmf was given x > n")
    }
    
    return(bgbb.pmf.General(params, 0, n, x))
}

bgbb.pmf.General <- function(params, n.cal, n.star, x.star) {
    
    max.length <- max(length(n.cal), length(n.star), length(x.star))
    
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    if (max.length%%length(n.star)) 
        warning("Maximum vector length not a multiple of the length of n.star")
    if (max.length%%length(x.star)) 
        warning("Maximum vector length not a multiple of the length of x.star")
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.pmf.General")
    
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    if (any(n.star < 0) || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not contain negative numbers.")
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    
    
    n.cal <- rep(n.cal, length.out = max.length)
    n.star <- rep(n.star, length.out = max.length)
    x.star <- rep(x.star, length.out = max.length)
    
    if (any(x.star > n.star)) {
        stop("bgbb.pmf.General was given x.star > n.star")
    }
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- rep(0, max.length)
    piece.1[x.star == 0] <- 1 - exp(lbeta(gamma, delta + n.cal[x.star == 0]) - lbeta(gamma, 
        delta))
    
    piece.2 <- exp(lchoose(n.star, x.star) + lbeta(alpha + x.star, beta + n.star - 
        x.star) - lbeta(alpha, beta) + lbeta(gamma, delta + n.cal + n.star) - lbeta(gamma, 
        delta))
    
    piece.3 = rep(0, max.length)
    rows.to.sum <- which(x.star <= n.star - 1)
    piece.3[rows.to.sum] <- unlist(sapply(rows.to.sum, function(index) {
        ii <- x.star[index]:(n.star[index] - 1)
        sum(exp(lchoose(ii, x.star[index]) + lbeta(alpha + x.star[index], beta + 
            ii - x.star[index]) - lbeta(alpha, beta) + lbeta(gamma + 1, delta + n.cal[index] + 
            ii) - lbeta(gamma, delta)))
    }))
    
    expectation <- piece.1 + piece.2 + piece.3
    return(expectation)
}

bgbb.Expectation <- function(params, n) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.Expectation")
    
    if (any(n < 0) || !is.numeric(n)) 
        stop("n must be numeric and may not contain negative numbers.")
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.one <- (alpha/(alpha + beta)) * (delta/(gamma - 1))
    piece.two <- exp(lgamma(gamma + delta) - lgamma(gamma + delta + n) + lgamma(1 + 
        delta + n) - lgamma(1 + delta))
    return(piece.one * (1 - piece.two))
}

bgbb.PlotFrequencyInCalibration <- function(params, rf.matrix, censor = NULL, plotZero = TRUE, 
    xlab = "Calibration period transactions", ylab = "Customers", title = "Frequency of Repeat Transactions") {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PlotFrequencyInCalibration")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotFrequencyInCalibration: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotFrequencyInCalibration: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotFrequencyInCalibration: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    max.n.cal <- max(n.cal)
    if (is.null(censor)) 
        censor <- max.n.cal
    total.custs <- sum(custs)
    actual.frequency <- rep(0, max.n.cal + 1)
    expected.frequency <- rep(0, max.n.cal + 1)
    
    for (ii in 0:max.n.cal) {
        actual.frequency[ii + 1] <- sum(custs[x == ii])
        expected.frequency[ii + 1] <- sum(unlist(sapply(unique(n.cal[n.cal >= ii]), 
            function(this.n.cal) {
                sum(custs[n.cal == this.n.cal]) * bgbb.pmf(params, this.n.cal, ii)
            })))
    }
    
    freq.comparison <- rbind(actual.frequency, expected.frequency)
    colnames(freq.comparison) <- 0:max.n.cal
    
    if (ncol(freq.comparison) <= censor) {
        censored.freq.comparison <- freq.comparison
    } else {
        ## Rename for easier coding
        fc <- freq.comparison
        ## Build censored freq comparison (cfc)
        cfc <- fc
        cfc <- cfc[, 1:(censor + 1)]
        cfc[1, (censor + 1)] <- sum(fc[1, (censor + 1):ncol(fc)])
        cfc[2, (censor + 1)] <- sum(fc[2, (censor + 1):ncol(fc)])
        
        censored.freq.comparison <- cfc
    }
    
    if (plotZero == FALSE) 
        censored.freq.comparison <- censored.freq.comparison[, -1]
    
    n.ticks <- ncol(censored.freq.comparison)
    if (plotZero == TRUE) {
        x.labels <- 0:(n.ticks - 1)
        if (censor < ncol(freq.comparison) - 1) 
            x.labels[n.ticks] <- paste(n.ticks - 1, "+", sep = "")
    } else {
        x.labels <- 1:(n.ticks)
        if (censor < ncol(freq.comparison) - 1) 
            x.labels[n.ticks] <- paste(n.ticks, "+", sep = "")
    }
    colnames(censored.freq.comparison) <- x.labels
    
    ylim <- c(0, ceiling(max(c(censored.freq.comparison[1, ], censored.freq.comparison[2, 
        ])) * 1.1))
    
    barplot(censored.freq.comparison, beside = TRUE, ylim = ylim, main = title, xlab = xlab, 
        ylab = ylab, col = 1:2)
    legend("top", legend = c("Actual", "Model"), col = 1:2, lwd = 2, cex = 0.75)
    
    return(censored.freq.comparison)
}

bgbb.PlotTrackingCum <- function(params, rf.matrix, actual.cum.repeat.transactions, 
    xlab = "Time", ylab = "Cumulative Transactions", xticklab = NULL, title = "Tracking Cumulative Transactions") {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PlotTrackingCum")
    if (any(actual.cum.repeat.transactions < 0) || !is.numeric(actual.cum.repeat.transactions)) 
        stop("actual.cum.repeat.transactions must be numeric and may not contain negative numbers.")
    
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotTrackingCum: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotTrackingCum: rf.matrix must have a column for the number of customers represented by each row, labelled \"custs\""))
    
    actual <- actual.cum.repeat.transactions
    n.periods <- length(actual)
    
    cust.birth.periods <- max(n.cal) - n.cal
    
    expected <- sapply(1:n.periods, function(interval) {
        if (interval <= min(cust.birth.periods)) 
            return(0)
        sum(bgbb.Expectation(params, interval - cust.birth.periods[cust.birth.periods <= 
            interval]) * custs[cust.birth.periods <= interval])
    })
    
    pur.comparison <- rbind(actual, expected)
    
    ylim <- c(0, max(c(actual, expected)) * 1.05)
    plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = ylim, 
        main = title)
    lines(expected, lty = 2, col = 2)
    if (is.null(xticklab) == FALSE) {
        if (ncol(pur.comparison) != length(xticklab)) {
            stop("Plot error, xticklab does not have the correct size in bgbb.PlotTrackingCum")
        }
        axis(1, at = 1:ncol(pur.comparison), labels = xticklab)
    }
    if (is.null(n.cal) == FALSE) {
        abline(v = max(n.cal), lty = 2)
    }
    legend("bottomright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1, 
        cex = 0.75)
    
    return(pur.comparison)
}

bgbb.PlotTrackingInc <- function(params, rf.matrix, actual.inc.repeat.transactions, 
    xlab = "Time", ylab = "Transactions", xticklab = NULL, title = "Tracking Incremental Transactions") {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PlotTrackingInc")
    if (any(actual.inc.repeat.transactions < 0) || !is.numeric(actual.inc.repeat.transactions)) 
        stop("actual.inc.repeat.transactions must be numeric and may not contain negative numbers.")
    
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotTrackingInc: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotTrackingInc: rf.matrix must have a column for the number of customers represented by each row, labelled \"custs\""))
    
    actual <- actual.inc.repeat.transactions
    n.periods <- length(actual)
    
    cust.birth.periods <- max(n.cal) - n.cal
    
    expected.cumulative <- sapply(1:n.periods, function(interval) {
        if (interval <= min(cust.birth.periods)) 
            return(0)
        sum(bgbb.Expectation(params, interval - cust.birth.periods[cust.birth.periods <= 
            interval]) * custs[cust.birth.periods <= interval])
    })
    
    expected <- dc.CumulativeToIncremental(expected.cumulative)
    
    pur.comparison <- rbind(actual, expected)
    
    ylim <- c(0, max(c(actual, expected)) * 1.05)
    plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = ylim, 
        main = title)
    lines(expected, lty = 2, col = 2)
    if (is.null(xticklab) == FALSE) {
        if (ncol(pur.comparison) != length(xticklab)) {
            stop("Plot error, xticklab does not have the correct size in bgbb.PlotTrackingInc")
        }
        axis(1, at = 1:ncol(pur.comparison), labels = xticklab)
    }
    if (is.null(n.cal) == FALSE) {
        abline(v = max(n.cal), lty = 2)
    }
    legend("topright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1, 
        cex = 0.75)
    
    return(pur.comparison)
}


bgbb.ConditionalExpectedTransactions <- function(params, n.cal, n.star, x, t.x) {
    
    max.length <- max(length(x), length(t.x), length(n.cal), length(n.star))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    if (max.length%%length(n.star)) 
        warning("Maximum vector length not a multiple of the length of n.star")
    
    dc.check.model.params(c("r", "alpha", "s", "beta"), params, "bgbb.ConditionalExpectedTransactions")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    if (any(n.star < 0) || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not contain negative numbers.")
    
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    n.star <- rep(n.star, length.out = max.length)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- 1/exp(bgbb.LL(params, x, t.x, n.cal))
    piece.2 <- exp(lbeta(alpha + x + 1, beta + n.cal - x) - lbeta(alpha, beta))
    piece.3 <- delta/(gamma - 1)
    piece.4 <- exp(lgamma(gamma + delta) - lgamma(1 + delta))
    piece.5 <- exp(lgamma(1 + delta + n.cal) - lgamma(gamma + delta + n.cal))
    piece.6 <- exp(lgamma(1 + delta + n.cal + n.star) - lgamma(gamma + delta + n.cal + 
        n.star))
    
    expected.frequency <- piece.1 * piece.2 * piece.3 * piece.4 * (piece.5 - piece.6)
    
    which.is.nan <- is.nan(expected.frequency)
    if (sum(which.is.nan) > 0) {
        error.msg.long <- paste("numerical error, parameters exploded in bgbb.ConditionalExpectedTransactions", 
            "params:", alpha, beta, gamma, delta, "n.cal:", n.cal[which.is.nan], 
            "n.star:", n.star[which.is.nan], "x:", x[which.is.nan], "t.x:", t.x[which.is.nan], 
            "...", "piece.1:", piece.1[which.is.nan], "piece.2:", piece.2[which.is.nan], 
            "piece.3:", piece.3[which.is.nan], "piece.4:", piece.4[which.is.nan], 
            "piece.5:", piece.5[which.is.nan], "piece.6:", piece.6[which.is.nan])
        stop(error.msg.long)
    }
    return(expected.frequency)
}

bgbb.PlotFreqVsConditionalExpectedFrequency <- function(params, n.star, rf.matrix, 
    x.star, trunc = NULL, xlab = "Calibration period transactions", ylab = "Holdout period transactions", 
    xticklab = NULL, title = "Conditional Expectation") {
    
    if (length(x.star) != nrow(rf.matrix)) 
        stop("x.star must have the same number of entries as rows in rf.matrix")
    if (!(length(n.star) == 1 || length(n.star) == nrow(rf.matrix))) 
        stop("n.star must be a single value or have as many entries as rows in rf.matrix")
    
    dc.check.model.params(c("r", "alpha", "s", "beta"), params, "bgbb.PlotFreqVsConditionalExpectedFrequency")
    
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    if (any(n.star < 0) || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not contain negative numbers.")
    
    n.star <- rep(n.star, length.out = nrow(rf.matrix))
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotFreqVsConditionalExpectedFrequency: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    if (is.null(trunc)) 
        trunc <- max(n.cal)
    
    if (trunc > max(n.cal)) {
        warning("The truncation number provided in bgbb.PlotFreqVsConditionalExpectedFrequency was greater than the maximum number of possible transactions. It has been reduced to ", 
            max(n.cal))
        trunc = max(n.cal)
    }
    
    actual.freq <- rep(0, max(n.cal) + 1)
    expected.freq <- rep(0, max(n.cal) + 1)
    bin.size <- rep(0, max(n.cal) + 1)
    
    for (ii in 0:max(n.cal)) {
        bin.size[ii + 1] <- sum(custs[x == ii])
        actual.freq[ii + 1] <- sum(x.star[x == ii])
        expected.freq[ii + 1] <- sum(bgbb.ConditionalExpectedTransactions(params, 
            n.cal[x == ii], n.star[x == ii], ii, t.x[x == ii]) * custs[x == ii])
    }
    
    comparison <- rbind(actual.freq/bin.size, expected.freq/bin.size, bin.size)
    colnames(comparison) <- paste("freq.", 0:max(n.cal), sep = "")
    
    if (is.null(xticklab) == FALSE) {
        x.labels <- xticklab
    } else {
        if ((trunc + 1) < ncol(comparison)) {
            x.labels <- 0:(trunc)
        } else {
            x.labels <- 0:(ncol(comparison) - 1)
        }
    }
    
    actual.freq <- comparison[1, 1:(trunc + 1)]
    expected.freq <- comparison[2, 1:(trunc + 1)]
    
    custs.in.plot <- sum(comparison[3, 1:(trunc + 1)])
    if (custs.in.plot < 0.9 * sum(custs)) {
        warning("Less than 90% of customers are represented in your plot (", custs.in.plot, 
            " of ", sum(custs), " are plotted).")
    }
    
    ylim <- c(0, ceiling(max(c(actual.freq, expected.freq)) * 1.1))
    plot(actual.freq, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, 
        ylab = ylab, main = title)
    lines(expected.freq, lty = 2, col = 2)
    
    axis(1, at = 1:(trunc + 1), labels = x.labels)
    legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    
    return(comparison)
}


bgbb.PlotRecVsConditionalExpectedFrequency <- function(params, n.star, rf.matrix, 
    x.star, trunc = NULL, xlab = "Calibration period recency", ylab = "Holdout period transactions", 
    xticklab = NULL, title = "Conditional Expected Transactions by Recency") {
    
    if (length(x.star) != nrow(rf.matrix)) 
        stop("x.star must have the same number of entries as rows in rf.matrix")
    if (!(length(n.star) == 1 || length(n.star) == nrow(rf.matrix))) 
        stop("n.star must be a single value or have as many entries as rows in rf.matrix")
    
    dc.check.model.params(c("r", "alpha", "s", "beta"), params, "bgbb.PlotRecVsConditionalExpectedFrequency")
    
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    if (any(n.star < 0) || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not contain negative numbers.")
    
    n.star <- rep(n.star, length.out = nrow(rf.matrix))
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    tryCatch(custs <- rf.matrix[, "custs"], error = function(e) stop("Error in bgbb.PlotRecVsConditionalExpectedFrequency: rf.matrix must have a column for the number of customers that have each combination of \"x\", \"t.x\", and \"n.cal\", labelled \"custs\""))
    
    if (is.null(trunc)) 
        trunc <- max(n.cal)
    
    if (trunc > max(n.cal)) {
        warning("The truncation number provided in bgbb.PlotRecVsConditionalExpectedFrequency was greater than the maximum number of possible transactions. It has been reduced to ", 
            max(n.cal))
        trunc = max(n.cal)
    }
    
    actual.freq <- rep(0, max(n.cal) + 1)
    expected.freq <- rep(0, max(n.cal) + 1)
    bin.size <- rep(0, max(n.cal) + 1)
    
    for (ii in 0:max(n.cal)) {
        bin.size[ii + 1] <- sum(custs[t.x == ii])
        actual.freq[ii + 1] <- sum(x.star[t.x == ii])
        expected.freq[ii + 1] <- sum(bgbb.ConditionalExpectedTransactions(params, 
            n.cal[t.x == ii], n.star[t.x == ii], x[t.x == ii], ii) * custs[t.x == 
            ii])
    }
    
    comparison <- rbind(actual.freq/bin.size, expected.freq/bin.size, bin.size)
    colnames(comparison) <- paste("rec.", 0:max(n.cal), sep = "")
    
    custs.in.plot <- sum(comparison[3, 1:(trunc + 1)])
    if (custs.in.plot < 0.9 * sum(custs)) {
        warning("Less than 90% of customers are represented in your plot (", custs.in.plot, 
            " of ", sum(custs), " are plotted).")
    }
    
    if (is.null(xticklab) == FALSE) {
        x.labels <- xticklab
    } else {
        if ((trunc + 1) < ncol(comparison)) {
            x.labels <- 0:(trunc)
        } else {
            x.labels <- 0:(ncol(comparison) - 1)
        }
    }
    
    actual.freq <- comparison[1, 1:(trunc + 1)]
    expected.freq <- comparison[2, 1:(trunc + 1)]
    
    ylim <- c(0, ceiling(max(c(actual.freq, expected.freq)) * 1.1))
    plot(actual.freq, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, 
        ylab = ylab, main = title)
    lines(expected.freq, lty = 2, col = 2)
    
    axis(1, at = 1:(trunc + 1), labels = x.labels)
    legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    
    return(comparison)
}

bgbb.PosteriorMeanLmProductMoment <- function(params, l, m, x, t.x, n.cal) {
    
    max.length <- max(length(x), length(t.x), length(n.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PosteriorMeanLmProduct")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    if (l < 0 || length(l) != 1 || !is.numeric(l)) 
        stop("l must be a single numeric value and may not be less than 0.")
    if (m < 0 || length(m) != 1 || !is.numeric(m)) 
        stop("m must be a single numeric value and may not be less than 0.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- exp(lbeta(alpha + l, beta) - lbeta(alpha, beta) + lbeta(gamma + m, 
        delta) - lbeta(gamma, delta))
    piece.2 <- exp(bgbb.LL(c(alpha + l, beta, gamma + m, delta), x, t.x, n.cal))
    piece.3 <- exp(bgbb.LL(params, x, t.x, n.cal))
    
    mean <- piece.1 * (piece.2/piece.3)
    
    return(mean)
}

bgbb.PosteriorMeanTransactionRate <- function(params, x, t.x, n.cal) {
    
    max.length <- max(length(x), length(t.x), length(n.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PosteriorMeanTransactionRate")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    
    mean.transaction.rate <- bgbb.PosteriorMeanLmProductMoment(params, 1, 0, x, t.x, 
        n.cal)
    return(mean.transaction.rate)
}

bgbb.rf.matrix.PosteriorMeanTransactionRate <- function(params, rf.matrix) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.rf.matrix.PosteriorMeanTransactionRate")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanTransactionRate: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanTransactionRate: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanTransactionRate: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    
    return(bgbb.PosteriorMeanTransactionRate(params, x, t.x, n.cal))
}

bgbb.PosteriorMeanDropoutRate <- function(params, x, t.x, n.cal) {
    
    max.length <- max(length(x), length(t.x), length(n.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PosteriorMeanDropoutRate")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    
    mean.dropout.rate <- bgbb.PosteriorMeanLmProductMoment(params, 0, 1, x, t.x, 
        n.cal)
    return(mean.dropout.rate)
}

bgbb.rf.matrix.PosteriorMeanDropoutRate <- function(params, rf.matrix) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.rf.matrix.PosteriorMeanDropoutRate")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanDropoutRate: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanDropoutRate: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.PosteriorMeanDropoutRate: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    
    return(bgbb.PosteriorMeanDropoutRate(params, x, t.x, n.cal))
}

bgbb.DERT <- function(params, x, t.x, n.cal, d) {
    
    max.length <- max(length(x), length(t.x), length(n.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    
    dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.DERT")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- exp(lbeta(alpha + x + 1, beta + n.cal - x) - lbeta(alpha, beta))
    piece.2 <- exp(lbeta(gamma, delta + n.cal + 1) - lbeta(gamma, delta))/(1 + d)
    piece.3 <- Re(hypergeo(1, delta + n.cal + 1, gamma + delta + n.cal + 1, 1/(1 + 
        d)))
    piece.4 <- exp(bgbb.LL(params, x, t.x, n.cal))
    
    dert <- piece.1 * piece.2 * (piece.3/piece.4)
    
    return(dert)
}

bgbb.rf.matrix.DERT <- function(params, rf.matrix, d) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.rf.matrix.DERT")
    
    tryCatch(x <- rf.matrix[, "x"], error = function(e) stop("Error in bgbb.rf.matrix.DERT: rf.matrix must have a frequency column labelled \"x\""))
    tryCatch(t.x <- rf.matrix[, "t.x"], error = function(e) stop("Error in bgbb.rf.matrix.DERT: rf.matrix must have a recency column labelled \"t.x\""))
    tryCatch(n.cal <- rf.matrix[, "n.cal"], error = function(e) stop("Error in bgbb.rf.matrix.DERT: rf.matrix must have a column for number of transaction opportunities in the calibration period, labelled \"n.cal\""))
    
    return(bgbb.DERT(params, x, t.x, n.cal, d))
}


bgbb.PAlive <- function(params, x, t.x, n.cal) {
    
    max.length <- max(length(x), length(t.x), length(n.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(n.cal)) 
        warning("Maximum vector length not a multiple of the length of n.cal")
    
    dc.check.model.params(c("r", "alpha", "s", "beta"), params, "bgbb.PAlive")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(n.cal < 0) || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not contain negative numbers.")
    
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    n.cal <- rep(n.cal, length.out = max.length)
    
    alpha <- params[1]
    beta <- params[2]
    gamma <- params[3]
    delta <- params[4]
    
    piece.1 <- exp(lbeta(alpha + x, beta + n.cal - x) - lbeta(alpha, beta) + lbeta(gamma, 
        delta + n.cal + 1) - lbeta(gamma, delta))
    piece.2 <- 1/exp(bgbb.LL(params, x, t.x, n.cal))
    
    p.alive <- piece.1 * piece.2
    
    return(p.alive)
}

bgbb.HeatmapHoldoutExpectedTrans <- function(params, n.cal, n.star, xlab = "Recency", 
    ylab = "Frequency", xticklab = NULL, title = "Heatmap of Conditional Expected Transactions") {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.HeatmapHoldoutExpectedTrans")
    if (n.cal < 0 || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not be negative.")
    if (n.star < 0 || !is.numeric(n.star)) 
        stop("n.star must be numeric and may not be negative.")
    
    heatmap.mx <- matrix(0, n.cal + 1, n.cal + 1)
    heatmap.mx[1, 1] <- bgbb.ConditionalExpectedTransactions(params, n.cal, n.star, 
        0, 0)
    for (xx in 1:n.cal) {
        for (tt in 1:n.cal) {
            if (xx <= tt) {
                expected.trans <- bgbb.ConditionalExpectedTransactions(params, n.cal, 
                  n.star, xx, tt)
                heatmap.mx[xx + 1, tt + 1] <- expected.trans
            }
        }
    }
    if (is.null(xticklab) == TRUE) {
        xticklab <- 0:n.cal
    }
    colnames(heatmap.mx) <- xticklab
    rownames(heatmap.mx) <- 0:n.cal
    heatmap(heatmap.mx, Rowv = NA, Colv = NA, col = gray(8:2/9), scale = "none", 
        ylab = ylab, xlab = xlab, main = title, verbose = TRUE)
    return(heatmap.mx)
}

bgbb.PlotFrequencyInHoldout <- function(params, n.cal, rf.matrix.holdout, censor = NULL, 
    plotZero = TRUE, title = "Frequency of Repeat Transactions", xlab = "Holdout period transactions", 
    ylab = "Customers") {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PlotFrequencyInHoldout")
    if (n.cal < 0 || !is.numeric(n.cal)) 
        stop("n.cal must be numeric and may not be negative.")
    
    tryCatch(x.star <- rf.matrix.holdout[, "x.star"], error = function(e) stop("Error in bgbb.PlotFrequencyInHoldout: rf.matrix.holdout must have a frequency column labelled \"x.star\""))
    tryCatch(n.star <- rf.matrix.holdout[, "n.star"], error = function(e) stop("Error in bgbb.PlotFrequencyInHoldout: rf.matrix.holdout must have a column with the number of transaction opportunities for each group, labelled \"n.star\""))
    tryCatch(custs <- rf.matrix.holdout[, "custs"], error = function(e) stop("Error in bgbb.PlotFrequencyInHoldout: rf.matrix.holdout must have a column for the number of customers represented by each row, labelled \"custs\""))
    
    max.n.star <- max(n.star)
    if (is.null(censor)) 
        censor <- max.n.star
    total.custs <- sum(custs)
    actual.frequency <- rep(0, max.n.star + 1)
    expected.frequency <- rep(0, max.n.star + 1)
    
    for (ii in 0:max.n.star) {
        actual.frequency[ii + 1] <- sum(custs[x.star == ii])
        expected.frequency[ii + 1] <- sum(unlist(sapply(unique(n.star[n.star >= ii]), 
            function(this.n.star) {
                sum(custs[n.star == this.n.star]) * bgbb.pmf.General(params, n.cal, 
                  this.n.star, ii)
            })))
    }
    
    freq.comparison <- rbind(actual.frequency, expected.frequency)
    colnames(freq.comparison) <- 0:max.n.star
    
    if (ncol(freq.comparison) <= censor) {
        censored.freq.comparison <- freq.comparison
    } else {
        ## Rename for easier coding
        fc <- freq.comparison
        ## Build censored freq comparison (cfc)
        cfc <- fc
        cfc <- cfc[, 1:(censor + 1)]
        cfc[1, (censor + 1)] <- sum(fc[1, (censor + 1):ncol(fc)])
        cfc[2, (censor + 1)] <- sum(fc[2, (censor + 1):ncol(fc)])
        
        if (plotZero == FALSE) {
            cfc <- cfc[, -1]
        }
        
        censored.freq.comparison <- cfc
    }
    
    if (plotZero == TRUE) {
        x.labels <- 0:(ncol(censored.freq.comparison) - 1)
    } else {
        x.labels <- 1:(ncol(censored.freq.comparison))
    }
    
    if (censor < ncol(freq.comparison) - 1) {
        x.labels[(censor + 1)] <- paste(censor, "+", sep = "")
    }
    colnames(censored.freq.comparison) <- x.labels
    
    barplot(censored.freq.comparison, beside = TRUE, main = title, xlab = xlab, ylab = ylab, 
        col = 1:2)
    legend("topright", legend = c("Actual", "Model"), col = 1:2, lwd = 2)
    
    return(censored.freq.comparison)
}

bgbb.PlotTransactionRateHeterogeneity <- function(params) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PlotTransactionRateHeterogeneity")
    
    alpha <- params[1]
    beta <- params[2]
    x.axis.ticks <- 0.01 * 0:100
    heterogeneity <- dbeta(x.axis.ticks, alpha, beta)
    plot(x.axis.ticks, heterogeneity, type = "l", xlab = "Transaction Rate", ylab = "Density", 
        main = "Heterogeneity in Transaction Rate")
    rate.mean <- round(alpha/(alpha + beta), 4)
    rate.var <- round((alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1)), 4)
    mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
    mtext(mean.var.label, side = 3)
    return(heterogeneity)
}

bgbb.PlotDropoutRateHeterogeneity <- function(params) {
    
    dc.check.model.params(c("alpha", "beta", "gamma", "delta"), params, "bgbb.PlotDropoutRateHeterogeneity")
    
    alpha <- params[3]
    beta <- params[4]
    x.axis.ticks <- 0.01 * 0:100
    heterogeneity <- dbeta(x.axis.ticks, alpha, beta)
    plot(x.axis.ticks, heterogeneity, type = "l", xlab = "Dropout rate", ylab = "Density", 
        main = "Heterogeneity in Dropout Rate")
    rate.mean <- round(alpha/(alpha + beta), 4)
    rate.var <- round((alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1)), 4)
    mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
    mtext(mean.var.label, side = 3)
    return(heterogeneity)
} 
