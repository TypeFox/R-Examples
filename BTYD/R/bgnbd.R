################################################## BG/NBD estimation, visualization functions

library(hypergeo)

bgnbd.cbs.LL <- function(params, cal.cbs) {
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.cbs.LL")
    
    tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in bgnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""))
    tryCatch(t.x <- cal.cbs[, "t.x"], error = function(e) stop("Error in bgnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""))
    tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in bgnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
    
    if ("custs" %in% colnames(cal.cbs)) {
        
        many_rows = function(vec, nreps) {
            return(rep(1, nreps) %*% t.default(vec))
        }
        
        custs <- cal.cbs[, "custs"]
        logvec = (1:length(custs)) * (custs > 1)
        logvec = logvec[logvec > 0]
        M = sum(logvec > 0)
        for (i in 1:M) {
            cal.cbs = rbind(cal.cbs, many_rows(cal.cbs[logvec[i], ], custs[logvec[i]] - 
                1))
        }
        x = cal.cbs[, "x"]
        t.x = cal.cbs[, "t.x"]
        T.cal = cal.cbs[, "T.cal"]
    }
    
    return(sum(bgnbd.LL(params, x, t.x, T.cal)))
}

bgnbd.LL <- function(params, x, t.x, T.cal) {
    
    beta.ratio = function(a, b, x, y) {
        exp(lgamma(a) + lgamma(b) - lgamma(a + b) - lgamma(x) - lgamma(y) + lgamma(x + 
            y))
    }
    
    max.length <- max(length(x), length(t.x), length(T.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(T.cal)) 
        warning("Maximum vector length not a multiple of the length of T.cal")
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.LL")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(T.cal < 0) || !is.numeric(T.cal)) 
        stop("T.cal must be numeric and may not contain negative numbers.")
    
    r = params[1]
    alpha = params[2]
    a = params[3]
    b = params[4]
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    T.cal <- rep(T.cal, length.out = max.length)
    
    A = r * log(alpha) + lgamma(r + x) - lgamma(r) - (r + x) * log(alpha + t.x)
    B = beta.ratio(a, b + x, a, b) * ((alpha + t.x)/(alpha + T.cal))^(r + x) + as.numeric((x > 
        0)) * beta.ratio(a + 1, b + x - 1, a, b)
    LL = sum(A + log(B))
    
    return(LL)
}

bgnbd.compress.cbs <- function(cbs, rounding = 3) {
    
    if (!("x" %in% colnames(cbs))) 
        stop("Error in bgnbd.compress.cbs: cbs must have a frequency column labelled \"x\"")
    if (!("t.x" %in% colnames(cbs))) 
        stop("Error in bgnbd.compress.cbs: cbs must have a recency column labelled \"t.x\"")
    if (!("T.cal" %in% colnames(cbs))) 
        stop("Error in bgnbd.compress.cbs: cbs must have a column for length of time observed labelled \"T.cal\"")
    
    orig.rows <- nrow(cbs)
    
    if (!("custs" %in% colnames(cbs))) {
        custs <- rep(1, nrow(cbs))
        cbs <- cbind(cbs, custs)
    }
    
    other.colnames <- colnames(cbs)[!(colnames(cbs) %in% c("x", "t.x", "T.cal"))]
    
    ## Round x, t.x and T.cal to the desired level
    cbs[, c("x", "t.x", "T.cal")] <- round(cbs[, c("x", "t.x", "T.cal")], rounding)
    
    ## Aggregate every column that is not x, t.x or T.cal by those columns. Do this by
    ## summing entries which have the same x, t.x and T.cal.
    cbs <- as.matrix(aggregate(cbs[, !(colnames(cbs) %in% c("x", "t.x", "T.cal"))], 
        by = list(x = cbs[, "x"], t.x = cbs[, "t.x"], T.cal = cbs[, "T.cal"]), sum))
    
    colnames(cbs) <- c("x", "t.x", "T.cal", other.colnames)
    final.rows <- nrow(cbs)
    message("Data reduced from ", orig.rows, " rows to ", final.rows, " rows.")
    return(cbs)
}

bgnbd.EstimateParameters <- function(cal.cbs, par.start = c(1, 3, 1, 3), max.param.value = 10000) {
    
    dc.check.model.params(c("r", "alpha", "a", "b"), par.start, "bgnbd.EstimateParameters")
    
    bgnbd.eLL <- function(params, cal.cbs, max.param.value) {
        params <- exp(params)
        params[params > max.param.value] = max.param.value
        return(-1 * bgnbd.cbs.LL(params, cal.cbs))
    }
    logparams = log(par.start)
    results = optim(logparams, bgnbd.eLL, cal.cbs = cal.cbs, max.param.value = max.param.value, 
        method = "L-BFGS-B")
    estimated.params <- exp(results$par)
    estimated.params[estimated.params > max.param.value] <- max.param.value
    return(estimated.params)
}

bgnbd.pmf <- function(params, t, x) {
    max.length <- max(length(t), length(x))
    if (max.length%%length(t)) 
        warning("Maximum vector length not a multiple of the length of t")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.pmf")
    
    if (any(t < 0) || !is.numeric(t)) 
        stop("t must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    t. <- rep(t, length.out = max.length)
    x <- rep(x, length.out = max.length)
    
    return(bgnbd.pmf.General(params, 0, t, x))
}

bgnbd.pmf.General <- function(params, t.start, t.end, x) {
    
    max.length = max(length(t.start), length(t.end), length(x))
    
    if (max.length%%length(t.start)) 
        warning("Maximum vector length not a multiple of the length of t.start")
    if (max.length%%length(t.end)) 
        warning("Maximum vector length not a multiple of the length of t.end")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.pmf.General")
    
    if (any(t.start < 0) || !is.numeric(t.start)) 
        stop("t.start must be numeric and may not contain negative numbers.")
    if (any(t.end < 0) || !is.numeric(t.end)) 
        stop("t.end must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    
    t.start = rep(t.start, length.out = max.length)
    t.end = rep(t.end, length.out = max.length)
    x = rep(x, length.out = max.length)
    
    if (any(t.start > t.end)) {
        stop("Error in bgnbd.pmf.General: t.start > t.end.")
    }
    r <- params[1]
    alpha <- params[2]
    a <- params[3]
    b <- params[4]
    equation.part.0 <- rep(0, max.length)
    t = t.end - t.start
    term3 = rep(0, max.length)
    term1 = beta(a, b + x)/beta(a, b) * gamma(r + x)/gamma(r)/factorial(x) * ((alpha/(alpha + 
        t))^r) * ((t/(alpha + t))^x)
    
    for (i in 1:max.length) {
        if (x[i] > 0) {
            ii = c(0:(x[i] - 1))
            summation.term = sum(gamma(r + ii)/gamma(r)/factorial(ii) * ((t[i]/(alpha + 
                t[i]))^ii))
            term3[i] = 1 - (((alpha/(alpha + t[i]))^r) * summation.term)
        }
    }
    term2 = as.numeric(x > 0) * beta(a + 1, b + x - 1)/beta(a, b) * term3
    
    return(term1 + term2)
}


bgnbd.ConditionalExpectedTransactions <- function(params, T.star, x, t.x, T.cal) {
    
    h2f1 <- function(a, b, c, z) {
        lenz <- length(z)
        j = 0
        uj <- 1:lenz
        uj <- uj/uj
        y <- uj
        lteps <- 0
        
        while (lteps < lenz) {
            lasty <- y
            j <- j + 1
            uj <- uj * (a + j - 1) * (b + j - 1)/(c + j - 1) * z/j
            y <- y + uj
            lteps <- sum(y == lasty)
        }
        return(y)
    }
    
    max.length <- max(length(T.star), length(x), length(t.x), length(T.cal))
    
    if (max.length%%length(T.star)) 
        warning("Maximum vector length not a multiple of the length of T.star")
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(T.cal)) 
        warning("Maximum vector length not a multiple of the length of T.cal")
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.ConditionalExpectedTransactions")
    
    if (any(T.star < 0) || !is.numeric(T.star)) 
        stop("T.star must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(T.cal < 0) || !is.numeric(T.cal)) 
        stop("T.cal must be numeric and may not contain negative numbers.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    T.cal <- rep(T.cal, length.out = max.length)
    
    r = params[1]
    alpha = params[2]
    a = params[3]
    b = params[4]
    term1 <- ((a + b + x - 1)/(a - 1))
    term2 <- 1 - ((alpha + T.cal)/(alpha + T.cal + T.star))^(r + x) * h2f1(r + x, 
        b + x, a + b + x - 1, T.star/(alpha + T.cal + T.star))
    term3 <- 1 + as.numeric(x > 0) * (a/(b + x - 1)) * ((alpha + T.cal)/(alpha + 
        t.x))^(r + x)
    out <- term1 * term2/term3
    return(out)
}

bgnbd.PAlive <- function(params, x, t.x, T.cal) {
    
    max.length <- max(length(x), length(t.x), length(T.cal))
    
    if (max.length%%length(x)) 
        warning("Maximum vector length not a multiple of the length of x")
    if (max.length%%length(t.x)) 
        warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length%%length(T.cal)) 
        warning("Maximum vector length not a multiple of the length of T.cal")
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PAlive")
    
    if (any(x < 0) || !is.numeric(x)) 
        stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x)) 
        stop("t.x must be numeric and may not contain negative numbers.")
    if (any(T.cal < 0) || !is.numeric(T.cal)) 
        stop("T.cal must be numeric and may not contain negative numbers.")
    
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    T.cal <- rep(T.cal, length.out = max.length)
    
    r = params[1]
    alpha = params[2]
    a = params[3]
    b = params[4]
    term1 = (a/(b + x - 1)) * ((alpha + T.cal)/(alpha + t.x))^(r + x)
    return(1/(1 + as.numeric(x > 0) * term1))
}

bgnbd.Expectation <- function(params, t) {
    
    h2f1 <- function(a, b, c, z) {
        lenz <- length(z)
        j = 0
        uj <- 1:lenz
        uj <- uj/uj
        y <- uj
        lteps <- 0
        
        while (lteps < lenz) {
            lasty <- y
            j <- j + 1
            uj <- uj * (a + j - 1) * (b + j - 1)/(c + j - 1) * z/j
            y <- y + uj
            lteps <- sum(y == lasty)
        }
        return(y)
    }
    
    if (any(t < 0) || !is.numeric(t)) 
        stop("t must be numeric and may not contain negative numbers.")
    r = params[1]
    alpha = params[2]
    a = params[3]
    b = params[4]
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.Expectation")
    
    term1 = (a + b - 1)/(a - 1)
    term2 = (alpha/(alpha + t))^r
    term3 = h2f1(r, b, a + b - 1, t/(alpha + t))
    output = term1 * (1 - term2 * term3)
    
    return(output)
}

bgnbd.PlotFrequencyInCalibration <- function(params, cal.cbs, censor, plotZero = TRUE, 
    xlab = "Calibration period transactions", ylab = "Customers", title = "Frequency of Repeat Transactions") {
    
    tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in bgnbd.PlotFrequencyInCalibration: cal.cbs must have a frequency column labelled \"x\""))
    tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in bgnbd.PlotFrequencyInCalibration: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotFrequencyInCalibration")
    if (censor > max(x)) 
        stop("censor too big (> max freq) in PlotFrequencyInCalibration.")
    
    x = cal.cbs[, "x"]
    T.cal = cal.cbs[, "T.cal"]
    n.x <- rep(0, max(x) + 1)
    ncusts = nrow(cal.cbs)
    for (ii in unique(x)) {
        # Get number of customers to buy n.x times, over the grid of all possible n.x
        # values (no censoring)
        n.x[ii + 1] <- sum(ii == x)
    }
    n.x.censor <- sum(n.x[(censor + 1):length(n.x)])
    n.x.actual <- c(n.x[1:censor], n.x.censor)  # This upper truncates at censor (ie. if censor=7, 8 categories: {0, 1, ..., 6, 7+}).
    T.value.counts <- table(T.cal)  # This is the table of counts of all time durations from customer birth to end of calibration period.
    T.values <- as.numeric(names(T.value.counts))  # These are all the unique time durations from customer birth to end of calibration period.
    n.T.values <- length(T.values)  # These are the number of time durations we need to consider.
    n.x.expected <- rep(0, length(n.x.actual))  # We'll store the probabilities in here.
    n.x.expected.all <- rep(0, max(x) + 1)  # We'll store the probabilities in here.
    
    for (ii in 0:max(x)) {
        # We want to run over the probability of each transaction amount.
        this.x.expected = 0
        for (T.idx in 1:n.T.values) {
            # We run over all people who had all time durations.
            T = T.values[T.idx]
            if (T == 0) 
                next
            n.T = T.value.counts[T.idx]  # This is the number of customers who had this time duration.
            prob.of.this.x.for.this.T = bgnbd.pmf(params, T, ii)
            expected.given.x.and.T = n.T * prob.of.this.x.for.this.T
            this.x.expected = this.x.expected + expected.given.x.and.T
        }
        n.x.expected.all[ii + 1] = this.x.expected
    }
    n.x.expected[1:censor] = n.x.expected.all[1:censor]
    n.x.expected[censor + 1] = sum(n.x.expected.all[(censor + 1):(max(x) + 1)])
    
    col.names <- paste(rep("freq", length(censor + 1)), (0:censor), sep = ".")
    col.names[censor + 1] <- paste(col.names[censor + 1], "+", sep = "")
    censored.freq.comparison <- rbind(n.x.actual, n.x.expected)
    colnames(censored.freq.comparison) <- col.names
    cfc.plot <- censored.freq.comparison
    if (plotZero == FALSE) 
        cfc.plot <- cfc.plot[, -1]
    n.ticks <- ncol(cfc.plot)
    if (plotZero == TRUE) {
        x.labels <- 0:(n.ticks - 1)
        x.labels[n.ticks] <- paste(n.ticks - 1, "+", sep = "")
    }
    ylim <- c(0, ceiling(max(cfc.plot) * 1.1))
    barplot(cfc.plot, names.arg = x.labels, beside = TRUE, ylim = ylim, main = title, 
        xlab = xlab, ylab = ylab, col = 1:2)
    legend("topright", legend = c("Actual", "Model"), col = 1:2, lwd = 2)
    return(censored.freq.comparison)
}


bgnbd.PlotFreqVsConditionalExpectedFrequency <- function(params, T.star, cal.cbs, 
    x.star, censor, xlab = "Calibration period transactions", ylab = "Holdout period transactions", 
    xticklab = NULL, title = "Conditional Expectation") {
    tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in bgnbd.PlotFreqVsConditionalExpectedFrequency: cal.cbs must have a frequency column labelled \"x\""))
    tryCatch(t.x <- cal.cbs[, "t.x"], error = function(e) stop("Error in bgnbd.PlotFreqVsConditionalExpectedFrequency: cal.cbs must have a recency column labelled \"t.x\""))
    tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in bgnbd.PlotFreqVsConditionalExpectedFrequency: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotFreqVsConditionalExpectedFrequency")
    if (censor > max(x)) 
        stop("censor too big (> max freq) in PlotFreqVsConditionalExpectedFrequency.")
    
    if (any(T.star < 0) || !is.numeric(T.star)) 
        stop("T.star must be numeric and may not contain negative numbers.")
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    
    n.bins = censor + 1
    transaction.actual = rep(0, n.bins)
    transaction.expected = rep(0, n.bins)
    bin.size = rep(0, n.bins)
    for (cc in 0:censor) {
        if (cc != censor) {
            this.bin = which(cc == x)
        } else if (cc == censor) {
            this.bin = which(x >= cc)
        }
        n.this.bin = length(this.bin)
        bin.size[cc + 1] = n.this.bin
        transaction.actual[cc + 1] = sum(x.star[this.bin])/n.this.bin
        transaction.expected[cc + 1] = sum(bgnbd.ConditionalExpectedTransactions(params, 
            T.star, x[this.bin], t.x[this.bin], T.cal[this.bin]))/n.this.bin
    }
    col.names = paste(rep("freq", length(censor + 1)), (0:censor), sep = ".")
    col.names[censor + 1] = paste(col.names[censor + 1], "+", sep = "")
    comparison = rbind(transaction.actual, transaction.expected, bin.size)
    colnames(comparison) = col.names
    if (is.null(xticklab) == FALSE) {
        x.labels = xticklab
    }
    if (is.null(xticklab) != FALSE) {
        if (censor < ncol(comparison)) {
            x.labels = 0:(censor)
            x.labels[censor + 1] = paste(censor, "+", sep = "")
        }
        if (censor >= ncol(comparison)) {
            x.labels = 0:(ncol(comparison))
        }
    }
    actual = comparison[1, ]
    expected = comparison[2, ]
    ylim = c(0, ceiling(max(c(actual, expected)) * 1.1))
    plot(actual, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, ylab = ylab, 
        main = title)
    lines(expected, lty = 2, col = 2)
    axis(1, at = 1:ncol(comparison), labels = x.labels)
    legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    return(comparison)
}

bgnbd.PlotRecVsConditionalExpectedFrequency <- function(params, cal.cbs, T.star, 
    x.star, xlab = "Calibration period recency", ylab = "Holdout period transactions", 
    xticklab = NULL, title = "Actual vs. Conditional Expected Transactions by Recency") {
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotRecVsConditionalExpectedFrequency")
    
    if (any(T.star < 0) || !is.numeric(T.star)) 
        stop("T.star must be numeric and may not contain negative numbers.")
    if (any(x.star < 0) || !is.numeric(x.star)) 
        stop("x.star must be numeric and may not contain negative numbers.")
    
    tryCatch(x <- cal.cbs[, "x"], error = function(e) stop("Error in bgnbd.PlotRecVsConditionalExpectedFrequency: cal.cbs must have a frequency column labelled \"x\""))
    tryCatch(t.x <- cal.cbs[, "t.x"], error = function(e) stop("Error in bgnbd.PlotRecVsConditionalExpectedFrequency: cal.cbs must have a recency column labelled \"t.x\""))
    tryCatch(T.cal <- cal.cbs[, "T.cal"], error = function(e) stop("Error in bgnbd.PlotRecVsConditionalExpectedFrequency: cal.cbs must have a column for length of time observed labelled \"T.cal\""))
    
    t.values <- sort(unique(t.x))
    n.recs <- length(t.values)
    transaction.actual <- rep(0, n.recs)
    transaction.expected <- rep(0, n.recs)
    rec.size <- rep(0, n.recs)
    
    for (tt in 1:n.recs) {
        this.t.x <- t.values[tt]
        this.rec <- which(t.x == this.t.x)
        n.this.rec <- length(this.rec)
        rec.size[tt] <- n.this.rec
        transaction.actual[tt] <- sum(x.star[this.rec])/n.this.rec
        transaction.expected[tt] <- sum(bgnbd.ConditionalExpectedTransactions(params, 
            T.star, x[this.rec], t.x[this.rec], T.cal[this.rec]))/n.this.rec
    }
    
    comparison <- rbind(transaction.actual, transaction.expected, rec.size)
    colnames(comparison) <- round(t.values, 3)
    
    bins <- seq(1, ceiling(max(t.x)))
    n.bins <- length(bins)
    actual <- rep(0, n.bins)
    expected <- rep(0, n.bins)
    bin.size <- rep(0, n.bins)
    
    x.labels <- NULL
    if (is.null(xticklab) == FALSE) {
        x.labels <- xticklab
    } else {
        x.labels <- 1:(n.bins)
    }
    point.labels <- rep("", n.bins)
    point.y.val <- rep(0, n.bins)
    for (ii in 1:n.bins) {
        if (ii < n.bins) {
            this.bin <- which(as.numeric(colnames(comparison)) >= (ii - 1) & as.numeric(colnames(comparison)) < 
                ii)
        } else if (ii == n.bins) {
            this.bin <- which(as.numeric(colnames(comparison)) >= ii - 1)
        }
        actual[ii] <- sum(comparison[1, this.bin])/length(comparison[1, this.bin])
        expected[ii] <- sum(comparison[2, this.bin])/length(comparison[2, this.bin])
        bin.size[ii] <- sum(comparison[3, this.bin])
    }
    
    ylim <- c(0, ceiling(max(c(actual, expected)) * 1.1))
    plot(actual, type = "l", xaxt = "n", col = 1, ylim = ylim, xlab = xlab, ylab = ylab, 
        main = title)
    lines(expected, lty = 2, col = 2)
    
    axis(1, at = 1:n.bins, labels = x.labels)
    legend("topleft", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    
    return(rbind(actual, expected, bin.size))
}

bgnbd.PlotTransactionRateHeterogeneity <- function(params, lim = NULL) {
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotTransactionRateHeterogeneity")
    shape <- params[1]
    rate <- params[2]
    rate.mean <- round(shape/rate, 4)
    rate.var <- round(shape/rate^2, 4)
    if (is.null(lim)) {
        lim = qgamma(0.99, shape = shape, rate = rate)
    }
    x.axis.ticks <- seq(0, lim, length.out = 100)
    heterogeneity <- dgamma(x.axis.ticks, shape = shape, rate = rate)
    plot(x.axis.ticks, heterogeneity, type = "l", xlab = "Transaction Rate", ylab = "Density", 
        main = "Heterogeneity in Transaction Rate")
    mean.var.label <- paste("Mean:", rate.mean, "    Var:", rate.var)
    mtext(mean.var.label, side = 3)
    return(rbind(x.axis.ticks, heterogeneity))
}

bgnbd.PlotDropoutRateHeterogeneity <- function(params, lim = NULL) {
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.PlotDropoutRateHeterogeneity")
    alpha_param = params[3]
    beta_param = params[4]
    beta_param.mean = round(alpha_param/(alpha_param + beta_param), 4)
    beta_param.var = round(alpha_param * beta_param/((alpha_param + beta_param)^2)/(alpha_param + 
        beta_param + 1), 4)
    if (is.null(lim)) {
        # get right end point of grid
        lim = qbeta(0.99, shape1 = alpha_param, shape2 = beta_param)
    }
    x.axis.ticks = seq(0, lim, length.out = 100)
    heterogeneity = dbeta(x.axis.ticks, shape1 = alpha_param, shape2 = beta_param)
    plot(x.axis.ticks, heterogeneity, type = "l", xlab = "Dropout Probability p", 
        ylab = "Density", main = "Heterogeneity in Dropout Probability")
    mean.var.label = paste("Mean:", beta_param.mean, "    Var:", beta_param.var)
    mtext(mean.var.label, side = 3)
    return(rbind(x.axis.ticks, heterogeneity))
}

bgnbd.ExpectedCumulativeTransactions <- function(params, T.cal, T.tot, n.periods.final) {
    
    dc.check.model.params(c("r", "alpha", "s", "beta"), params, "bgnbd.ExpectedCumulativeTransactions")
    
    if (any(T.cal < 0) || !is.numeric(T.cal)) 
        stop("T.cal must be numeric and may not contain negative numbers.")
    
    if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
        stop("T.cal must be a single numeric value and may not be negative.")
    if (length(n.periods.final) > 1 || n.periods.final < 0 || !is.numeric(n.periods.final)) 
        stop("n.periods.final must be a single numeric value and may not be negative.")
    
    intervals <- seq(T.tot/n.periods.final, T.tot, length.out = n.periods.final)
    
    cust.birth.periods <- max(T.cal) - T.cal
    
    expected.transactions <- sapply(intervals, function(interval) {
        if (interval <= min(cust.birth.periods)) 
            return(0)
        sum(bgnbd.Expectation(params, interval - cust.birth.periods[cust.birth.periods <= 
            interval]))
    })
    
    return(expected.transactions)
}


bgnbd.PlotTrackingCum <- function(params, T.cal, T.tot, actual.cu.tracking.data, 
    xlab = "Week", ylab = "Cumulative Transactions", xticklab = NULL, title = "Tracking Cumulative Transactions") {
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.Plot.PlotTrackingCum")
    
    if (any(T.cal < 0) || !is.numeric(T.cal)) 
        stop("T.cal must be numeric and may not contain negative numbers.")
    if (any(actual.cu.tracking.data < 0) || !is.numeric(actual.cu.tracking.data)) 
        stop("actual.cu.tracking.data must be numeric and may not contain negative numbers.")
    
    if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
        stop("T.cal must be a single numeric value and may not be negative.")
    
    actual <- actual.cu.tracking.data
    expected <- bgnbd.ExpectedCumulativeTransactions(params, T.cal, T.tot, length(actual))
    
    cu.tracking.comparison <- rbind(actual, expected)
    
    ylim <- c(0, max(c(actual, expected)) * 1.05)
    plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = ylim, 
        main = title)
    lines(expected, lty = 2, col = 2)
    if (is.null(xticklab) == FALSE) {
        if (ncol(cu.tracking.comparison) != length(xticklab)) {
            stop("Plot error, xticklab does not have the correct size")
        }
        axis(1, at = 1:ncol(cu.tracking.comparison), labels = xticklab)
    } else {
        axis(1, at = 1:length(actual), labels = 1:length(actual))
    }
    abline(v = max(T.cal), lty = 2)
    
    legend("bottomright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    
    return(cu.tracking.comparison)
}

bgnbd.PlotTrackingInc <- function(params, T.cal, T.tot, actual.inc.tracking.data, 
    xlab = "Week", ylab = "Transactions", xticklab = NULL, title = "Tracking Weekly Transactions") {
    
    dc.check.model.params(c("r", "alpha", "a", "b"), params, "bgnbd.Plot.PlotTrackingCum")
    
    if (any(T.cal < 0) || !is.numeric(T.cal)) 
        stop("T.cal must be numeric and may not contain negative numbers.")
    if (any(actual.inc.tracking.data < 0) || !is.numeric(actual.inc.tracking.data)) 
        stop("actual.inc.tracking.data must be numeric and may not contain negative numbers.")
    
    if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot)) 
        stop("T.cal must be a single numeric value and may not be negative.")
    
    actual <- actual.inc.tracking.data
    expected <- dc.CumulativeToIncremental(bgnbd.ExpectedCumulativeTransactions(params, 
        T.cal, T.tot, length(actual)))
    
    ylim <- c(0, max(c(actual, expected)) * 1.05)
    plot(actual, type = "l", xaxt = "n", xlab = xlab, ylab = ylab, col = 1, ylim = ylim, 
        main = title)
    lines(expected, lty = 2, col = 2)
    if (is.null(xticklab) == FALSE) {
        if (length(actual) != length(xticklab)) {
            stop("Plot error, xticklab does not have the correct size")
        }
        axis(1, at = 1:length(actual), labels = xticklab)
    } else {
        axis(1, at = 1:length(actual), labels = 1:length(actual))
    }
    abline(v = max(T.cal), lty = 2)
    
    legend("topright", legend = c("Actual", "Model"), col = 1:2, lty = 1:2, lwd = 1)
    
    return(rbind(actual, expected))
}


 
