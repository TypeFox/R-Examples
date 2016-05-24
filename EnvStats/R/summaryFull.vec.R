summaryFull.vec <-
function (x, stats = "all", trim = 0.1, sd.method = "sqrt.unbiased", 
    geo.sd.method = "sqrt.unbiased", skew.list = list(), kurtosis.list = list(), 
    cv.list = list(), digits = max(3, getOption("digits") - 3), 
    digit.type = "signif", x.name = NULL, stats.in.rows = TRUE) 
{
    digit.type <- match.arg(digit.type, c("signif", "round"))
    if (is.null(x.name)) 
        x.name <- deparse(substitute(x))
    x <- as.vector(unlist(x))
    if (length(x) < 1 || !is.numeric(x)) {
        x <- numeric(0)
        warning("length(x) < 1 and/or x not numeric.")
    }
    include.vec <- c(include.n = any(is.element(stats, c("all", 
        "for.non.pos", "n"))), include.n.miss = any(is.element(stats, 
        c("all", "for.non.pos", "n.miss"))), include.n.total = any(is.element(stats, 
        c("all", "for.non.pos", "n.total"))), include.mean = any(is.element(stats, 
        c("all", "for.non.pos", "mean"))), include.median = any(is.element(stats, 
        c("all", "for.non.pos", "median"))), include.trimmed.mean = any(is.element(stats, 
        c("all", "for.non.pos", "trimmed.mean"))), include.geo.mean = any(is.element(stats, 
        c("all", "geo.mean"))), include.skew = any(is.element(stats, 
        c("all", "for.non.pos", "skew"))), include.kurtosis = any(is.element(stats, 
        c("all", "for.non.pos", "kurtosis"))), include.min = any(is.element(stats, 
        c("all", "for.non.pos", "min"))), include.max = any(is.element(stats, 
        c("all", "for.non.pos", "max"))), include.range = any(is.element(stats, 
        c("all", "for.non.pos", "range"))), include.1st.quart = any(is.element(stats, 
        c("all", "for.non.pos", "1st.quart"))), include.3rd.quart = any(is.element(stats, 
        c("all", "for.non.pos", "3rd.quart"))), include.sd = any(is.element(stats, 
        c("all", "for.non.pos", "sd"))), include.geo.sd = any(is.element(stats, 
        c("all", "geo.sd"))), include.iqr = any(is.element(stats, 
        c("all", "for.non.pos", "iqr"))), include.mad = any(is.element(stats, 
        c("all", "for.non.pos", "mad"))), include.cv = any(is.element(stats, 
        c("all", "cv"))))
    n.stats <- sum(include.vec)
    if (n.stats == 0) 
        stop("No summary statistics specified.")
    n.total <- length(x)
    if (n.total == 0) {
        sum.vec <- c(n = 0, n.miss = 0, n.total = 0, rep(NA, 
            n.stats - 3))
        warning("'x' is an empty object")
    }
    else {
        wna <- which.na(x)
        n.miss <- length(wna)
        if (n.miss == n.total) {
            sum.vec <- c(n = 0, n.miss, n.total, rep(NA, n.stats - 
                3))
            warning("All values of 'x' are missing")
        }
        else {
            if (n.miss) 
                x <- x[-wna]
            if ((n.total - n.miss) == 1) {
                sum.vec <- c(n = 1, n.miss = n.miss, n.total = n.total, 
                  mean.x = x, median.x = x, tr.mean = x, geo.mean.x = if (include.vec["include.geo.mean"]) geoMean(x) else NA, 
                  skew.x = NA, kurtosis.x = NA, min.x = x, max.x = x, 
                  range.x = 0, Q1 = x, Q3 = x, sd.x = NA, geo.sd.x = NA, 
                  iqr = 0, mad.x = 0, cv.x = NA)
            }
            else {
                n.x <- length(x)
                s <- summary(x, digits = digits)
                names.s <- names(s)
                s <- as.vector(s)
                names(s) <- names.s
                mean.x <- s["Mean"]
                min.x <- s["Min."]
                max.x <- s["Max."]
                range.x <- max.x - min.x
                Q1 <- s["1st Qu."]
                Q3 <- s["3rd Qu."]
                iqr <- Q3 - Q1
                tr.mean <- mean(x, trim = trim)
                geo.mean.x <- if (include.vec["include.geo.mean"]) 
                  geoMean(x)
                else NA
                sd.x <- sd(x)
                if (sd.method != "sqrt.unbiased") 
                  sd.x <- sqrt((n.x - 1)/n.x) * sd.x
                geo.sd.x <- if (include.vec["include.geo.sd"]) 
                  geoSD(x, sqrt.unbiased = geo.sd.method == "sqrt.unbiased")
                else NA
                skew.x <- do.call("skewness", c(list(x = x), 
                  skew.list))
                kurtosis.x <- do.call("kurtosis", c(list(x = x), 
                  kurtosis.list))
                cv.x <- do.call("cv", c(list(x = x), cv.list))
                sum.vec <- c(n.x, n.miss, n.total, mean.x, s["Median"], 
                  tr.mean, geo.mean.x, skew.x, kurtosis.x, min.x, 
                  max.x, range.x, Q1, Q3, sd.x, geo.sd.x, iqr, 
                  mad(x), cv.x)
                sum.vec <- do.call(digit.type, list(sum.vec, 
                  digits = digits))
            }
        }
    }
    name.vec <- c("N", "NA's", "N.Total", "Mean", "Median", paste(trunc(trim * 
        100), "% Trimmed Mean", sep = ""), "Geometric Mean", 
        "Skew", "Kurtosis", "Min", "Max", "Range", "1st Quartile", 
        "3rd Quartile", "Standard Deviation", "Geometric Standard Deviation", 
        "Interquartile Range", "Median Absolute Deviation", "Coefficient of Variation")
    names(sum.vec) <- name.vec
    sum.vec <- sum.vec[include.vec]
    ret.val <- t(sum.vec)
    dimnames(ret.val)[[1]] <- x.name
    if (stats.in.rows) {
        ret.val <- t(ret.val)
    }
    ret.val
}
