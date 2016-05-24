print.intervalEstimate <-
function (x, conf.cov.sig.digits = .Options$digits, limits.sig.digits = .Options$digits, 
    ...) 
{
    ncxn <- nchar(x$name)
    if (!is.null(x$parameter)) {
        cat(x$name, " Interval for:", space(33 - ncxn - 14), 
            x$parameter, "\n\n", sep = "")
    }
    if (x$name == "Tolerance") {
        cat(x$name, " Interval Coverage:", space(33 - ncxn - 
            19), paste(format(100 * x$coverage, digits = conf.cov.sig.digits), 
            "%", sep = ""), "\n\n", sep = "")
        cat("Coverage Type:", space(19), x$coverage.type, "\n\n", 
            sep = "")
    }
    if (!is.null(x$method)) 
        cat(x$name, " Interval Method:", space(33 - ncxn - 17), 
            x$method, "\n\n", sep = "")
    if (!is.null(x$normal.transform.power)) 
        cat("Normal Transform Power:", space(10), x$normal.transform.power, 
            "\n\n", sep = "")
    cat(x$name, " Interval Type:", space(33 - ncxn - 15), x$type, 
        "\n\n", sep = "")
    if (!is.null(x$conf.level)) {
        cat("Confidence Level:", space(16), paste(format(100 * 
            x$conf.level, digits = conf.cov.sig.digits), "%", 
            sep = ""), "\n\n", sep = "")
    }
    if (!is.null(x$limit.ranks)) 
        cat(paste(x$name, " Limit Rank(s):", space(33 - ncxn - 
            16), sep = ""), x$limit.ranks, "\n\n", sep = " ")
    if (x$name == "Prediction") {
        rule <- x$rule
        standard.rule <- is.null(rule) || rule == "k.of.m"
        names.x <- names(x)
        if (is.element("r", names.x) && x$r > 1) {
            simultaneous <- TRUE
            string1 <- paste("\n(per Sampling Occasion):", space(9), 
                sep = "")
            string2 <- string1
        }
        else {
            simultaneous <- FALSE
            string1 <- paste(":", space(9), sep = "")
            string2 <- paste(":", space(13), sep = "")
        }
        if ((!is.null(x$n.sum) && x$n.sum == 1) || (!is.null(x$n.mean) && 
            x$n.mean == 1) || (!is.null(x$n.geomean) && x$n.geomean == 
            1) || (!is.null(x$n.transmean) && x$n.transmean == 
            1) || (!is.null(x$n.median) && x$n.median == 1)) {
            if (standard.rule) {
                if (is.element("m", names.x) && x$k < x$m) {
                  cat("Minimum Number of\nFuture Observations\n", 
                    "Interval Should Contain", string1, x$k, 
                    "\n\n", sep = "")
                  cat("Total Number of\nFuture Observations", 
                    string2, x$m, "\n\n", sep = "")
                }
                else {
                  cat("Number of Future Observations:", space(3), 
                    x$k, "\n\n", sep = "")
                }
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
            }
            else {
                cat("Maximum Number of\nFuture Observations", 
                  string2, x$m, "\n\n", sep = "")
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
            }
        }
        else {
            if (!is.null(x$n.sum)) {
                if (is.element("m", names.x) && x$k < x$m) {
                  cat("Minimum Number of\nFuture Sums\n", "Interval Should Contain", 
                    string1, x$k, "\n\n", sep = "")
                  if (!simultaneous) 
                    string2 <- paste(":", space(21), sep = "")
                  cat("Total Number of\nFuture Sums", string2, 
                    x$m, "\n\n", sep = "")
                }
                else {
                  cat("Number of Future Sums:", space(11), x$k, 
                    "\n\n", sep = "")
                }
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
                cat("Sample Size for Sums:", space(12), x$n.sum, 
                  "\n\n", sep = "")
            }
            else if (!is.null(x$n.mean)) {
                if (is.element("m", names.x) && x$k < x$m) {
                  cat("Minimum Number of\nFuture Means\n", "Interval Should Contain", 
                    string1, x$k, "\n\n", sep = "")
                  if (!simultaneous) 
                    string2 <- paste(":", space(20), sep = "")
                  cat("Total Number of\nFuture Means", string2, 
                    x$m, "\n\n", sep = "")
                }
                else {
                  cat("Number of Future Means:", space(10), x$k, 
                    "\n\n", sep = "")
                }
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
                cat("Sample Size for Means:", space(11), x$n.mean, 
                  "\n\n", sep = "")
            }
            else if (!is.null(x$n.geomean)) {
                if (is.element("m", names.x) && x$k < x$m) {
                  cat("Minimum Number of\nFuture Geometric Means\n", 
                    "Interval Should Contain", string1, x$k, 
                    "\n\n", sep = "")
                  if (!simultaneous) 
                    string2 <- paste(":", space(10), sep = "")
                  cat("Total Number of\nFuture Geometric Means", 
                    string2, x$m, "\n\n", sep = "")
                }
                else {
                  cat("Number of Future\nGeometric Means:", space(17), 
                    x$k, "\n\n", sep = "")
                }
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
                cat("Sample Size for\nGeometric Means:", space(17), 
                  x$n.geomean, "\n\n", sep = "")
            }
            else if (!is.null(x$n.transmean)) {
                if (is.element("m", names.x) && x$k < x$m) {
                  cat("Minimum Number of\nFuture Transformed Means\n", 
                    "Interval Should Contain", string1, x$k, 
                    "\n\n", sep = "")
                  if (!simultaneous) 
                    string2 <- paste(":", space(8), sep = "")
                  cat("Total Number of\nFuture Transformed Means", 
                    string2, x$m, "\n\n", sep = "")
                }
                else {
                  cat("Number of Future\nTransformed Means:", 
                    space(15), x$k, "\n\n", sep = "")
                }
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
                cat("Sample Size for\nTransformed Means:", space(15), 
                  x$n.transmean, "\n\n", sep = "")
            }
            else if (!is.null(x$n.median)) {
                if (is.element("m", names.x) && x$k < x$m) {
                  cat("Minimum Number of\nFuture Medians\n", 
                    "Interval Should Contain", string1, x$k, 
                    "\n\n", sep = "")
                  if (!simultaneous) 
                    string2 <- paste(":", space(18), sep = "")
                  cat("Total Number of\nFuture Medians", string2, 
                    x$m, "\n\n", sep = "")
                }
                else {
                  cat("Number of Future Medians:", space(5), 
                    x$k, "\n\n", sep = "")
                }
                if (simultaneous) {
                  cat("Number of Future\nSampling Occasions:", 
                    space(14), x$r, "\n\n", sep = "")
                }
                cat("Sample Size for Medians:", space(9), x$n.median, 
                  "\n\n", sep = "")
            }
        }
    }
    cat(x$name, " Interval:", space(33 - ncxn - 10), paste(paste(format(names(x$limits), 
        justify = "left"), format(x$limits, digits = limits.sig.digits, 
        nsmall = 0), sep = " = "), collapse = paste("\n", space(33), 
        sep = "")), "\n\n", sep = "")
    invisible(x)
}
