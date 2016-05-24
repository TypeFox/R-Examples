ciTableProp <-
function (n1 = 10, p1.hat = c(0.1, 0.2, 0.3), n2 = n1, p2.hat.minus.p1.hat = c(0.2, 
    0.1, 0), sample.type = "two.sample", ci.type = "two.sided", 
    conf.level = 0.95, digits = 2, ci.method = "score", correct = TRUE, 
    tol = 10^-(digits + 1)) 
{
    sample.type <- match.arg(sample.type, c("two.sample", "one.sample"))
    ci.type <- match.arg(ci.type, c("two.sided", "lower", "upper"))
    alternative <- switch(ci.type, two.sided = "two.sided", lower = "greater", 
        upper = "less")
    ci.method <- match.arg(ci.method, c("score", "exact"))
    if (sample.type == "two.sample" & ci.method != "score") 
        stop("Exact method not available for two-sample confidence interval")
    if (!is.vector(n1, mode = "numeric") || length(n1) != 1 || 
        !is.finite(n1) || n1 != trunc(n1) || n1 < 2) 
        stop("\"n1\" must be a positive integer greater than 1")
    if (!is.vector(p1.hat, mode = "numeric") || length(p1.hat) < 
        1 || any(!is.finite(p1.hat)) || any(p1.hat <= 0) || any(p1.hat >= 
        1)) 
        stop(paste("\"p1.hat\" must be a numeric vector with at least one element,", 
            "and all elements must be greater than 0 and less than 1"))
    if (length(conf.level) != 1 || !is.numeric(conf.level) || 
        conf.level <= 0 || conf.level >= 100) 
        stop("\"conf.level\" must be a numeric scalar greater than 0 and less than 100")
    if (length(digits) != 1 || !is.numeric(digits) || digits != 
        trunc(digits) || digits < 0) 
        stop("\"digits\" must be a positive integer")
    if (sample.type == "two.sample") {
        if (!is.vector(n2, mode = "numeric") || length(n2) != 
            1 || !is.finite(n2) || n2 != trunc(n2) || n2 < 2) 
            stop("\"n2\" must be a positive integer greater than 1")
        if (!is.vector(n2, mode = "numeric") || length(p2.hat.minus.p1.hat) < 
            1 || any(!is.finite(p2.hat.minus.p1.hat))) 
            stop(paste("\"p2.hat.minus.p1.hat\" must be a numeric vector with at least one element,", 
                "and no missing (NA), undefined (NaN), or non-finite (-Inf, Inf) values"))
        x1 <- unique(round((p1.hat * n1), 0))
        p1.hat <- x1/n1
        n.row <- length(p1.hat)
        n.col <- length(p2.hat.minus.p1.hat)
        x1.rep <- rep(x1, n.col)
        p2.hat.rep <- as.vector(outer(p1.hat, p2.hat.minus.p1.hat, 
            FUN = "+"))
        x2.rep <- round((p2.hat.rep * n2), 0)
        p2.hat.rep <- x2.rep/n2
        if (any(p2.hat.rep <= 0) || any(p2.hat.rep >= 1)) 
            stop(paste("One or more values of \"p2.hat.minus.p1.hat\"", 
                "is/are too small or too big", "for the given value(s) of \"p1.hat\""))
        n <- length(x1.rep)
        LCL <- rep(as.numeric(NA), n)
        UCL <- rep(as.numeric(NA), n)
        Diff <- rep(as.numeric(NA), n)
        o.warn <- options(warn = -1)
        for (i in 1:n) {
            test.list <- prop.test(x = c(x1.rep[i], x2.rep[i]), 
                n = c(n1, n2), alternative = alternative, conf.level = conf.level, 
                correct = correct)
            LCL[i] <- -test.list$conf.int[2]
            UCL[i] <- -test.list$conf.int[1]
            Diff[i] <- diff(test.list$estimate)
        }
        options(o.warn)
        LCL <- format(round(LCL, digits = digits))
        UCL <- format(round(UCL, digits = digits))
        df <- data.frame(matrix(paste("[", LCL, ", ", UCL, "]", 
            sep = ""), nrow = n.row, ncol = n.col), stringsAsFactors = FALSE)
        row.names(df) <- paste("P1.hat=", format(round(p1.hat, 
            digits)), sep = "")
        Diff.mat <- matrix(Diff, nrow = n.row, ncol = n.col)
        if (all(abs(diff(Diff.mat)) < tol)) {
            names(df) <- paste("Diff=", round(Diff.mat[1, ], 
                digits), sep = "")
        }
        else {
            df <- data.frame(df, round(Diff.mat, digits), stringsAsFactors = FALSE)
            col.index <- as.vector(t(matrix(c((1:n.col) + n.col, 
                1:n.col), ncol = 2)))
            df <- df[, col.index]
            names(df) <- rep(c("Diff", "CI"), n.col)
        }
    }
    else {
        x <- unique(round((p1.hat * n1), 0))
        p.hat <- x/n1
        n <- length(x)
        LCL <- rep(as.numeric(NA), n)
        UCL <- rep(as.numeric(NA), n)
        if (ci.method == "score") {
            o.warn <- options(warn = -1)
            for (i in 1:n) {
                test.list <- prop.test(x = x[i], n = n1, alternative = alternative, 
                  conf.level = conf.level, correct = correct)
                LCL[i] <- test.list$conf.int[1]
                UCL[i] <- test.list$conf.int[2]
            }
            options(o.warn)
            LCL <- format(round(LCL, digits = digits))
            UCL <- format(round(UCL, digits = digits))
        }
        else {
            for (i in 1:n) {
                test.list <- binom.test(x = x[i], n = n1, alternative = alternative, 
                  conf.level = conf.level)
                LCL[i] <- test.list$conf.int[1]
                UCL[i] <- test.list$conf.int[2]
            }
            LCL <- format(round(LCL, digits = digits))
            UCL <- format(round(UCL, digits = digits))
        }
        df <- data.frame(matrix(paste("[", LCL, ", ", UCL, "]", 
            sep = ""), nrow = 1, ncol = length(p.hat)), stringsAsFactors = FALSE)
        row.names(df) <- ""
        names(df) <- paste("P.hat=", p.hat, sep = "")
    }
    df
}
