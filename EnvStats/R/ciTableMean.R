ciTableMean <-
function (n1 = 10, n2 = n1, diff.or.mean = 2:0, SD = 1:3, sample.type = "two.sample", 
    ci.type = "two.sided", conf.level = 0.95, digits = 1) 
{
    sample.type <- match.arg(sample.type, c("two.sample", "one.sample"))
    ci.type <- match.arg(ci.type, c("two.sided", "lower", "upper"))
    if (!is.vector(n1, mode = "numeric") || length(n1) != 1 || 
        !is.finite(n1) || n1 != trunc(n1) || n1 < 2) 
        stop("\"n1\" must be a positive integer greater than 1")
    if (sample.type == "two.sample") {
        if (!is.vector(n2, mode = "numeric") || length(n2) != 
            1 || !is.finite(n2) || n2 != trunc(n2) || n2 < 2) 
            stop("\"n2\" must be a positive integer greater than 1")
    }
    if (!is.vector(diff.or.mean, mode = "numeric") || length(diff.or.mean) < 
        1 || any(!is.finite(diff.or.mean))) 
        stop(paste("\"diff.or.mean\" must be a numeric vector with at least one element", 
            "and no missing (NA), undefined (NaN), or non-finite (-Inf, Inf) values"))
    if (!is.vector(SD, mode = "numeric") || length(SD) < 1 || 
        any(!is.finite(SD))) 
        stop(paste("\"SD\" must be a numeric vector with at least one element", 
            "and no missing (NA), undefined (NaN), or non-finite (-Inf, Inf) values"))
    if (length(conf.level) != 1 || !is.numeric(conf.level) || 
        conf.level <= 0 || conf.level >= 100) 
        stop("\"conf.level\" must be a numeric scalar greater than 0 and less than 100")
    if (length(digits) != 1 || !is.numeric(digits) || digits != 
        trunc(digits) || digits < 0) 
        stop("\"digits\" must be a positive integer")
    if (ci.type == "two.sided") {
        if (sample.type == "one.sample") {
            hw <- qt(p = (1 + conf.level)/2, df = n1 - 1) * SD/sqrt(n1)
        }
        else {
            hw <- qt(p = (1 + conf.level)/2, df = n1 + n2 - 2) * 
                SD * sqrt(1/n1 + 1/n2)
        }
        LCL <- format(round(t(outer(diff.or.mean, hw, FUN = "-")), 
            digits = digits))
        UCL <- format(round(t(outer(diff.or.mean, hw, FUN = "+")), 
            digits = digits))
    }
    else if (ci.type == "lower") {
        if (sample.type == "one.sample") {
            hw <- qt(p = conf.level, df = n1 - 1) * SD/sqrt(n1)
        }
        else {
            hw <- qt(p = conf.level, df = n1 + n2 - 2) * SD * 
                sqrt(1/n1 + 1/n2)
        }
        LCL <- format(round(t(outer(diff.or.mean, hw, FUN = "-")), 
            digits = digits))
        UCL <- rep(Inf, length(LCL))
    }
    else {
        if (sample.type == "one.sample") {
            hw <- qt(p = conf.level, df = n1 - 1) * SD/sqrt(n1)
        }
        else {
            hw <- qt(p = conf.level, df = n1 + n2 - 2) * SD * 
                sqrt(1/n1 + 1/n2)
        }
        UCL <- format(round(t(outer(diff.or.mean, hw, FUN = "+")), 
            digits = digits))
        LCL <- rep(-Inf, length(UCL))
    }
    df <- data.frame(matrix(paste("[", LCL, ", ", UCL, "]", sep = ""), 
        nrow = length(SD), ncol = length(diff.or.mean)), stringsAsFactors = FALSE)
    row.names(df) <- paste("SD=", SD, sep = "")
    if (sample.type == "two.sample") 
        names(df) <- paste("Diff=", diff.or.mean, sep = "")
    else names(df) <- paste("Mean=", diff.or.mean, sep = "")
    df
}
