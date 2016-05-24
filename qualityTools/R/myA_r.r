.myADTest = function(x, distribution, ...) {                                                                     ####   .MYADTESTS-FUNCTION
    #require(MASS, quietly = TRUE)
    if (missing(distribution)) 
        distribution = "normal"
    data.name = names(x)
    if (is.data.frame(x)) 
        x = x[, 1]
    dots = list(...)
    parameter = NULL
    smaller = NULL
    pFun = NULL
    tableValue = FALSE
    A = 0
    x <- sort(x[complete.cases(x)])
    n = length(x)
    if (n < 8) 
        stop("sample size must be greater than 7")
    if (n > 40) 
        warning("sample size is greater than 40")
    if (is.character(distribution)) {
        pFun = .charToDistFunc(distribution, type = "p")
        distribution = tolower(distribution)
        if (is.null(pFun)) 
            stop(paste(deparse(substitute(distribution)), " is not supported!"))
    }
    else {
        pFun = match.fun(distribution)
    }
#    if (identical(distribution, "log-normal")) {                               ####
#        x = log(x)                                                             ####
#        distribution = "normal"                                                ####
#    }                                                                          ####
    if (length(dots) == 0) {
        fittedDistr = MASS::fitdistr(x, distribution)
        parameter = fittedDistr$estimate
        if (distribution == "normal") {
            parameter["mean"] = mean(x)
            parameter["sd"] = sd(x)
        }
        p = do.call(pFun, c(list(x), as.list(parameter)))
    }
    else {
        p = pFun(x, ...)
    }
    h = (2 * seq(1:n) - 1) * (log(p) + log(1 - rev(p)))
    A = -n - mean(h)
    AA = (1 + 0.75/n + 2.25/n^2) * A
    if (AA < 0.2) {
        pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
    }
    else if (AA < 0.34) {
        pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
    }
    else if (AA < 0.6) {
        pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
    }
    else {
        pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
    }
    if (identical(distribution, "cauchy")) {
        pval = NA
    }
    if (identical(distribution, "beta")) {
        pval = NA
    }
    if (identical(distribution, "chi-squared")) {
        pval = NA
    }
    if (identical(distribution, "f")) {
        pval = NA
    }
    if (identical(distribution, "t")) {
        pval = NA
    }
    if (identical(distribution, "geometric")) {
        pval = NA
    }
    if (identical(distribution, "poisson")) {
        pval = NA
    }
    if (identical(distribution, "negative-binomial")) {
        pval = NA
    }
    if (identical(distribution, "weibull")) {
        AWei = A * (1 + 1/sqrt(n))
        tableValue = TRUE
        smaller = TRUE
        if (AWei < 0.474) {
            pval = 0.25
            smaller = FALSE
        }
        if (AWei >= 0.474) 
            pval = 0.25
        if (AWei >= 0.637) 
            pval = 0.1
        if (AWei >= 0.757) 
            pval = 0.05
        if (AWei >= 0.877) 
            pval = 0.025
        if (AWei >= 1.038) 
            pval = 0.01
    }
    if (identical(distribution, "exponential")) {
        AExp = A * (1 + 0.6/n)
        pval = NA
        if (0.95 < AExp) {
            pval = exp(0.731 - 3.009 * AExp + 0.15 * AExp^2)
        }
        if (0.51 < AExp & AExp < 0.95) {
            pval = exp(0.9209 - 3.353 * AExp + 0.3 * AExp^2)
        }
        if (0.26 < AExp & AExp < 0.51) {
            pval = 1 - exp(-6.1327 + 20.218 * AExp - 18.663 * AExp^2)
        }
        if (AExp < 0.26) {
            pval = 1 - exp(-12.2204 + 67.459 * AExp - 110.3 * AExp^2)
        }
    }
    if (identical(distribution, "logistic")) {
        ALogist = A * (1 + 0.25/n)
        tableValue = TRUE
        smaller = TRUE
        if (ALogist < 0.426) {
            pval = 0.25
            smaller = FALSE
        }
        if (ALogist >= 0.426) {
            pval = 0.25
        }
        if (ALogist >= 0.563) {
            pval = 0.1
        }
        if (ALogist >= 0.66) {
            pval = 0.05
        }
        if (ALogist >= 0.769) {
            pval = 0.025
        }
        if (ALogist >= 0.906) {
            pval = 0.01
        }
        if (ALogist >= 1.1) {
            pval = 0.005
        }
    }
    if (identical(distribution, "gamma")) {
        tableValue = TRUE
        gammaDF = data.frame(c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, Inf), c(0.486, 0.477, 0.475, 
            0.473, 0.472, 0.472, 0.471, 0.471, 0.471, 0.47, 0.47, 0.47), c(0.657, 0.643, 0.639, 0.637, 
            0.635, 0.635, 0.634, 0.633, 0.633, 0.632, 0.632, 0.631), c(0.786, 0.768, 0.762, 0.759, 
            0.758, 0.757, 0.755, 0.754, 0.754, 0.754, 0.753, 0.752), c(0.917, 0.894, 0.886, 0.883, 
            0.881, 0.88, 0.878, 0.877, 0.876, 0.876, 0.875, 0.873), c(1.092, 1.062, 1.052, 1.048, 
            1.045, 1.043, 1.041, 1.04, 1.039, 1.038, 1.037, 1.035), c(1.227, 1.19, 1.178, 1.173, 
            1.17, 1.168, 1.165, 1.164, 1.163, 1.162, 1.161, 1.159))
        names(gammaDF) = c("m", 0.75, 0.9, 0.95, 0.975, 0.99, 0.995)
        critCheck <- gammaDF[min(which(gammaDF$m >= parameter["shape"])), 2:length(gammaDF)] > A
        if (any(critCheck)) {
            firPos <- min(which(critCheck))
        }
        else {
            firPos <- length(gammaDF)
        }
        if (firPos == 1) {
            pValue <- 1 - as.numeric(names(gammaDF)[2])
            pval = pValue
            pValue <- paste(">", pValue)
            smaller = FALSE
        }
        else {
            pValue <- 1 - as.numeric(names(gammaDF)[firPos])
            pval = pValue
            pValue <- paste("<=", pValue)
            smaller = TRUE
        }
    }
    out = list()
    out$data.name = data.name
    out$statistic = as.vector(data.frame(A = A))
    out$parameter = parameter
    out$p.value = as.vector(data.frame(p = pval))
    out$smaller = smaller
    out$tableValue = tableValue
    out$conf.int = NULL
    out$estimate = NULL
    temp = NULL
    if (is.character(distribution)) 
        temp = as.vector(distribution)
    else temp = deparse(substitute(distribution))
    names(temp) = "distribution"
    out$null.value = temp
    out$method = paste("Anderson Darling Test for", temp, "distribution")
    class(out) = "adtest"
    return(out)
}
print.adtest = function(x, digits = 4, quote = TRUE, prefix = "", ...) {
    cat("\n")
    cat(strwrap(x$method, prefix = "\t"), sep = "\n")
    cat("\n")
    cat("data: ", x$data.name, "\n")
    out <- character()
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 4))))
    if (!is.null(x$parameter)) 
        out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 3))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = digits)
        if (x$tableValue) {
            if (x$smaller) 
                out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("<=", fp)))
            else out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "=") fp else paste(">", fp)))
        }
        else {
            out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
        }
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    cat("alternative hypothesis: ")
    if (!is.null(x$null.value)) {
        if (length(x$null.value) == 1) {
            cat("true", names(x$null.value), "is not equal to", x$null.value, "\n")
        }
        else {
            cat(x$alternative, "\nnull values:\n")
            print(x$null.value, ...)
        }
    }
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), "percent confidence interval:\n", format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
    }
    if (!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
} 
