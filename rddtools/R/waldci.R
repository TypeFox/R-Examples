#' Confint allowing vcov
#' 
#' Version of vcov allowing for confint
#' @param x Object of class lm or else
#' @param parm  specification of which parameters are to be given confidence intervals, see confint
#' @param level the confidence level required, see confint()
#' @param vcov. Specific covariance function to pass to coeftest. See help of sandwich
#' @param df Degrees of freedom
#' @param \ldots Further arguments


waldci <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, ...) {
    UseMethod("waldci")
}

waldci.default <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, ...) {
    ## use S4 methods if loaded
    coef0 <- if ("stats4" %in% loadedNamespaces()) 
        stats4::coef else coef
    vcov0 <- if ("stats4" %in% loadedNamespaces()) 
        stats4::vcov else vcov
    
    ## extract coefficients and standard errors
    est <- coef0(x)
    if (is.null(vcov.)) 
        se <- vcov0(x) else {
        if (is.function(vcov.)) 
            se <- vcov.(x) else se <- vcov.
    }
    se <- sqrt(diag(se))
    
    ## match using names and compute t/z statistics
    if (!is.null(names(est)) && !is.null(names(se))) {
        anames <- names(est)[names(est) %in% names(se)]
        est <- est[anames]
        se <- se[anames]
    }
    
    ## process level
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    
    ## get quantile from central limit theorem
    if (is.null(df)) {
        df <- try(df.residual(x), silent = TRUE)
        if (inherits(df, "try-error")) 
            df <- NULL
    }
    if (is.null(df)) 
        df <- 0
    fac <- if (is.finite(df) && df > 0) 
        qt(a, df = df) else qnorm(a)
    
    ## set up confidence intervals
    ci <- cbind(est + fac[1] * se, est + fac[2] * se)
    colnames(ci) <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3L), "%")
    
    ## process parm
    if (is.null(parm)) 
        parm <- seq_along(est)
    # if(is.character(parm)) parm <- which(parm %in% names(est))
    if (is.character(parm)) 
        parm <- which(names(est) %in% parm)
    ci <- ci[parm, , drop = FALSE]
    return(ci)
}


## copy of stats:::format.perc
format.perc <- function(probs, digits) paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")

waldci.rdd_reg_np <- function(x, level = 0.95, vcov. = NULL, df = Inf, ...) {
    
    inf_met <- infType(x)  ## def in Misc.R
    if (inf_met == "se") {
        if (!is.null(vcov.) | !is.infinite(df)) {
            warning("Arg 'vcov.' and 'df' only work for rdd_reg with inf='lm'")
        }
        ## code recycled from stats::confint.default
        co <- rdd_coef(x, allInfo = TRUE)
        a <- (1 - level)/2
        a <- c(a, 1 - a)
        fac <- qnorm(a)
        pct <- format.perc(a, 3)  ## import!!
        ci <- array(NA, dim = c(1, 2L), dimnames = list("D", pct))
        ci[] <- co[, "Estimate"] + co[, "Std. Error"] %o% fac
        return(ci)
    } else {
        waldci.default(x$RDDslot$model, parm = "D", level = level, vcov. = vcov., df = df, ...)
    }
}




waldci.glm <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = Inf, ...) waldci.default(x, parm = parm, level = level, 
    vcov. = vcov., df = df, ...)

waldci.mlm <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, ...) {
    ## obtain vcov
    v <- if (is.null(vcov.)) 
        vcov(x) else if (is.function(vcov.)) 
        vcov.(x) else vcov.
    
    ## nasty hack: replace coefficients so that their names match the vcov() method
    x$coefficients <- structure(as.vector(x$coefficients), .Names = colnames(vcov(x)))
    
    ## call default method
    waldci.default(x, parm = parm, level = level, vcov. = v, df = df, ...)
}

waldci.survreg <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = Inf, ...) {
    if (is.null(vcov.)) 
        v <- vcov(x) else {
        if (is.function(vcov.)) 
            v <- vcov.(x) else v <- vcov.
    }
    if (length(x$coefficients) < NROW(x$var)) {
        x$coefficients <- c(x$coefficients, `Log(scale)` = log(x$scale))
    }
    waldci.default(x, parm = parm, level = level, vcov. = v, df = df, ...)
}


if (FALSE) {
    
    library(sandwich)
    library(lmtest)
    
    reg <- lm(freeny)
    
    ### Regular
    all(confint(reg) == waldci(reg))
    confint(reg)
    co_reg <- coeftest(reg)
    co_reg[, 1] + qnorm(0.975) * co_reg[, 2]
    co_reg[, 1] + qt(0.975, df = reg[["df.residual"]]) * co_reg[, 2]
    
    ## vcovHC
    waldci(reg, vcov. = vcovHC)
    co <- coeftest(reg, vcov. = vcovHC)
    co[, 1] + qnorm(0.975) * co[, 2]
    co[, 1] + qt(0.975, df = reg[["df.residual"]]) * co[, 2]
    
} 
