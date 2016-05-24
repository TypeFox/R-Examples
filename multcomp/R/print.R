
# $Id: print.R 299 2010-08-25 14:25:31Z thothorn $

### print methods
print.glht <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n\t", "General Linear Hypotheses\n\n")
    if (!is.null(x$type))
        cat("Multiple Comparisons of Means:", x$type, "Contrasts\n\n\n")
    beta <- coef(x)
    lh <- matrix(beta, ncol = 1)
    colnames(lh) <- "Estimate"
    alt <- switch(x$alternative,
                  "two.sided" = "==", "less" = ">=", "greater" = "<=")
    rownames(lh) <- paste(names(beta), alt, x$rhs)
    cat("Linear Hypotheses:\n")
    print(lh, digits = digits)
    cat("\n")
    invisible(x)
}

print.summary.glht <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
    if (!is.null(x$type))
        cat("Multiple Comparisons of Means:", x$type, "Contrasts\n\n\n")
    call <- if (isS4(x$model)) x$model@call else x$model$call
    if (!is.null(call)) {
        cat("Fit: ")
        print(call)
        cat("\n")
    }

    ### <FIXME>: compute coefmat in summary.glht for easier access???
    pq <- x$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    error <- attr(pq$pvalues, "error")
    pname <- switch(x$alternativ,
        "less" = paste("Pr(<", ifelse(x$df == 0, "z", "t"), ")", sep = ""),
        "greater" = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""),
        "two.sided" = paste("Pr(>|", ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
    colnames(mtests) <- c("Estimate", "Std. Error",
        ifelse(x$df == 0, "z value", "t value"), pname)
    type <- pq$type

    ### print p values according to simulation precision
    if (!is.null(error) && error > .Machine$double.eps) {
        sig <- which.min(abs(1 / error - (10^(1:10))))
        sig <- 1 / (10^sig)
    } else {
        sig <- .Machine$double.eps
    }
    cat("Linear Hypotheses:\n")
    alt <- switch(x$alternative,
                  "two.sided" = "==", "less" = ">=", "greater" = "<=")
    rownames(mtests) <- paste(rownames(mtests), alt, x$rhs)
    ### </FIXME>

    printCoefmat(mtests, digits = digits, 
                 has.Pvalue = TRUE, P.values = TRUE, eps.Pvalue = sig)

    switch(type, 
        "univariate" = cat("(Univariate p values reported)"),
        "single-step" = cat("(Adjusted p values reported -- single-step method)"),
        "Shaffer" = cat("(Adjusted p values reported -- Shaffer method)"),
        "Westfall" = cat("(Adjusted p values reported -- Westfall method)"),
        cat("(Adjusted p values reported --", type, "method)")
    )
    cat("\n\n")
    invisible(x)                    
}

print.confint.glht <- function(x, digits = max(3, getOption("digits") - 3), 
                              ...) 
{
    xtmp <- x
    cat("\n\t", "Simultaneous Confidence Intervals\n\n")
    if (!is.null(x$type))
        cat("Multiple Comparisons of Means:", x$type, "Contrasts\n\n\n")
    level <- attr(x$confint, "conf.level")
    attr(x$confint, "conf.level") <- NULL
    cat("Fit: ")
    if (isS4(x$model)) {
        print(x$model@call)
    } else {
        print(x$model$call)
    }
    cat("\n")
    error <- attr(x$confint, "error")
    if (!is.null(error) && error > .Machine$double.eps)
        digits <- min(digits, which.min(abs(1 / error - (10^(1:10)))))
    cat("Quantile =", round(attr(x$confint, "calpha"), digits))
    cat("\n")
    if (attr(x, "type") == "adjusted") {
        cat(paste(level * 100, 
                  "% family-wise confidence level\n", sep = ""), "\n\n")
    } else {
        cat(paste(level * 100, 
                  "% confidence level\n", sep = ""), "\n\n")
    }
    cat("Linear Hypotheses:\n")
    alt <- switch(x$alternative,
                  "two.sided" = "==", "less" = ">=", "greater" = "<=")
    rownames(x$confint) <- paste(rownames(x$confint), alt, x$rhs)
    print(format(x$confint, nsmall = digits, digits = digits), quote = FALSE)
    cat("\n")
    invisible(xtmp)
}

print.contrMat <- function(x, digits = max(3, getOption("digits") - 3), ...) {

    xtmp <- x
    cat("\n\t", "Multiple Comparisons of Means:", attr(x, "type"), "Contrasts\n\n")
    attr(x, "type") <- NULL
    class(x) <- "matrix"  
    print(x, digits = digits)
    invisible(xtmp)
}

print.summary.gtest <- function(x, 
    digits = max(3, getOption("digits") - 3), ...) {

    print.glht(x, digits = digits)
    cat("Global Test:\n")
    if (x$test$type == "Chisq") {
        pr <- data.frame(x$test$SSH, x$test$df[1], x$test$pval)
        names(pr) <- c("Chisq", "DF", "Pr(>Chisq)")
    }
    if (x$test$type == "F") {
        pr <- data.frame(x$test$fstat, x$test$df[1], x$test$df[2], 
                         x$test$pval)
        names(pr) <- c("F", "DF1", "DF2", "Pr(>F)")
    }
    print(pr, digits = digits)
    invisible(x)
}
