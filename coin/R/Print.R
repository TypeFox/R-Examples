setMethod("show",
    signature = "IndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(
            statistic = setNames(object@statistic@teststatistic, nm = "c"),
            p.value = object@distribution@pvalue(object@statistic@teststatistic),
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "MaxTypeIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(
            statistic = setNames(object@statistic@teststatistic, nm = "maxT"),
            p.value = object@distribution@pvalue(object@statistic@teststatistic),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "QuadTypeIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(
            statistic = setNames(object@statistic@teststatistic, nm = "chi-squared"),
            p.value = object@distribution@pvalue(object@statistic@teststatistic),
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        if (length(object@distribution@parameters) == 1 &&
              names(object@distribution@parameters) == "df")
            RET$parameter <- setNames(object@distribution@parameters[[1]], nm = "df")
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTest",
    definition = function(object) {
        distname <- switch(class(object@distribution),
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")

        RET <- list(
            statistic = setNames(object@statistic@teststatistic, nm = "Z"),
            p.value = object@distribution@pvalue(object@statistic@teststatistic),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method)
        )
        if (length(object@nullvalue) > 0)
            RET$null.value = setNames(object@nullvalue, nm = object@parameter)
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

setMethod("show",
    signature = "ScalarIndependenceTestConfint",
    definition = function(object) {
        distname <- switch(class(object@distribution),
            "AsymptNullDistribution" = "Asymptotic",
            "ApproxNullDistribution" = "Approximative",
            "ExactNullDistribution" = "Exact")
        ci <- confint(object, level = object@conf.level)

        RET <- list(
            statistic = setNames(object@statistic@teststatistic, nm = "Z"),
            p.value = object@distribution@pvalue(object@statistic@teststatistic),
            alternative = object@statistic@alternative,
            data.name = varnames(object@statistic),
            method = paste(distname, object@method),
            conf.int = ci$conf.int,
            estimate = ci$estimate
        )
        if (length(object@nullvalue) > 0)
            RET$null.value = setNames(object@nullvalue, nm = object@parameter)
        if (length(object@estimates) > 0)
            RET <- c(RET, object@estimates)
        class(RET) <- "htest"
        print(RET)
        invisible(RET)
    }
)

print.ci <- function(x, ...) {
    if(!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")),
            "percent confidence interval:\n",
            format(c(x$conf.int[1], x$conf.int[2])), "\n")
    }
    if(!is.null(x$estimate)) {
        cat("sample estimates:\n")
        print(x$estimate, ...)
    }
    cat("\n")
    invisible(x)
}

print.MCp <- function(x, ...) {
    p <- x
    attributes(p) <- NULL
    print(p)
    ci <- list(conf.int = attr(x, "conf.int"))
    class(ci) <- "ci"
    print(ci)
}

print.cutpoint <- function(x, ...) {
    cat(paste0("  ", dQuote("best"), " cutpoint: ", x$label, "\n"))
    if (!is.null(x$covariable))
        cat(paste0("       covariable: ", x$covariable, "\n"))
}
