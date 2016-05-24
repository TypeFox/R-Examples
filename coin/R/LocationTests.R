### permutation test without transformations
oneway_test <- function(object, ...) UseMethod("oneway_test")

oneway_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("oneway_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

oneway_test.IndependenceProblem <- function(object, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            return(TRUE)
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear and two-sample tests
    args$teststat <- if (is_ordered_x(object) || twosamp) "scalar"
                     else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        object@method <- "Two-Sample Fisher-Pitman Permutation Test"
        object@nullvalue <- 0
    } else
        object@method <- "K-Sample Fisher-Pitman Permutation Test"

    return(object)
}


### OK, OK, here is the most prominent one ...
wilcox_test <- function(object, ...) UseMethod("wilcox_test")

wilcox_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("wilcox_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

wilcox_test.IndependenceProblem <- function(object,
    conf.int = FALSE, conf.level = 0.95, ...) {

    args <- setup_args(
        teststat = "scalar",
        ytrafo = function(data)
            trafo(data, numeric_trafo = rank_trafo),
        check = function(object) {
            if (!is_2sample(object))
                stop(sQuote("object"),
                     " does not represent a two-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            return(TRUE)
        }
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Wilcoxon-Mann-Whitney Test"
    object@nullvalue <- 0

    if (conf.int && has_distribution(args)) {
        object <- new("ScalarIndependenceTestConfint", object)
        object@confint <- function(level)
            confint_location(object@statistic, object@distribution,
                             level = level)
        object@conf.level <- conf.level
    }

    return(object)
}


### Kruskal-Wallis test
kruskal_test <- function(object, ...) UseMethod("kruskal_test")

kruskal_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("kruskal_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

kruskal_test.IndependenceProblem <- function(object, ...) {

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = rank_trafo),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            return(TRUE)
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <- if (is_ordered_x(object)) "scalar"
                     else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else
        object@method <- "Kruskal-Wallis Test"

    return(object)
}


### normal quantiles (van der Waerden) test
normal_test <- function(object, ...) UseMethod("normal_test")

normal_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("normal_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

normal_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                normal_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            return(TRUE)
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear and two-sample tests
    args$teststat <- if (is_ordered_x(object) || twosamp) "scalar"
                     else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        object@method <- "Two-Sample van der Waerden (Normal Quantile) Test"
        object@nullvalue <- 0
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_location(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample van der Waerden (Normal Quantile) Test"

    return(object)
}


### median test
median_test <- function(object, ...) UseMethod("median_test")

median_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("median_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

median_test.IndependenceProblem <- function(object,
    mid.score = c("0", "0.5", "1"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                median_trafo(y, mid.score = mid.score)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            return(TRUE)
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear and two-sample tests
    args$teststat <- if (is_ordered_x(object) || twosamp) "scalar"
                     else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        object@method <- "Two-Sample Brown-Mood Median Test"
        object@nullvalue <- 0
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_location(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Brown-Mood Median Test"

    return(object)
}


### Savage test
savage_test <- function(object, ...) UseMethod("savage_test")

savage_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("savage_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

savage_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                savage_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            return(TRUE)
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear and two-sample tests
    args$teststat <- if (is_ordered_x(object) || twosamp) "scalar"
                     else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (twosamp) {
        object@method <- "Two-Sample Savage Test"
        object@nullvalue <- 0
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_location(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Savage Test"

    return(object)
}
