### Taha test
taha_test <- function(object, ...) UseMethod("taha_test")

taha_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("taha_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

taha_test.IndependenceProblem <- function(object,
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                rank_trafo(y)^2),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            if (is_ordered_x(object))
                stop(sQuote(colnames(object@x)), " is an ordered factor")
            return(TRUE)
        }
    )
    ## set test statistic to scalar for two-sample test
    args$teststat <- if (twosamp) "scalar" else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        object@method <- "Two-Sample Taha Test"
        object@parameter <- "ratio of scales"
        object@nullvalue <- 1
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_scale(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Taha Test"

    return(object)
}


### Klotz Test
klotz_test <- function(object, ...) UseMethod("klotz_test")

klotz_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("klotz_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

klotz_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                klotz_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            if (is_ordered_x(object))
                stop(sQuote(colnames(object@x)), " is an ordered factor")
            return(TRUE)
        }
    )
    ## set test statistic to scalar for two-sample test
    args$teststat <- if (twosamp) "scalar" else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        object@method <- "Two-Sample Klotz Test"
        object@parameter <- "ratio of scales"
        object@nullvalue <- 1
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_scale(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Klotz Test"

    return(object)
}


### Mood Test
mood_test <- function(object, ...) UseMethod("mood_test")

mood_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("mood_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

mood_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                mood_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            if (is_ordered_x(object))
                stop(sQuote(colnames(object@x)), " is an ordered factor")
            return(TRUE)
        }
    )
    ## set test statistic to scalar for two-sample test
    args$teststat <- if (twosamp) "scalar" else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        object@method <- "Two-Sample Mood Test"
        object@parameter <- "ratio of scales"
        object@nullvalue <- 1
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_scale(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Mood Test"

    return(object)
}


### Ansari-Bradley test
ansari_test <- function(object, ...) UseMethod("ansari_test")

ansari_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("ansari_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

ansari_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                ansari_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (!is_numeric_y(object))
                stop(sQuote(colnames(object@y)), " is not a numeric variable")
            if (is_ordered_x(object))
                stop(sQuote(colnames(object@x)), " is an ordered factor")
            return(TRUE)
        }
    )
    ## set test statistic to scalar for two-sample test
    args$teststat <- if (twosamp) "scalar" else "quadratic"
    ## swap alternative in one-sample case
    ## (a *large* test statistic implies that sample 1 is *less* variable)
    if (twosamp) {
        alternative <- match.arg(args$alternative,
                                 c("two.sided", "less", "greater"))
        if (alternative == "less")
            args$alternative <- "greater"
        else if (alternative == "greater")
            args$alternative <- "less"
    }

    object <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        object@method <- "Two-Sample Ansari-Bradley Test"
        object@parameter <- "ratio of scales"
        object@nullvalue <- 1
        object@statistic@alternative <- alternative
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_scale(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Ansari-Bradley Test"

    return(object)
}


### Fligner-Killeen test
fligner_test <- function(object, ...) UseMethod("fligner_test")

fligner_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("fligner_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

fligner_test.IndependenceProblem <- function(object,
    ties.method = c("mid-ranks", "average-scores"),
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                fligner_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (is_ordered_x(object))
                stop(sQuote(colnames(object@x)), " is an ordered factor")
            return(TRUE)
        }
    )
    ## set test statistic to scalar for two-sample test
    args$teststat <- if (twosamp) "scalar" else "quadratic"

    ## eliminate location differences (see 'stats/R/fligner.test')
    if (!is_numeric_y(object))
        stop(sQuote(colnames(object@y)), " is not a numeric variable")
    object@y[[1]] <- object@y[[1]] -
        tapply(object@y[[1]], object@x[[1]], median)[object@x[[1]]]

    object <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        object@method <- "Two-Sample Fligner-Killeen Test"
        object@parameter <- "ratio of scales"
        object@nullvalue <- 1
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_scale(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Fligner-Killeen Test"

    return(object)
}


### Conover-Iman (1978) test
conover_test <- function(object, ...) UseMethod("conover_test")

conover_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("conover_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

conover_test.IndependenceProblem <- function(object,
    conf.int = FALSE, conf.level = 0.95, ...) {

    twosamp <- is_2sample(object)

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                rank_trafo(abs(y))^2),
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
            if (is_ordered_x(object))
                stop(sQuote(colnames(object@x)), " is an ordered factor")
            return(TRUE)
        }
    )
    ## set test statistic to scalar for two-sample test
    args$teststat <- if (twosamp) "scalar" else "quadratic"

    ## eliminate location differences
    if (!is_numeric_y(object))
        stop(sQuote(colnames(object@y)), " is not a numeric variable")
    object@y[[1]] <- object@y[[1]] -
        tapply(object@y[[1]], object@x[[1]], mean)[object@x[[1]]]

    object <- do.call("independence_test", c(list(object = object), args))

    if (twosamp) {
        object@method <- "Two-Sample Conover-Iman Test"
        object@parameter <- "ratio of scales"
        object@nullvalue <- 1
        if (conf.int && has_distribution(args)) {
            object <- new("ScalarIndependenceTestConfint", object)
            object@confint <- function(level)
                confint_scale(object@statistic, object@distribution,
                                 level = level)
            object@conf.level <- conf.level
        }
    } else
        object@method <- "K-Sample Conover-Iman Test"

    return(object)
}
