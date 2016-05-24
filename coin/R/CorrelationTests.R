### Spearman test
spearman_test <- function(object, ...) UseMethod("spearman_test")

spearman_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("spearman_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

spearman_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate", "none"), ...) {

    args <- setup_args(
        teststat = "scalar",
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate", "none")
        ),
        xtrafo = function(data)
            trafo(data, numeric_trafo = rank_trafo),
        ytrafo = function(data)
            trafo(data, numeric_trafo = rank_trafo),
        check = function(object) {
            if (!is_corr(object))
                stop(sQuote("object"),
                     " does not represent a univariate correlation problem")
            return(TRUE)
        }
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Spearman Correlation Test"
    object@parameter <- "rho"
    object@nullvalue <- 0

    return(object)
}


### Fisher-Yates correlation test (based on van der Waerden scores)
fisyat_test <- function(object, ...) UseMethod("fisyat_test")

fisyat_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("fisyat_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

fisyat_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate", "none"),
    ties.method = c("mid-ranks", "average-scores"), ...) {

    args <- setup_args(
        teststat = "scalar",
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate", "none")
        ),
        xtrafo = function(data)
            trafo(data, numeric_trafo = function(x)
                normal_trafo(x, ties.method = ties.method)),
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                normal_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_corr(object))
                stop(sQuote("object"),
                     " does not represent a univariate correlation problem")
            return(TRUE)
        }
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Fisher-Yates (Normal Quantile) Correlation Test"
    object@parameter <- "rho"
    object@nullvalue <- 0

    return(object)
}


## quadrant test (Hajek, Sidak & Sen, pp. 124--125)
quadrant_test <- function(object, ...) UseMethod("quadrant_test")

quadrant_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("quadrant_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

quadrant_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate", "none"),
    mid.score = c("0", "0.5", "1"), ...) {
    ## <FIXME> in principle is "exact" also possible, unless mid.score == "0.5",
    ## since the data is effectively reduced to a 2x2 table.  But...
    ## </FIXME>

    args <- setup_args(
        teststat = "scalar",
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate", "none")
        ),
        xtrafo = function(data)
            trafo(data, numeric_trafo = function(x)
                median_trafo(x, mid.score = mid.score)),
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                median_trafo(y, mid.score = mid.score)),
        check = function(object) {
            if (!is_corr(object))
                stop(sQuote("object"),
                     " does not represent a univariate correlation problem")
            return(TRUE)
        }
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Quadrant Test"
    object@parameter <- "rho"
    object@nullvalue <- 0

    return(object)
}


## Koziol-Nemec test
koziol_test <- function(object, ...) UseMethod("koziol_test")

koziol_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("koziol_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

koziol_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate", "none"),
    ties.method = c("mid-ranks", "average-scores"), ...) {

    args <- setup_args(
        teststat = "scalar",
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate", "none")
        ),
        xtrafo = function(data)
            trafo(data, numeric_trafo = function(x)
                koziol_trafo(x, ties.method = ties.method)),
        ytrafo = function(data)
            trafo(data, numeric_trafo = function(y)
                koziol_trafo(y, ties.method = ties.method)),
        check = function(object) {
            if (!is_corr(object))
                stop(sQuote("object"),
                     " does not represent a univariate correlation problem")
            return(TRUE)
        }
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Koziol-Nemec Test"
    object@parameter <- "rho"
    object@nullvalue <- 0

    return(object)
}
