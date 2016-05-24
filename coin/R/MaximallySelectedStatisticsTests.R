### generalized maximally selected statistics
maxstat_test <- function(object, ...) UseMethod("maxstat_test")

maxstat_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("maxstat_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

maxstat_test.table <- function(object, ...) {

    object <- table2IndependenceProblem(object)
    object <- do.call("maxstat_test", c(list(object = object), list(...)))
    return(object)
}

maxstat_test.IndependenceProblem <- function(object,
    teststat = c("maximum", "quadratic"),
    distribution = c("asymptotic", "approximate", "none"),
    minprob = 0.1, maxprob = 1 - minprob, ...) {

    args <- setup_args(
        teststat = match.arg(teststat),
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate", "none")
        ),
        xtrafo = function(data)
            trafo(data,
                  numeric_trafo = function(x)
                      maxstat_trafo(x, minprob = minprob, maxprob = maxprob),
                  factor_trafo = function(x)
                      fmaxstat_trafo(x, minprob = minprob, maxprob = maxprob),
                  ordered_trafo = function(x)
                      ofmaxstat_trafo(x, minprob = minprob, maxprob = maxprob))
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Generalized Maximally Selected Statistics"

    ## estimate cutpoint
    wm <- which.max(apply(abs(statistic(object, "standardized")), 1, max))
    whichvar <- attr(object@statistic@xtrans, "assign")[wm]
    maxcontr <- object@statistic@xtrans[, wm]
    if (is.factor(object@statistic@x[[whichvar]])) {
        cp <- levels(object@statistic@x[[whichvar]][maxcontr > 0][, drop = TRUE])
        cp0 <- levels(object@statistic@x[[whichvar]][maxcontr == 0][, drop = TRUE])
        lab <- paste0("{", paste0(cp, collapse = ", "), "} vs. {",
                      paste0(cp0, collapse = ", "), "}")
    } else {
        cp <- max(object@statistic@x[[whichvar]][maxcontr > 0])
        lab <- paste0("<= ", format(cp, digits = getOption("digits")))
    }
    if (ncol(object@statistic@x) > 1)
        estimate <- list(covariable = colnames(object@statistic@x)[whichvar],
                         cutpoint = cp, label = lab)
    else
        estimate <- list(cutpoint = cp, label = lab)
    class(estimate) <- c("cutpoint", "list")
    object@estimates <- list(estimate = estimate)

    return(object)
}
