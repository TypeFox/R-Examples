### a generic test procedure for classical (and not so classical) tests
symmetry_test <- function(object, ...) UseMethod("symmetry_test")

symmetry_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("symmetry_test", "SymmetryProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

symmetry_test.table <- function(object, ...) {

    object <- table2df_sym(object)
    object <- new("SymmetryProblem", x = object["conditions"],
                  y = object["response"])
    object <- do.call("symmetry_test", c(list(object = object), list(...)))
    return(object)
}

symmetry_test.SymmetryProblem <- function(object,
    teststat = c("maximum", "quadratic", "scalar"),
    distribution = c("asymptotic", "approximate", "exact", "none"),
    alternative = c("two.sided", "less", "greater"),
    xtrafo = trafo, ytrafo = trafo, scores = NULL, check = NULL, paired = FALSE,
    ...) {

    addargs <- list(...)
    if (length(addargs) > 0L)
        warning("additional arguments ",
                paste0(names(addargs), collapse = ", "),
                " will be ignored")

    teststat <- match.arg(teststat)
    alternative <- match.arg(alternative)
    distribution <- check_distribution_arg(distribution)

    ## convert factors to ordered and attach scores if requested
    if (!is.null(scores))
        object <- setscores(object, scores)

    ## transform data if requested and setup a test problem
    object <- new("IndependenceTestProblem", object, xtrafo = xtrafo,
                  ytrafo = ytrafo, ...)

    if (!is.null(check)) {
        if (is.function(check)) {
            if (!check(object))
                stop(sQuote("check"), " failed")
        } else {
            stop(sQuote("check"), " is not a function")
        }
    }

    ## check type of test statistic and alternative
    if (!is_scalar(object)) {
        if (teststat == "scalar") {
            warning("Length linear statistic > 1, using ",
                    sQuote("maximum"), "-type test statistic")
            teststat <- "maximum"
        }
    } else {
        if (teststat == "maximum") teststat <- "scalar"
    }
    if (alternative != "two.sided" && teststat == "quadratic")
        warning(sQuote("alternative"), " is ignored for ",
                teststat, " test statistics")

    ## compute linear statistic, conditional expectation and
    ## conditional covariance
    object <- new("IndependenceLinearStatistic", object, varonly = TRUE)
###         varonly = class(distribution) == "approximate" && teststat == "maximum")

    ## compute test statistic and corresponding null distribution
    ## return object inheriting from class `IndependenceTest'
    switch(teststat,
        "scalar" = {
            object <- new("ScalarIndependenceTestStatistic", object,
                          alternative = alternative, paired = paired)
            new("ScalarIndependenceTest", statistic = object,
                distribution = distribution(object),
                method = "General Symmetry Test", call = match.call())
        },
        "maximum" = {
            object <- new("MaxTypeIndependenceTestStatistic", object,
                          alternative = alternative)
            new("MaxTypeIndependenceTest", statistic = object,
                distribution = distribution(object),
                method = "General Symmetry Test", call = match.call())
        },
        "quadratic" = {
            object <- new("QuadTypeIndependenceTestStatistic", object,
                          paired = paired)
            new("QuadTypeIndependenceTest", statistic = object,
                distribution = distribution(object),
                method = "General Symmetry Test", call = match.call())
        }
    )
}
