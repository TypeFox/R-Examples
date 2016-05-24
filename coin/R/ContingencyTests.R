### Pearson's chi-squared test
chisq_test <- function(object, ...) UseMethod("chisq_test")

chisq_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("chisq_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

chisq_test.table <- function(object, ...) {

    object <- table2IndependenceProblem(object)
    object <- do.call("chisq_test", c(list(object = object), list(...)))
    return(object)
}

chisq_test.IndependenceProblem <- function(object, ...) {

    check <- function(object) {
        if (!is_contingency(object))
            stop(sQuote("object"),
                 " does not represent a contingency problem")
        if (nlevels(object@block) != 1)
            stop(sQuote("object"), " contains blocks: use ",
                 sQuote("cmh_test"), " instead")
        return(TRUE)
    }
    n <- sum(object@weights)

    args <- setup_args()
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quadratic"
    ## distribution must be checked
    args$distribution <- check_distribution_arg(args$distribution)
    ## alternative is needed later
    args$alternative <- match.arg(args$alternative,
                                  c("two.sided", "less", "greater"))

    ## transform data if requested and setup a test problem
    object <- new("IndependenceTestProblem", object, args$xtrafo, args$ytrafo)

    if (!check(object))
        stop(sQuote("check"), " failed")

    object <- new("IndependenceLinearStatistic", object, varonly = FALSE)

    ## use the classical chisq statistic based on Pearson
    ## residuals (O - E)^2 / E
    ## see Th. 3.1 and its proof in Strasser & Weber (1999).
    object <-
        if (args$teststat == "scalar") {
            object <-
                new("ScalarIndependenceTestStatistic", object, args$alternative)
            object@teststatistic <- object@teststatistic * sqrt(n / (n - 1))
            object@standardizedlinearstatistic <-
                object@standardizedlinearstatistic * sqrt(n / (n - 1))
            object@covariance <-
                new("CovarianceMatrix", covariance(object) * (n - 1) / n)
            new("ScalarIndependenceTest", statistic = object,
                distribution = args$distribution(object))
        } else {
            if (args$alternative != "two.sided")
                warning(sQuote("alternative"),
                        " is ignored for quadratic test statistics")
            object <- new("QuadTypeIndependenceTestStatistic", object)
            object@teststatistic <- object@teststatistic * n / (n - 1)
            object@standardizedlinearstatistic <-
                object@standardizedlinearstatistic * sqrt(n / (n - 1))
            object@covariance <-
                new("CovarianceMatrix", covariance(object) * (n - 1) / n)
            object@covarianceplus <- MPinv(covariance(object))$MPinv
            new("QuadTypeIndependenceTest", statistic = object,
                distribution = args$distribution(object))
        }

    if (is_doubly_ordered(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else if (is_singly_ordered(object@statistic))
        object@method <- "Generalized Pearson Chi-Squared Test"
    else
        object@method <- "Pearson Chi-Squared Test"

    object@call <- match.call()

    return(object)
}


### generalized Cochran-Mantel-Haenzel test
cmh_test <- function(object, ...) UseMethod("cmh_test")

cmh_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("cmh_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

cmh_test.table <- function(object, ...) {

    object <- table2IndependenceProblem(object)
    object <- do.call("cmh_test", c(list(object = object), list(...)))
    return(object)
}

cmh_test.IndependenceProblem <- function(object, ...) {

    args <- setup_args(
        check = function(object) {
            if (!is_contingency(object))
                stop(sQuote("object"),
                     " does not represent a contingency problem")
            return(TRUE)
        }
    )
    ## convert factors to ordered and attach scores if requested
    if (!is.null(args$scores)) {
        object <- setscores(object, args$scores)
        args$scores <- NULL
    }
    ## set test statistic to scalar for linear-by-linear tests
    args$teststat <-
        if ((is.ordered(object@x[[1]]) && is.ordered(object@y[[1]])) ||
                ((is.ordered(object@x[[1]]) && nlevels(object@y[[1]]) == 2) ||
                 (is.ordered(object@y[[1]]) && nlevels(object@x[[1]]) == 2)))
            "scalar"
        else "quadratic"

    object <- do.call("independence_test", c(list(object = object), args))

    if (is_doubly_ordered(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else
        object@method <- "Generalized Cochran-Mantel-Haenszel Test"

    return(object)
}


### linear-by-linear association test
lbl_test <- function(object, ...) UseMethod("lbl_test")

lbl_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("lbl_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

lbl_test.table <- function(object, ...) {

    object <- table2IndependenceProblem(object)
    object <- do.call("lbl_test", c(list(object = object), list(...)))
    return(object)
}

lbl_test.IndependenceProblem <- function(object,
    distribution = c("asymptotic", "approximate", "none"), ...) {

    ## convert factors to ordered
    object@x[] <- lapply(object@x, function(x)
        if (is.factor(x) && nlevels(x) > 2) ordered(x) else x)
    object@y[] <- lapply(object@y, function(y)
        if (is.factor(y) && nlevels(y) > 2) ordered(y) else y)

    args <- setup_args(
        teststat = "scalar",
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate", "none")
        ),
        check = function(object) {
            if (!is_doubly_ordered(object))
                stop(sQuote("object"),
                     " does not represent a problem with ordered data")
            return(TRUE)
        }
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "Linear-by-Linear Association Test"

    return(object)
}
