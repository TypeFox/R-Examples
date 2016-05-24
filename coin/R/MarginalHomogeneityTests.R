### marginal homogeneity test
mh_test <- function(object, ...) UseMethod("mh_test")

mh_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("mh_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

mh_test.table <- function(object, ...) {

    object <- table2df_sym(object)
    object <- new("SymmetryProblem", x = object["conditions"],
                  y = object["response"])
    object <- do.call("mh_test", c(list(object = object), list(...)))
    return(object)
}

mh_test.SymmetryProblem <- function(object, ...) {

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

    object <- do.call("symmetry_test", c(list(object = object), args))

    if (is_ordered(object@statistic))
        object@method <- "Marginal Homogeneity Test for Ordered Data"
    else
        object@method <- "Marginal Homogeneity Test"

    return(object)
}
