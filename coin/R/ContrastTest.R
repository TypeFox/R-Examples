### contrast test
contrast_test <- function(object, ...) UseMethod("contrast_test")

contrast_test.formula <- function(formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    ft("contrast_test", "IndependenceProblem", formula, data, subset, weights,
       frame = parent.frame(), ...)
}

contrast_test.IndependenceProblem <- function(object,
    cmatrix, distribution = c("asymptotic", "approximate"), ...) {

    if (!(ncol(object@x) == 1 && is.factor(object@x[[1]])))
        stop(sQuote("object@x"), " is not univariate or a factor")

    if  (!is.matrix(cmatrix) || nrow(cmatrix) != nlevels(object@x[[1]]))
        stop(sQuote("cmatrix"), " is not a matrix with ",
             nlevels(object@x), " rows")

    if (is.null(colnames(cmatrix)))
        colnames(cmatrix) <- paste0("C", 1:ncol(cmatrix))

    args <- setup_args(
        teststat = "maximum",
        distribution = check_distribution_arg(
            distribution, values = c("asymptotic", "approximate")
        ),
        xtrafo = function(data)
            trafo(data) %*% cmatrix
    )

    object <- do.call("independence_test", c(list(object = object), args))

    object@method <- "General Contrast Test"

    return(object)
}
