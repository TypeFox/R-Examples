### sign test
sign_test <- function(object, ...) UseMethod("sign_test")

sign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    object <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(object$bl)) {
        if (is.Surv(object$y[[1]]))
            stop(sQuote("y"), " is not a numeric variable")
        if (is.Surv(object$x[[1]]))
            stop(sQuote("x"), " is not a numeric variable")
        object <- list(y = data.frame(y = c(object$y[[1]], object$x[[1]])),
                       x = data.frame(x = gl(2, length(object$x[[1]]))),
                       bl = factor(rep.int(1:length(object$x[[1]]), 2)))
    }
    object <- new("SymmetryProblem", x = object$x, y = object$y,
                  block = object$bl)
    object <- do.call("sign_test", c(list(object = object), list(...)))
    return(object)
}

sign_test.SymmetryProblem <- function(object, ...) {

    y <- object@y[[1]]
    x <- object@x[[1]]

    if (!is_numeric_y(object))
        stop(sQuote("y"), " is not a numeric variable")
    if (is_2sample(object))
        diffs <- tapply(1:length(y), object@block, function(b)
            y[b][x[b] == levels(x)[1]] - y[b][x[b] == levels(x)[2]])
    else
        stop(sQuote("object"),
             " does not represent a paired two-sample problem",
             " (maybe the grouping variable is not a factor?)")

    abs_diffs <- abs(diffs)
    if (all(abs_diffs < .Machine$double.eps))
        stop("all pairwise differences equal zero")

    diffs <- diffs[abs_diffs > 0]
    n <- length(diffs)

    object <- new("SymmetryProblem",
                  x = data.frame(x = factor(rep.int(0:1, n),
                                            labels = c("pos", "neg"))),
                  y = data.frame(y = as.vector(rbind(as.numeric(diffs > 0),
                                                     as.numeric(diffs < 0)))),
                  block = gl(n, 2))

    args <- setup_args(teststat = "scalar", paired = TRUE)

    object <- do.call("symmetry_test", c(list(object = object), args))

    object@method <- "Sign Test"
    object@nullvalue <- 0

    return(object)
}


### Wilcoxon signed-rank test
wilcoxsign_test <- function(object, ...) UseMethod("wilcoxsign_test")

wilcoxsign_test.formula <- function(formula, data = list(), subset = NULL, ...)
{
    object <- formula2data(formula, data, subset, frame = parent.frame(), ...)
    if (is.null(object$bl)) {
        if (is.Surv(object$y[[1]]))
            stop(sQuote("y"), " is not a numeric variable")
        if (is.Surv(object$x[[1]]))
            stop(sQuote("x"), " is not a numeric variable")
        object <- list(y = data.frame(y = c(object$y[[1]], object$x[[1]])),
                       x = data.frame(x = gl(2, length(object$x[[1]]))),
                       bl = factor(rep.int(1:length(object$x[[1]]), 2)))
    }
    object <- new("SymmetryProblem", x = object$x, y = object$y,
                  block = object$bl)
    object <- do.call("wilcoxsign_test", c(list(object = object), list(...)))
    return(object)
}

wilcoxsign_test.SymmetryProblem <- function(object,
    zero.method = c("Pratt", "Wilcoxon"), ...) {

    zero.method <- match.arg(zero.method)

    y <- object@y[[1]]
    x <- object@x[[1]]

    if (!is_numeric_y(object))
        stop(sQuote("y"), " is not a numeric variable")
    if (is_2sample(object))
        diffs <- tapply(1:length(y), object@block, function(b)
            y[b][x[b] == levels(x)[1]] - y[b][x[b] == levels(x)[2]])
    else
        stop(sQuote("object"),
             " does not represent a paired two-sample problem",
             " (maybe the grouping variable is not a factor?)")

    abs_diffs <- abs(diffs)
    if (all(abs_diffs < .Machine$double.eps))
        stop("all pairwise differences equal zero")

    pos_abs_diffs <- abs_diffs > 0
    if (zero.method == "Pratt") {
        rank_abs_diffs <- rank(abs_diffs)
        pos <- (rank_abs_diffs * (diffs > 0))[pos_abs_diffs]
        neg <- (rank_abs_diffs * (diffs < 0))[pos_abs_diffs]
    } else {
        diffs <- diffs[pos_abs_diffs]
        abs_diffs <- abs_diffs[pos_abs_diffs]
        rank_abs_diffs <- rank(abs_diffs)
        pos <- rank_abs_diffs * (diffs > 0)
        neg <- rank_abs_diffs * (diffs < 0)
    }
    n <- length(pos)

    object <- new("SymmetryProblem",
                  x = data.frame(x = factor(rep.int(0:1, n),
                                            labels = c("pos", "neg"))),
                  y = data.frame(y = as.vector(rbind(pos, neg))),
                  block = gl(n, 2))

    args <- setup_args(teststat = "scalar", paired = TRUE)

    object <- do.call("symmetry_test", c(list(object = object), args))

    if (zero.method == "Pratt")
        object@method <- "Wilcoxon-Pratt Signed-Rank Test"
    else
        object@method <- "Wilcoxon Signed-Rank Test"
    object@nullvalue <- 0

    return(object)
}


### Friedman Test
friedman_test <- function(object, ...) UseMethod("friedman_test")

friedman_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("friedman_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

friedman_test.SymmetryProblem <- function(object, ...) {

    args <- setup_args(
        ytrafo = function(data)
            trafo(data, numeric_trafo = rank_trafo,
                  block = object@block),
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

    object <- do.call("symmetry_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Page Test"
    else
        object@method <- "Friedman Test"

    return(object)
}


### Quade Test
quade_test <- function(object, ...) UseMethod("quade_test")

quade_test.formula <- function(formula, data = list(), subset = NULL, ...) {

    ft("quade_test", "SymmetryProblem", formula, data, subset,
       frame = parent.frame(), ...)
}

quade_test.SymmetryProblem <- function(object, ...) {

    args <- setup_args(
        check = function(object) {
            if (!is_Ksample(object))
                stop(sQuote("object"),
                     " does not represent a K-sample problem",
                     " (maybe the grouping variable is not a factor?)")
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

    if (!is_numeric_y(object))
        stop(sQuote(colnames(object@y)), " is not a numeric variable")
    y <- split(object@y[[1]], object@block)
    R <- lapply(y, function(y) rank(y) - (length(y) + 1) / 2)
    Q <- rank(vapply(y, function(y) max(y) - min(y), NA_real_, USE.NAMES = FALSE))
    object@y[[1]] <- unsplit(lapply(seq_along(Q), function(i) Q[i] * R[[i]]),
                             object@block)

    object <- do.call("symmetry_test", c(list(object = object), args))

    if (is_ordered_x(object@statistic))
        object@method <- "Linear-by-Linear Association Test"
    else
        object@method <- "Quade Test"

    return(object)
}
