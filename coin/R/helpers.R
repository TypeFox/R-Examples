asymptotic <- function(maxpts = 25000, abseps = 0.001, releps = 0) {
    function(object)
        AsymptNullDistribution(object, maxpts = maxpts, abseps = abseps,
                               releps = releps)
}

approximate <- function(B = 10000, parallel = c("no", "multicore", "snow"),
                        ncpus = 1, cl = NULL) {
    parallel <- match.arg(parallel)
    function(object)
        ApproxNullDistribution(object, B = B, parallel = parallel,
                               ncpus = ncpus, cl = cl)
}

exact <- function(algorithm = c("auto", "shift", "split-up"), fact = NULL) {
    algorithm <- match.arg(algorithm)
    function(object)
        ExactNullDistribution(object, algorithm = algorithm, fact = fact)
}

LinearStatistic <- function(x, y, weights) {
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_LinearStatistic", x, y, weights, PACKAGE = "coin")
}

ExpectCovarInfluence <- function(y, weights) {
    storage.mode(y) <- "double"
    storage.mode(weights) <- "double"
    .Call("R_ExpectCovarInfluence", y, weights, PACKAGE = "coin")
}

expectvaronly <- function(x, y, weights) {
    indx <- rep.int(seq_len(nrow(x)), weights)
    x <- x[indx, , drop = FALSE]
    y <- y[indx, , drop = FALSE]
    n <- nrow(x)
    Ey <- colMeans(y)
    Vy <- rowMeans((t(y) - Ey)^2)

    rSx <- colSums(x)
    rSx2 <- colSums(x^2)
    ## in case rSx _and_ Ey are _both_ vectors
    E <- .Call("R_kronecker", Ey, rSx, PACKAGE = "coin")
    V <- n / (n - 1) * .Call("R_kronecker", Vy, rSx2, PACKAGE = "coin")
    V <- V - 1 / (n - 1) * .Call("R_kronecker", Vy, rSx^2, PACKAGE = "coin")
    list(E = drop(E), V = matrix(V, nrow = 1L))
}

ExpectCovarLinearStatistic <- function(x, y, weights, varonly = FALSE) {
    if (varonly) {
        ev <- expectvaronly(x, y, weights)
### <FIXME> This conflicts with the "new" initialize methods that were brought
###         over from 'party' in r927
###         new("ExpectCovar", expectation = ev$E, covariance = ev$V)
### </FIXME>
        ec <- new("ExpectCovar")
        ec@expectation <- ev$E
        ec@covariance <- ev$V
        ec
    } else {
        storage.mode(x) <- "double"
        storage.mode(y) <- "double"
        storage.mode(weights) <- "double"
        expcovinf <- ExpectCovarInfluence(y, weights)
        .Call("R_ExpectCovarLinearStatistic", x, y, weights, expcovinf,
               PACKAGE = "coin")
    }
}

pmvn <- function(lower, upper, mean, corr, ..., conf.int = TRUE) {
    p <- if (length(corr) > 1L)
             pmvnorm(lower = lower, upper = upper, mean = mean,
                     corr = corr, ...)
         else
             pmvnorm(lower = lower, upper = upper, mean = mean,
                     sigma = 1L, ...)
    if (conf.int) {
        error <- attr(p, "error")
        attributes(p) <- NULL
        ci <- c(max(0, p - error), min(p + error, 1))
        attr(ci, "conf.level") <- 0.99
        attr(p, "conf.int") <- ci
    } else
        attributes(p) <- NULL
    p
}

qmvn <- function(p, mean, corr, ...) {
    q <- if (length(corr) > 1L)
             qmvnorm(p = p, mean = mean, corr = corr,
                     tail = "both.tails", ...)$quantile
         else
             qmvnorm(p = p, mean = mean, sigma = 1,
                     tail = "both.tails", ...)$quantile
    attributes(q) <- NULL
    q
}

### copied from package MASS
MPinv <- function (X, tol = eps())
{
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    RET <- if (all(Positive))
               Xsvd$v %*% (1 / Xsvd$d * t(Xsvd$u))
           else if (!any(Positive))
               array(0, dim(X)[2L:1L])
           else
               Xsvd$v[, Positive, drop = FALSE] %*%
                 ((1 / Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
    list(MPinv = RET, rank = sum(Positive))
}

copyslots <- function(source, target) {
    slots <- names(getSlots(class(source)))
    slots <- slots[(slots %in% names(getSlots(class(target))))]
    if (length(slots) == 0)
        stop("no common slots to copy to")
    for (s in slots)
        eval(parse(text = paste("target@", s, " <- source@", s)))
    return(target)
}

ft <- function(test, class, formula, data = list(), subset = NULL,
    weights = NULL, ...) {

    object <- formula2data(formula, data, subset, weights = weights, ...)
    object <- new(class, x = object$x, y = object$y, block = object$bl,
                  weights = object$w)
    args <- list(...)
    args$frame <- NULL

    ## warn users of weighted rank tests
    if (test %in% ranktests && !is.null(object@weights) &&
        !is_unity(object@weights))
        warning("Rank transformation doesn't take weights into account")

    do.call(test, c(list(object = object), args))
}

ranktests <-
    c("wilcox_test", "kruskal_test", "normal_test", "median_test",
      "savage_test", "taha_test", "klotz_test", "mood_test", "ansari_test",
      "fligner_test", "conover_test", "logrank_test", "quade_test",
      "friedman_test", "wilcoxsign_test", "spearman_test", "fisyat_test",
      "quadrant_test", "koziol_test")

formula2data <- function(formula, data, subset, weights = NULL, ...) {
    no_weights <- is.null(weights)

    dat <- ModelEnvFormula(
        formula = formula,
        data = data,
        subset = subset,
        other = if (no_weights) list() else list(weights = weights),
        na.action = na.omit,
        designMatrix = FALSE, responseMatrix = FALSE,
        ...
    )

    ## rhs of formula
    if (has(dat, "input"))
        x <- dat@get("input")
    else
        stop("missing right hand side of formula")

    ## ~ x + y is allowed
    if (has(dat, "response"))
        y <- dat@get("response")
    else {
        if (ncol(x) == 2L) {
            y <- x[2L]
            x <- x[1L]
        } else
            stop("missing left hand side of formula")
    }

    ## y ~ x | block or ~ y + x | block
    if (has(dat, "blocks")) {
        block <- dat@get("blocks")
        attr(block[[1L]], "blockname") <- colnames(block)
    } else
        block <- NULL

    list(x = x, y = y, block = block, bl = block[[1L]],
         w = if (no_weights) NULL else dat@get("weights")[[1L]])
}

setscores <- function(x, scores) {

    if (is.null(scores)) return(x)

    varnames <- names(scores)
    if (!is.list(scores) || is.null(varnames))
       stop(sQuote("scores"), " is not a named list")

    missing <- varnames[!varnames %in% c(colnames(x@x), colnames(x@y))]
    if (length(missing) > 0L)
        stop("Variable(s) ", paste(missing, sep = ", "),
             " not found in ", sQuote("x"))

    for (var in varnames) {
        if (!is.null(x@x[[var]])) {
            if (!is.factor(x@x[[var]]))
                stop(sQuote(var), " is not a factor")
            if (nlevels(x@x[[var]]) != length(scores[[var]]))
                stop("scores for variable ", sQuote(var), " don't match")
            x@x[[var]] <- ordered(x@x[[var]], levels = levels(x@x[[var]]))
            if (nlevels(x@x[[var]]) == 2)
                scores[[var]] <- 0:1      # must be 0:1 for exact p-values
            attr(x@x[[var]], "scores") <- scores[[var]]
        }
        if (!is.null(x@y[[var]])) {
            if (!is.factor(x@y[[var]]))
                stop(sQuote(var), " is not a factor")
            if (nlevels(x@y[[var]]) != length(scores[[var]]))
                stop("scores for variable ", sQuote(var), " don't match")
            x@y[[var]] <- ordered(x@y[[var]], levels = levels(x@y[[var]]))
            if (nlevels(x@y[[var]]) == 2)
                scores[[var]] <- 0:1      # must be 0:1 for exact p-values
            attr(x@y[[var]], "scores") <- scores[[var]]
        }
    }
    return(x)
}

### user-supplied trafo functions may return a vector or matrices
### with NROW being equal for the x and y variables
check_trafo <- function(tx, ty) {

    if (!(is.numeric(tx) || is.logical(tx)))
        stop(sQuote("xtrafo"), " does not return a numeric or logical vector")
    if (!(is.numeric(ty) || is.logical(ty)))
        stop(sQuote("ytrafo"), " does not return a numeric or logical vector")
    if (NROW(tx) != NROW(ty))
        stop("Dimensions of returns of ", sQuote("xtrafo"), " and ",
             sQuote("ytrafo"), " don't match")
    if (!is.matrix(tx)) tx <- matrix(tx, ncol = 1L)
    if (!is.matrix(ty)) ty <- matrix(ty, ncol = 1L)
    storage.mode(tx) <- "double"
    storage.mode(ty) <- "double"
    list(xtrafo = tx, ytrafo = ty)
}

table2df <- function(x) {
    if (!is.table(x))
        stop(sQuote("x"), " is not of class ", sQuote("table"))
    x <- as.data.frame(x)
    freq <- x[["Freq"]]
    x <- x[rep.int(seq_len(nrow(x)), freq), , drop = FALSE]
    rownames(x) <- seq_len(nrow(x))
    return(x[, colnames(x) != "Freq"])
}

table2df_sym <- function(x) {
    x <- table2df(x)
    lx <- levels(x[[1L]])
    if (!all(vapply(x, function(x) all(levels(x) == lx), NA)))
        stop("table ", sQuote("x"), " does not represent a symmetry problem")
    data.frame(conditions = factor(rep.int(colnames(x),
                                           rep.int(nrow(x), ncol(x)))),
               response = factor(unlist(x, recursive = FALSE,
                                        use.names = FALSE),
                                 labels = lx))
}

table2IndependenceProblem <- function(object) {

    df <- as.data.frame(object)
    if (ncol(df) == 3L)
        new("IndependenceProblem",
            x = df[1L], y = df[2L], block = NULL, weights = df[["Freq"]])
    else if (ncol(df) == 4L) {
        attr(df[[3L]], "blockname") <- colnames(df)[3L]
        new("IndependenceProblem",
            x = df[1L], y = df[2L], block = df[[3L]], weights = df[["Freq"]])
    } else
        stop(sQuote("object"), " is not a two- or three-way contingency table")
}

is_ytrafo <- function()
    any(vapply(sys.calls(), function(i)
            identical(as.character(i)[1], "ytrafo"), NA))

is_factor_y <- function(object)
    ncol(object@y) == 1L && is.factor(object@y[[1L]])

is_factor_x <- function(object)
    ncol(object@x) == 1L && is.factor(object@x[[1L]])

is_ordered_y <- function(object)
    ncol(object@y) == 1L && is.ordered(object@y[[1L]])

is_ordered_x <- function(object)
    ncol(object@x) == 1L && is.ordered(object@x[[1L]])

is_unordered_y <- function(object)
    ncol(object@y) == 1L && is.factor(object@y[[1L]]) && !is.ordered(object@y[[1L]])

is_unordered_x <- function(object)
    ncol(object@x) == 1L && is.factor(object@x[[1L]]) && !is.ordered(object@x[[1L]])

is_numeric_y <- function(object)
    ncol(object@y) == 1L && is.numeric(object@y[[1L]]) && !is.Surv(object@y[[1L]])

is_numeric_x <- function(object)
    ncol(object@x) == 1L && is.numeric(object@x[[1L]]) && !is.Surv(object@x[[1L]])

is_censored_y <- function(object)
    ncol(object@y) == 1L && is.Surv(object@y[[1L]])

is_censored_x <- function(object)
    ncol(object@x) == 1L && is.Surv(object@x[[1L]])

is_2sample <- function(object)
    ncol(object@x) == 1L && nlevels(object@x[[1L]]) == 2L

is_Ksample <- is_factor_x

is_corr <- function(object)
    (is_numeric_y(object) || is_censored_y(object)) &&
        (is_numeric_x(object) || is_censored_x(object))

is_contingency <- function(object)
    is_factor_y(object) && is_factor_x(object)

is_contingency_2x2 <- function(object)
    is_contingency(object) && nlevels(object@y[[1L]]) == 2L &&
        nlevels(object@x[[1L]]) == 2L

is_singly_ordered <- function(object) {
    ## NOTE: unordered factors with exactly 2 levels are regarded as ordered
    (is_ordered_y(object) &&
     (is_numeric_x(object) || is_censored_x(object) ||
      (is_unordered_x(object) && nlevels(object@x[[1L]]) > 2L))) ||
    (is_ordered_x(object) &&
     (is_numeric_y(object) || is_censored_y(object) ||
      (is_unordered_y(object) && nlevels(object@y[[1L]]) > 2L)))
}

is_doubly_ordered <- function(object) {
    ## NOTE: unordered factors with exactly 2 levels are regarded as ordered
    (is_ordered_y(object) &&
     (is_ordered_x(object) ||
      (is_unordered_x(object) && nlevels(object@x[[1L]]) == 2L))) ||
    (is_ordered_x(object) &&
     (is_ordered_y(object) ||
      (is_unordered_y(object) && nlevels(object@y[[1L]]) == 2L)))
}

is_ordered <- function(object)
    is_singly_ordered(object) || is_doubly_ordered(object)

is_completeblock <- function(object)
    all(table(object@x[[1L]], object@block) == 1L)

is_scalar <- function(object)
    ncol(object@xtrans) == 1L && ncol(object@ytrans) == 1L

is_integer <- function(x, fact = NULL) {
    if (is.null(fact))
        fact <- c(1, 2, 10, 100, 1000, 10000, 100000)
    f <- vapply(fact, function(f) max(abs(round(x * f) - (x * f))) < eps(), NA)
    if (RET <- any(f))
        attr(RET, "fact") <- min(fact[f])
    RET
}

is_monotone <- function (x)
    all(x == cummax(x)) || all(x == cummin(x))

isequal <- function(a, b) {
    attributes(a) <- NULL
    attributes(b) <- NULL
    if (!isTRUE(all.equal(a, b))) {
        print(a, digits = 10)
        print(b, digits = 10)
        return(FALSE)
    } else {
        return(TRUE)
    }
}

has_distribution <- function(args)
    !(is.character(args$distribution) && args$distribution == "none")

check_distribution_arg <- function(distribution,
    values = c("asymptotic", "approximate", "exact", "none")) {
    if (is.character(distribution)) {
        distribution <- match.arg(distribution, values)
        if (distribution == "none")
            function(object) new("NullDistribution")
        else
            eval(parse(text = paste0(distribution, "()")))
    } else
        distribution
}

setup_args <- function(...) {
    cl <- match.call(independence_test.IndependenceProblem,
                     call = sys.call(sys.parent()), expand.dots = FALSE)
    ## find default arguments and values
    args <- formals(independence_test.IndependenceProblem)
    args$object <- args$... <- NULL
    nm <- names(args)
    ## replace default values with user-specified values
    for (i in nm[nm %in% names(cl)])
        args[[i]] <- cl[[i]]
    ## override default and user-specified values
    for (i in nm[nm %in% names(list(...))])
        args[[i]] <- list(...)[[i]]
    lapply(args, eval.parent)
}

statnames <- function(object) {
    nc <- ncol(object@ytrans)
    nr <- ncol(object@xtrans)
    dn <- list(colnames(object@xtrans),
               colnames(object@ytrans))
    if (is.null(dn[[1L]])) {
        if (nr == 1L) {
            dn[[1L]] <- ""
        } else {
            dn[[1L]] <- paste0("X", seq_len(nr))
        }
    }
    if (is.null(dn[[2L]])) {
        if (nc == 1L) {
            dn[[2L]] <- ""
        } else {
            dn[[2L]] <- paste0("Y", seq_len(nc))
        }
    }
    list(dimnames = dn,
         names = paste(rep.int((dn[[1L]]), nc),
                       rep.int((dn[[2L]]), rep.int(nr, nc)),
                       sep = ifelse(dn[[1L]] == "" || dn[[2L]] == "", "", ":")))
}

varnames <- function(object) {
    yordered <- vapply(object@y, is.ordered, NA)
    ynames <- paste0(colnames(object@y), ifelse(yordered, " (ordered)", ""),
                     collapse = ", ")

    if (length(object@x) == 1) {
        if (is.ordered(object@x[[1]])) {
            xnames <- paste0(
                colnames(object@x), " (",
                paste0(levels(object@x[[1]]), collapse = " < "),
                ")"
            )
        } else {
            if (is.factor(object@x[[1]])) {
                xnames <- paste0(
                    colnames(object@x), " (",
                    paste0(levels(object@x[[1]]), collapse = ", "),
                    ")"
                )
            } else {
                xnames <- colnames(object@x)
            }
        }
    } else {
        xordered <- vapply(object@x, is.ordered, NA)
        xnames <- paste0(colnames(object@x), ifelse(xordered, "(ordered)", ""),
                         collapse = ", ")
    }

    if (nlevels(object@block) > 1) {
        bn <- attr(object@block, "blockname")
        if (is.null(bn))
            bn <- "block"
        xnames <- paste(xnames, paste("\n\t stratified by", bn))
    }

    if (nchar(xnames) > options("width")$width / 2)
        paste(ynames, "by\n\t", xnames, collapse = "")
    else
        paste(ynames, "by", xnames, collapse = "")
}

eps <- function() sqrt(.Machine$double.eps)

EQ <- function(x, y)
    abs(x - y) < eps()

GE <- function(x, y)
    x > y | abs(x - y) < eps()

LE <- function(x, y)
    x < y | abs(x - y) < eps()

### don't use! never!
get_weights <- function(object) object@statistic@weights
get_xtrans <- function(object) object@statistic@xtrans
get_ytrans <- function(object) object@statistic@ytrans

is_unity <- function(x)
    max(abs(x - 1.0)) < eps()

setRownames <- function (object, nm) {
    rownames(object) <- nm
    object
}

setColnames <- function (object, nm) {
    colnames(object) <- nm
    object
}

n_decimal_digits <- function(x)
    nchar(sub("^[[:digit:]]*[.]", "", format(min(x), scientific = FALSE)))

if(getRversion() < "2.15.0")
    paste0 <- function(...) paste(..., sep = "")

if(getRversion() < "3.1.0") {
    cospi <- function(x) cos(pi * x)
    anyNA <- function(x) any(is.na(x))
}

if(getRversion() < "3.2.0") {
    isNamespaceLoaded <- function(name) !is.null(.getNamespace(name))
    lengths <- function(x, use.names = TRUE)
        vapply(x, length, NA_integer_, USE.NAMES = use.names)
}
