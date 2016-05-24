#' A buggy version of plyr's count function
#' 
#' An old version of the \code{plyr} package's \code{count} function that
#' fails when you pass it a \code{factor} input.
#' @param df A data frame or an atomic input.
#' @param vars Variables in \code{df} to count unique values of
#' @param wt_var Optional variable to weight by.  See \code{\link[plyr]{count}}.
#' @return A data frame with label and freq columns.
#' @note In case the ``buggy'' part of the name didn't give it away, 
#' this is not suitable for real world usage!
#' @seealso
#' \code{\link[plyr]{count}} and \code{\link[base]{table}}
#' @examples
#' \dontrun{
#' buggy_count(factor()) #oops!
#' }
#' @export
buggy_count <- function (df, vars = NULL, wt_var = NULL) 
{
    if (is.vector(df)) {
        df <- data.frame(x = df)
    }
    if (!is.null(vars)) {
        vars <- as.quoted(vars)
        df2 <- quickdf(eval.quoted(vars, df))
    }
    else {
        df2 <- df
    }
    id <- ninteraction(df2, drop = TRUE)
    u_id <- !duplicated(id)
    labels <- df2[u_id, , drop = FALSE]
    labels <- labels[order(id[u_id]), , drop = FALSE]
    if (is.null(wt_var) && "freq" %in% names(df)) {
        message("Using freq as weighting variable")
        wt_var <- "freq"
    }
    if (!is.null(wt_var)) {
        wt_var <- plyr::as.quoted(wt_var)
        if (length(wt_var) > 1) {
            stop("wt_var must be a single variable", call. = FALSE)
        }
        wt <- plyr::eval.quoted(wt_var, df)[[1]]
        freq <- plyr::vaggregate(wt, id, sum, .default = 0)
    }
    else {
        freq <- tabulate(id, attr(id, "n"))
    }
    plyr::unrowname(data.frame(labels, freq))
}

# These next functions are taken directly from plyr.  They exist only to prevent having to call plyr's internals, which may change..

ninteraction <- function (.variables, drop = FALSE) 
{
    lengths <- vapply(.variables, length, integer(1))
    .variables <- .variables[lengths != 0]
    if (length(.variables) == 0) {
        n <- nrow(.variables) %||% 0L
        return(structure(seq_len(n), n = n))
    }
    if (length(.variables) == 1) {
        return(id_var(.variables[[1]], drop = drop))
    }
    ids <- rev(lapply(.variables, id_var, drop = drop))
    p <- length(ids)
    ndistinct <- vapply(ids, attr, "n", FUN.VALUE = numeric(1), 
        USE.NAMES = FALSE)
    n <- prod(ndistinct)
    if (n > 2^31) {
        char_id <- do.call("paste", c(ids, sep = "\r"))
        res <- match(char_id, unique(char_id))
    }
    else {
        combs <- c(1, cumprod(ndistinct[-p]))
        mat <- do.call("cbind", ids)
        res <- c((mat - 1L) %*% combs + 1L)
    }
    attr(res, "n") <- n
    if (drop) {
        id_var(res, drop = TRUE)
    }
    else {
        structure(as.integer(res), n = attr(res, "n"))
    }
}

as.quoted <- function (x, env = parent.frame()) 
{ 
    UseMethod("as.quoted")
}

as.quoted.call <- function (x, env = parent.frame()) 
{
    structure(as.list(x)[-1], env = env, class = "quoted")
}

as.quoted.character <- function (x, env = parent.frame()) 
{
    structure(lapply(x, function(x) parse(text = x)[[1]]), env = env, 
        class = "quoted")
}

as.quoted.factor <- function (x, env = parent.frame()) 
{
    as.quoted(as.character(x), env)
}

as.quoted.formula <- function (x, env = parent.frame()) 
{
    simplify <- function(x) {
        if (length(x) == 2 && x[[1]] == as.name("~")) {
            return(simplify(x[[2]]))
        }
        if (length(x) < 3) 
            return(list(x))
        op <- x[[1]]
        a <- x[[2]]
        b <- x[[3]]
        if (op == as.name("+") || op == as.name("*") || op == 
            as.name("~")) {
            c(simplify(a), simplify(b))
        }
        else if (op == as.name("-")) {
            c(simplify(a), bquote(-.(x), list(x = simplify(b))))
        }
        else {
            list(x)
        }
    }
    structure(simplify(x), env = env, class = "quoted")
}

as.quoted.namef <- function (x, env = parent.frame()) 
{
    structure(list(x), env = env, class = "quoted")
}

as.quoted.NULL <- function (x, env = parent.frame()) 
{
    structure(list(), env = env, class = "quoted")
}

as.quoted.numericf <- function (x, env = parent.frame()) 
{
    structure(x, env = env, class = c("quoted", "numeric"))
}

as.quoted.quotedf <- function (x, env = parent.frame()) 
{
    x
}

quickdf <- function (list) 
{
    rows <- unique(unlist(lapply(list, NROW)))
    stopifnot(length(rows) == 1)
    names(list) <- make_names(list, "X")
    class(list) <- "data.frame"
    attr(list, "row.names") <- c(NA_integer_, -rows)
    list
}

eval.quoted <- function (exprs, envir = NULL, enclos = NULL, try = FALSE) 
{
    if (is.numeric(exprs)) 
        return(envir[exprs])
    if (!is.null(envir) && !is.list(envir) && !is.environment(envir)) {
        stop("envir must be either NULL, a list, or an environment.")
    }
    qenv <- if (is.quoted(exprs)) 
        attr(exprs, "env")
    else parent.frame()
    if (is.null(envir)) 
        envir <- qenv
    if (is.data.frame(envir) && is.null(enclos)) 
        enclos <- qenv
    if (try) {
        results <- lapply(exprs, failwith(NULL, eval, quiet = TRUE), 
            envir = envir, enclos = enclos)
    }
    else {
        results <- lapply(exprs, eval, envir = envir, enclos = enclos)
    }
    names(results) <- names(exprs)
    results
}

is.quoted <- function (x) 
{
    inherits(x, "quoted")
}

failwith <- function (default = NULL, f, quiet = FALSE) 
{
    f <- match.fun(f)
    function(...) plyr::try_default(f(...), default, quiet = quiet)
}

`%||%` <- function (a, b) 
{
    if (is.null(a)) b else a
}

id_var <- function (x, drop = FALSE) 
{
    if (length(x) == 0) 
        return(structure(integer(), n = 0L))
    if (!is.null(attr(x, "n")) && !drop) 
        return(x)
    if (is.factor(x) && !drop) {
        id <- as.integer(addNA(x, ifany = TRUE))
        n <- length(levels(x))
    }
    else {
        levels <- sort(unique(x), na.last = TRUE)
        id <- match(x, levels)
        n <- max(id)
    }
    structure(id, n = n)
}

make_names <- function (x, prefix = "X") 
{
    nm <- names(x)
    if (is.null(nm)) {
        nm <- rep.int("", length(x))
    }
    n <- sum(nm == "", na.rm = TRUE)
    nm[nm == ""] <- paste(prefix, seq_len(n), sep = "")
    nm
}
