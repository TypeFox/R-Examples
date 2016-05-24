dist <-
function(x, y = NULL, method = NULL, ...,
         diag = FALSE, upper = FALSE, pairwise = FALSE,
         by_rows = TRUE, convert_similarities = TRUE,
         auto_convert_data_frames = TRUE)
{

### PARAMETER HANDLING
    ## convenience hack to allow dist(x, "method")
    if ((is.function(y) || is.character(y)) && is.null(method)) {
        method <- y
        y <- NULL
    }

    ## transform data frame into matrix iff all columns are atomic and either numeric (integer or double) or logical or complex
    is.n_l_c <- function(x)
        all(sapply(x, is.numeric)) ||
        all(sapply(x, is.logical)) ||
        all(sapply(x, is.complex))
    if (is.data.frame(x) && auto_convert_data_frames && is.n_l_c(x))
        x <- as.matrix(x)
    if (is.data.frame(y) && !is.null(y) && auto_convert_data_frames && is.n_l_c(y))
        y <- as.matrix(y)

    ## vector handling
    if (is.vector(x) && is.atomic(x))
        x <- as.matrix(x)
    if (!is.null(y) && is.vector(y) && is.atomic(y))
        y <- as.matrix(y)

    ## method lookup
    reg_entry <- NULL
    if (is.null(method))
        method <- if (is.data.frame(x))
            "Gower"
        else if (is.logical(x))
            "Jaccard"
        else
            "Euclidean"
    if (!is.function(method))
        reg_entry <- if (inherits(method, "proxy_registry_entry"))
            method
        else
            pr_DB$get_entry(method)

    ## some checks
    if (!is.data.frame(x) && !is.matrix(x) && !is.list(x))
        stop("Can only handle data frames, vectors, matrices, and lists!")
    if ( is.data.frame(x) && !by_rows)
        stop("Cannot transpose mixed data frames")
    if (!is.null(y)) {
        if (is.data.frame(x) && !is.data.frame(y)
            || is.matrix(x) && !is.matrix(y)
            || is.list(x) && !is.list(y))
            stop("x and y must be of same type.")
        if (is.matrix(x) && is.matrix(y) || is.data.frame(x) && is.data.frame(y))
            if (by_rows && (ncol(x) != ncol(y)))
                stop("x and y must be conform in columns.")
            else if (!by_rows && (nrow(x) != nrow(y)))
                stop("x and y must be conform in rows.")
    }

### PREPROCESS
    params <- list(...)
    if (!is.null(reg_entry)) {
        if(!is.na(reg_entry$PREFUN)) {
            tmp <- do.call(reg_entry$PREFUN,
                           c(list(x, y, pairwise, params, reg_entry)))
            if (!is.null(tmp)) {
                x <- tmp$x
                y <- tmp$y
                pairwise <- tmp$pairwise
                params <- tmp$p
                reg_entry <- tmp$reg_entry
            }
        }
        method <- reg_entry$FUN
    }

    ## helper function for calling the C-level loops
    .proxy_external <- function(CFUN, x, y)
        do.call(".External",
                c(list(CFUN, x, y, pairwise,
                       if (!is.function(method)) get(method) else method),
                  params
                  )
                )

    result <-
### PASS-THROUGH-cases
        if (!is.null(reg_entry) && !reg_entry$loop) {
            if (!by_rows && !is.list(x)) {
                x <- t(x)
                if (!is.null(y))
                    y <- t(y)
            }
            if (reg_entry$C_FUN)
                do.call(".Call", c(list(method), list(x), list(y), pairwise, params,
                                   list(PACKAGE = reg_entry$PACKAGE)))
            else {   ## user functions need not implement pairwise
                if (!is.null(reg_entry$PACKAGE))
                    do.call(method, c(list(x), list(y), params), envir = asNamespace(reg_entry$PACKAGE))
                else
                    do.call(method, c(list(x), list(y), params))
            }
        } else if (is.null(y)) {
### LOOP WORKHORSE for auto-proximities
            ## transpose data for column-wise loop
            if (!by_rows && !is.list(x))
                x <- t(x)
            if (is.list(x) && !is.null(reg_entry) && reg_entry$abcd)
                x <- do.call("rbind", x)
            if (is.matrix(x) && !is.null(reg_entry) && reg_entry$abcd)
                ## binary matrix
                .proxy_external(R_apply_dist_binary_matrix, x != 0, NULL)
            else if (is.matrix(x))
                ## real, integer matrix
                .proxy_external(R_apply_dist_matrix, x, NULL)
            else if (is.list(x) && !(is.data.frame(x) && by_rows))
                ## list
                .proxy_external(R_apply_dist_list, x, NULL)
            else ## data frame (by rows)
                .proxy_external(R_apply_dist_data_frame, x, NULL)

        } else {
### LOOP WORKHORSE for cross-proximities
            ## transpose data for column-wise loop
            if (!by_rows && !is.list(x)) {
                x <- t(x)
                y <- t(y)
            }
            if (is.list(x) && !is.null(reg_entry) && reg_entry$abcd)
            {
                x <- do.call("rbind", x)
                y <- do.call("rbind", x)
            }
            if (is.matrix(x) && !is.null(reg_entry) && reg_entry$abcd)
                ## binary matrices
                .proxy_external(R_apply_dist_binary_matrix, x != 0, y != 0)
            else if (is.matrix(x))
                ## real, integer matrices
                .proxy_external(R_apply_dist_matrix, x, y)
            else if (is.list(x) && !(is.data.frame(x) && by_rows))
                ## lists
                .proxy_external(R_apply_dist_list, x, y)
            else ## data frames (by rows)
                .proxy_external(R_apply_dist_data_frame, x, y)

        }

### set col/rownames for cross-proximity-objects (if needed)
    if (is.matrix(result) && is.null(dimnames(result)))
        if (is.list(x) && !is.data.frame(x)) {
            rownames(result) <- names(x)
            colnames(result) <- names(y)
        } else if (by_rows) {
            rownames(result) <- rownames(x)
            colnames(result) <- rownames(y)
        } else {
            rownames(result) <- colnames(x)
            colnames(result) <- colnames(y)
        }

### POSTPROCESS
    if (!is.null(reg_entry)) {
        if (!is.na(reg_entry$POSTFUN))
            result <- do.call(reg_entry$POSTFUN, c(list(result, params)))
        if (!reg_entry$distance &&
            !(is.logical(convert_similarities) && !convert_similarities)) {
            result <- if (is.function(convert_similarities) ||
                          is.character(convert_similarities))
                do.call(convert_similarities, list(result))
            else if (is.null(reg_entry$convert))
                pr_simil2dist(result)
            else
                do.call(reg_entry$convert, list(result))
        }
        method <- reg_entry$names[1]
    }

### RETURN DIST-OBJECT
    result <-
        if (is.matrix(result))
            structure(result, class = "crossdist")
        else
        if (inherits(result, "dist"))
            structure(result, Diag = diag, Upper = upper)
        else
            structure(result, class = "pairdist")
    structure(result,
              method = if (is.character(method))
                            method
                       else
                          if (missing(method))
                               deparse(substitute(y))
                          else deparse(substitute(method)),
              call = match.call())
}

simil <-
function(x, y = NULL, method = NULL, ...,
         diag = FALSE, upper = FALSE, pairwise = FALSE,
         by_rows = TRUE, convert_distances = TRUE,
         auto_convert_data_frames = TRUE)
{
    ## convenience to allow dists(x, "method")
    if ((is.function(y) || is.character(y)) && is.null(method)) {
        method <- y
        y <- NULL
    }
    if (is.null(method))
        method <- if (is.data.frame(x))
            "Gower"
        else if (is.logical(x))
            "Jaccard"
        else
            "correlation"

    ret <- dist(x, y, method, ..., diag = diag, upper = upper, pairwise = pairwise,
                by_rows = by_rows, convert_similarities = FALSE,
                auto_convert_data_frames = auto_convert_data_frames)

    ## possibly convert to similarity
    reg_entry <- pr_DB$get_entry(attr(ret, "method"), stop_if_missing = FALSE)
    if (!is.null(reg_entry)) {
        if (reg_entry$distance &&
            !(is.logical(convert_distances) && !convert_distances)) {
            ret <- if (is.function(convert_distances) ||
                       is.character(convert_distances))
                do.call(convert_distances, list(ret))
            else if (is.null(reg_entry$convert))
                pr_simil2dist(ret)
            else
                do.call(reg_entry$convert, list(ret))
        }
    }

    class(ret) <- unique(c(if (inherits(ret, "crossdist")) "crosssimil" else "simil",
                           class(ret)))
    ret
}

as.matrix.simil <-
function (x, ...)
{
    ret <- NextMethod(x, ...)
    diag(ret) <- 1
    ret
}

# note that a simil object must always also be a dist
# object for method dispatch

as.simil <-
function(x, FUN = NULL)
{
    if (inherits(x, c("simil", "crosssimil")))
        x
    else if (inherits(x, c("dist", "crossdist"))) {
        class(x) <- if (inherits(x, "dist"))
            c("simil", class(x))
        else
            c("crosssimil", setdiff(class(x), "crossdist"))
        if (!is.null(FUN))
            FUN(x)
        else {
            reg_entry <- NULL
            if (!is.null(attr(x, "method")))
                reg_entry <- pr_DB$get_entry(attr(x, "method"),
                                             stop_if_missing = FALSE)
            if (!is.null(reg_entry) && !is.null(reg_entry$convert))
                do.call(reg_entry$convert, list(x))
            else
                pr_dist2simil(x)
        }
    } else
        structure(stats::as.dist(x), class = c("simil", "dist"))
}

as.dist <-
function(x, FUN = NULL)
{
    if (inherits(x, c("simil", "crosssimil"))) {
        class(x) <- if (inherits(x, "simil"))
            setdiff(class(x), "simil")
        else
            c("crossdist", setdiff(class(x), "crosssimil"))
        if (!is.null(FUN))
            FUN(x)
        else {
            reg_entry <- NULL
            if (!is.null(attr(x, "method")))
                reg_entry <- pr_DB$get_entry(attr(x, "method"),
                                             stop_if_missing = FALSE)
            if (!is.null(reg_entry) && !is.null(reg_entry$convert))
                do.call(reg_entry$convert, list(x))
            else
                pr_simil2dist(x)
        }
    } else if (inherits(x, c("dist", "crossdist")))
        x
    else
        stats::as.dist(x)
}

## as we do not know if the object is the result of some
## user-defined transformation the values of
## s(x,x) are not defined.

## we need to copy stats::as.matrix.dist() since the use of ::: is deprecated:

as.matrix <-
function(x, ...)
    base::as.matrix(x, ...)

as_matrix_dist <-
function (x, ...)
{
    size <- attr(x, "Size")
    df <- matrix(0, size, size)
    df[row(df) > col(df)] <- x
    df <- df + t(df)
    labels <- attr(x, "Labels")
    dimnames(df) <- if (is.null(labels))
        list(seq_len(size), seq_len(size))
    else list(labels, labels)
    df
}


as.matrix.simil <-
function(x, diag = NA, ...) {
    x <- as_matrix_dist(x)
    diag(x) <- diag
    x
}

## however, it seems reasonable to assume that d(x,x)=0,
## which is also the default in stats.

as.matrix.dist <-
function(x, diag = 0, ...) {
    x <- as_matrix_dist(x)
    diag(x) <- diag
    x
}

print.crossdist <-
print.crosssimil <-
function (x, digits = getOption("digits"),
          justify = "none", right = TRUE, ...)
{
    if (length(x) > 0) {
        m <- as.matrix(x)
        cf <- format(m, digits = digits, justify = justify)
        print(cf, quote = FALSE, right = right, ...)
    } else {
        cat(data.class(x), "(0)\n", sep = "")
    }
    invisible(x)
}

print.pairdist <-
function(x, ...)
{
    print(as.vector(x), ...)
    invisible(x)
}

pr_simil2dist <-
function(x)
    1 - x

pr_dist2simil <-
function(x)
    1 / (1 + x)

###
