### * plot function
plot.gset <-
function(x, type = NULL, ylim = NULL,
         xlab = "Universe", ylab = "Membership Grade", ...)
{
    if (cset_is_empty(x, na.rm = TRUE))
        return(invisible(x))

    if (is.null(type))
        type <- if(cset_is_set_or_fuzzy_set(x, na.rm = TRUE) &&
                   .domain_is_numeric(x))
            "l"
        else
            "barplot"

    if (!cset_is_fuzzy_multiset(x) && (type != "barplot"))
        return(plot(tuple(x), ...))

    if (is.null(ylim) && !cset_is_multiset(x, na.rm = TRUE))
        ylim <- c(0,1)

    m <- cset_memberships(x)
    if (cset_is_fuzzy_multiset(x)) {
        maxlen <- max(sapply(m, length), na.rm = TRUE)
        m <- sapply(m, .expand_membership, len = maxlen)
    }
    barplot(m,
            ylim = ylim,
            names.arg = LABELS(x),
            beside = TRUE,
            ...)

    invisible(x)
}

plot.cset <-
function(x, ...)
    plot.gset(x, ...)

plot.set <-
function(x, ...)
{
    if (all(sapply(x, is.cset)))
        plot(as.tuple(x), ...)
    else
        plot.gset(x, ...)
}

plot.tuple <-
function(x, type = "l", ylim = NULL,
         xlab = "Universe", ylab = "Membership Grade", col = 1,
         continuous = TRUE, ...)
{
    l <- as.list(x)

    ## expand charfun generators
    for (i in seq_along(l))
        if (is.charfun_generator(l[[i]]))
            l[[i]] <- gset(charfun = l[[i]](), universe = .expand(NULL))

    ## check for fuzzy multisets
    ind <- sapply(l, cset_is_fuzzy_multiset)
    if (all(ind))
        stop("Cannot plot tuple of fuzzy multisets.")
    if (any(ind)) {
        warning("All fuzzy multisets ignored.")
        l <- l[!ind]
    }

    ## ylim
    ylim <- if (is.null(ylim) &&
                !any(sapply(l, cset_is_multiset, na.rm = TRUE)))
        c(0, 1)
    else
        c(0, max(unlist(lapply(l, .get_memberships)), na.rm = TRUE))

    ## find common domain
    names(l) <- NULL ## remove labels
    universe <- do.call(set_union,
                        c(lapply(l, cset_support),
                          lapply(lapply(l, cset_universe), as.set))
                        )
    NUMMODE <- .domain_is_numeric(universe)

    ## prepare plot region
    m <- if (NUMMODE) {
        universe <- cset(universe, matchfun = match)
        sort(unique(unlist(universe)))
    } else
        seq_along(universe)

    plot(m, rep.int(0, length(universe)),
         ylim = ylim,
         type = "n",
         axes = FALSE,
         xlab = xlab,
         ylab = ylab,
         ...)
    axis(2)
    if (NUMMODE)
        axis(1, labels = TRUE)
    else
        axis(1, at = m, labels = LABELS(universe))

    ## plot sets
    lines(as.tuple(l), type = type, col = col, universe = universe, ...)
    invisible(x)
}

plot.charfun_generator <-
function(x, universe = NULL, ...)
    plot(gset(charfun = x(), universe = .expand(universe)), ...)

## lines

lines.tuple <-
function(x, col = 1, universe = NULL, ...)
{
    l <- as.list(x)
    col <- rep(col, length.out = length(x))

    ## expand charfun generators
    for (i in seq_along(l))
        if (is.charfun_generator(l[[i]]))
            l[[i]] <- gset(l[[i]](), universe = .expand(NULL))

    ## check for fuzzy multisets
    ind <- sapply(l, cset_is_fuzzy_multiset)
    if (all(ind))
        stop("Cannot plot tuple of fuzzy multisets.")
    if (any(ind)) {
        warning("All fuzzy multisets ignored.")
        l <- l[!ind]
    }

    ## try to deduce universe
    names(l) <- NULL ## remove labels
    universe <- do.call(set_union,
                        c(lapply(l, cset_support),
                          lapply(lapply(l, cset_universe), as.set))
                        )

    ## call workhorse
    for (i in seq_along(l))
        lines(l[[i]], col = col[i], universe = universe, ...)

    invisible(x)
}

lines.set <-
function(x, ...)
{
    if (all(sapply(x, is.cset)))
        lines(as.tuple(x), ...)
    else
        lines.gset(x, ...)
}

lines.cset <-
function(x, ...)
    lines.gset(x, ...)

lines.gset <-
function(x, type = "l", col = 1, continuous = TRUE, universe = NULL, ...)
{
    universe <- .expand(universe)
    if (.domain_is_numeric(universe)) {
        ndom <- sort(unique(unlist(cset(universe, matchfun = match))))
        y <- rep.int(0, length(ndom))
        s <- cset(x, matchfun = match)
        .match <- matchfun(function(x, y) isTRUE(all.equal(x, y)))
        matches <- .match(unlist(s), ndom)
        isna <- is.na(matches)
        if(all(isna)) return(invisible(x))
        y[matches[!isna]] <- cset_memberships(s)[!isna]

        ## remove discontinuities
        if (continuous) {
            sup <- range(which(y > 0))
            zeros <- (y == 0) &
            seq_along(y) >= sup[1] &
                  seq_along(y) <= sup[2]
            y <- y[!zeros]
            ndom <- ndom[!zeros]
        }

        ## make sure that first and last values of support
        ## connect straight to bottom
        if (!type %in% c("p", "h")) {
            .insert <-
                function(X, i, value)
                    c(X[1:(i-1)], value, X[i:length(X)])

            ind <- min(which(y > 0))
            if (ind > 1L) {
                y <- .insert(y, ind, 0)
                ndom <- .insert(ndom, ind, ndom[ind])
            }

            ind <- max(which(y > 0))
            if (ind < length(y)) {
                y <- .insert(y, ind + 1, 0)
                ndom <- .insert(ndom, ind + 1, ndom[ind])
            }
        }

        lines(ndom, y, col = col, type = type, ...)
    } else {
        lines(.exact_match(cset_support(x), as.list(universe)),
              cset_memberships(x),
              col = col,
              type = type,
              ...)
    }

    invisible(x)
}

lines.charfun_generator <-
function(x, universe = NULL, ...)
    lines(gset(charfun = x(), universe = .expand(universe)), ...)
