
## also available in lattice, but not exported
is.characterOrExpression <- function (x) is.character(x) || is.expression(x)

.arrayIndices <- function(d, i)
    ## Suppose we have an array 'x' with dimension 'd'.  We can index
    ## 'x' in two different ways: x[i] or x[i_1, i_2, ..., i_d].
    ## Here, we are given 'i', and want to compute i_1, i_2, ..., i_d.
{
    ## Here's what we are doing:
    ## For length(d) == 3, note that (for 0-based indexing)
    ## 
    ## i1 = (i mod d[1])
    ## i2 = (i mod d[1] * d[2]) div d[1]
    ## i3 = (i mod d[1] * d[2] * d[3]) div d[1] * d[2]
    n <- length(d)
    ans <- vector(mode = "list", length = n)
    for (k in seq_along(ans))
    {
        ans[[k]] <- 1L + (((i-1L) %% prod(head(d, k)))
                          %/%
                          prod(head(d, k-1L)))
    }
    ans
}


combineLimits <-
    function(x, margin.x = 2L, margin.y = 1L,
             extend = TRUE, adjust.labels = TRUE)
{
    if (length(dim(x)) == 1L)
        warning("Only one conditioning variable; nothing interesting will happen.")
    indices <- .arrayIndices(dim(x), seq_len(prod(dim(x))))
    ## For regular `numeric' scales, all we need is to modify
    ## $[xy].scales.  But for `factor' scales, we need to leave
    ## $[xy].scales alone, and instead modify $[xy].used.at and $[xy].num.limit
    modifyLimits <- function(limits, margin, ext)
    {
        limits <- array(limits, dim = dim(x))
        for (i in seq_len(prod(dim(x))))
        {
            ## index.combine <- index.entry <- Rows(indices, i)
            index.combine <- Rows(indices, i)
            index.combine[margin] <- list(TRUE)
            ## limits[[i]] <-
            ##     range(do.call("[", c(list(limits), index.combine)), finite = TRUE)
            
            li <- unlist(do.call("[", c(list(limits), index.combine)))
            limits[[i]] <- if(all(is.na(li))) li else range(li, finite = TRUE)
        }
        if (ext) lapply(limits, lattice:::extend.limits)
        else limits
    }
    modifyUsed <- function(used.at, margin)
    {
        used.at <- array(used.at, dim = dim(x))
        for (i in seq_len(prod(dim(x))))
        {
            index.combine <- Rows(indices, i)
            index.combine[margin] <- list(TRUE)
            li <- unlist(do.call("[", c(list(used.at), index.combine)))
            used.at[[i]] <- sort(unique(li))
        }
        used.at
    }
    if (x$x.scales$relation != "free" && x$y.scales$relation != "free")
        warning("Function only has effect for scales with 'relation=\"free\"'.")
    if (x$x.scales$relation == "free" && length(margin.x))
    {
        if (is.characterOrExpression(x$x.limits[[1]]))
        {
            x$x.used.at <- modifyUsed(x$x.used.at, margin.x)
            x$x.num.limit <- modifyLimits(x$x.num.limit, margin.x, ext = FALSE)
        }
        else
            x$x.limits <- modifyLimits(x$x.limits, margin.x, ext = extend)
    }
    if (x$y.scales$relation == "free" && length(margin.y))
    {
        if (is.characterOrExpression(x$y.limits[[1]]))
        {
            x$y.used.at <- modifyUsed(x$y.used.at, margin.y)
            x$y.num.limit <- modifyLimits(x$y.num.limit, margin.y, ext = FALSE)
        }
        else
            x$y.limits <- modifyLimits(x$y.limits, margin.y, ext = extend)
    }
    if (adjust.labels)
    {
        ## Drop all but left/bottom-most labels, and set space to 0
        ## for those.  Needs to know layout, and will set it unless
        ## already set.
        npackets <- prod(dim(x))
        par.settings <- if (is.null(x$par.settings)) list() else x$par.settings
        if (is.null(x$layout))
            x$layout <-
                if (length(dim(x)) == 1L) c(dim(x), 1)
                else dim(x)[1:2]
        else if (!isTRUE(all.equal(x$layout[1:2], dim(x)[1:2])))
        {
            warning("'layout' does not match dimensions; displayed scales may be wrong.")
        }
        if (any(is.na(x$layout) | x$layout == 0))
            stop("'layout' must explicitly determine number of rows and columns")
        if (x$x.scales$relation == "free" && length(margin.x))
        {
            ## change x-scales
            if (is.list(x$x.scales$at))
            {
                warning("Explicit per-panel tick mark locations ignored")
                x$x.scales$at <- FALSE
            }
            page.at <- 
                if (x$as.table)
                    rep(list(NULL, x$x.scales$at),
                        c(x$layout[1] * (x$layout[2]-1), x$layout[1]))
                else
                    rep(list(x$x.scales$at, NULL),
                        c(x$layout[1], x$layout[1] * (x$layout[2]-1)))
            x$x.scales$at <- rep(page.at, length.out = npackets)
            par.settings <-
                if (x$as.table)
                    modifyList(par.settings,
                               list(layout.heights =
                                    list(axis.panel = rep(c(0, 1), c(x$layout[2]-1, 1)))))
                else 
                    modifyList(par.settings,
                               list(layout.heights =
                                    list(axis.panel = rep(c(1, 0), c(1, x$layout[2]-1)))))
        }
        if (x$y.scales$relation == "free" && length(margin.y))
        {
            ## change y-scales
            if (is.list(x$y.scales$at))
            {
                warning("Explicit per-panel tick mark locations ignored")
                x$y.scales$at <- FALSE
            }
            page.at <- rep(list(TRUE, NULL), c(1, x$layout[1]-1))
            x$y.scales$at <- rep(page.at, length.out = npackets)
            par.settings <-
                modifyList(par.settings,
                           list(layout.widths =
                                list(axis.panel = rep(c(1, 0), c(1, x$layout[1]-1)))))
        }
        x$par.settings <- par.settings
    }
    x
}



