#
#  tree/R/treemisc.R Copyright (C) 1994-2013 B. D. Ripley
#  miscellaneous support routines for package 'tree'.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#


.noGenerics <- TRUE

.onUnload <- function(libpath) library.dynam.unload("tree", libpath)

tree_env <- new.env()

cv.tree <- function(object, rand, FUN = prune.tree, K = 10, ...)
{
    if(!inherits(object, "tree")) stop("not legitimate tree")
    m <- model.frame(object)
    extras <- match.call(expand.dots = FALSE)$...
    FUN <- deparse(substitute(FUN))
    init <- do.call(FUN, c(list(object), extras))
    if(missing(rand)) rand <- sample(K, length(m[[1L]]), replace = TRUE)
    cvdev <- 0
    for(i in unique(rand)) {
        tlearn <- tree(model = m[rand != i,  , drop = FALSE])
        plearn <- do.call(FUN, c(list(tlearn, newdata =
                                      m[rand == i, , drop = FALSE],
                                      k = init$k), extras))
        cvdev <- cvdev + plearn$dev
    }
    init$dev <- cvdev
    init
}

data.tree <- function(tree)
{
    oc <- tree$call
    while(oc[[1L]] != "tree") oc <- eval(oc[[2L]])$call
    if(is.null(oc$data)) NULL
    else {
        warning("retrieving data from ", deparse(oc$data))
        eval(oc$data)
    }
}

descendants <- function(nodes, include = TRUE)
{
    n <- length(nodes)
    if(n == 1L) return(matrix(TRUE, 1L, 1L))
    ind <- 1L:n
    desc <- matrix(FALSE, n, n)
    if(include) diag(desc) <- TRUE
    parents <- match((nodes %/% 2L), nodes)
    lev <- floor(log(nodes, base = 2))
    desc[1L, 2L:n] <- TRUE
    for(i in max(lev):2L) {
        desc[cbind(ind[parents[lev == i]], ind[lev == i])] <- TRUE
        parents[lev == i] <- parents[parents[lev == i]]
        lev[lev == i] <- i - 1L
    }
    desc
}

deviance.tree <- function(object, detail = FALSE, ...)
{
    if(!inherits(object, "tree")) stop("not legitimate tree")
    frame <- object$frame
    if(detail) frame$dev
    else sum(frame$dev[frame$var == "<leaf>"])
}

labels.tree <- function(object, pretty = TRUE, collapse = TRUE, ...)
{
    if(!inherits(object, "tree")) stop("not legitimate tree")
    frame <- object$frame
    xlevels <- attr(object, "xlevels")
    var <- as.character(frame$var)
    splits <- matrix(sub("^>", " > ", sub("^<", " < ", frame$splits)),, 2L)
    lt <- c(letters, 0:5) # max 32 levels
    if(!is.null(pretty)) {
        if(pretty) xlevels <- lapply(xlevels, abbreviate, minlength=pretty)
        for(i in grep("^:", splits[, 1L],))
            for(j in 1L:2L) {
                sh <- splits[i, j]
                nc <- nchar(sh)
                sh <- substring(sh, 2L:nc, 2L:nc)
                xl <- xlevels[[var[i]]][match(sh, lt)]
                splits[i, j] <- paste0(": ", paste(as.vector(xl), collapse=","))

            }
    }
    if(!collapse) return(array(paste0(var, splits), dim(splits)))
    node <- as.integer(row.names(frame))
    parent <- match((node %/% 2L), node)
    odd <- as.logical(node %% 2L)
    node[odd] <- paste0(var[parent[odd]], splits[parent[odd], 2L])
    node[!odd] <- paste0(var[parent[!odd]], splits[parent[!odd], 1L])
    node[1L] <- "root"
    node
}

misclass.tree <- function(tree, detail = FALSE)
{
    if(!inherits(tree, "tree")) stop("not legitimate tree")
    if(is.null(attr(tree, "ylevels")))
        stop("misclassification error rate is appropriate for factor responses only")
    if(is.null(y <- tree$y)) y <- model.response(model.frame(tree))
    if(is.null(wts <- tree$weights)) wts <- model.weights(model.frame(tree))
    if(is.null(wts)) wts <- rep(1, length(y))
    frame <- tree$frame
    if(detail) {
        which <- descendants(as.integer(row.names(frame)))
        tmp <- as.vector((which[, tree$where] *
                          outer(frame$yval, y, "!=")) %*% wts)
        names(tmp) <- row.names(tree$frame)
        tmp
    }
    else sum(wts*(frame$yval[tree$where] != y))
}

model.frame.tree <- function(formula, ...)
{
    if(!is.null(m <- formula$model)) return(m)
    oc <- formula$call
    if(substring(deparse(oc[[1L]]), 1L, 7L) == "predict") {
        m <- eval(oc$newdata)
        if(is.null(attr(m, "terms"))) {
            object <- eval(oc$object)
            m <- model.frame(object$terms, m, na.pass)
        }
        return(m)
    }
    while(!deparse(oc[[1L]]) %in% c("tree", "tree::tree", "tree:::tree"))
        oc <- eval(oc[[2L]])$call
    oc$subset <- names(formula$where)
    oc$method <- "model.frame"
    eval(oc)
}

node.match <- function(nodes, nodelist, leaves, print.it = TRUE)
{
    node.index <- match(nodes, nodelist, nomatch = 0L)
    bad <- nodes[node.index == 0L]
    if(length(bad) > 0L & print.it)
        warning("supplied nodes ", paste(bad, collapse = ","),
                " are not in this tree")
    good <- nodes[node.index > 0L]
    if(!missing(leaves) && any(leaves <- leaves[node.index])) {
        warning("supplied nodes ",
                paste(good[leaves], collapse = ","), " are leaves")
        node.index[node.index > 0L][!leaves]
    }
    else node.index[node.index > 0L]
}

partition.tree <- function(tree, label = "yval", add = FALSE, ordvars, ...)
{
    ptXlines <- function(x, v, xrange, xcoord = NULL, ycoord = NULL, tvar, i = 1L)
    {
        if(v[i] == "<leaf>") {
            y1 <- (xrange[1L] + xrange[3L])/2
            y2 <- (xrange[2L] + xrange[4L])/2
            return(list(xcoord = xcoord, ycoord = c(ycoord, y1, y2), i = i))
        }
        if(v[i] == tvar[1L]) {
            xcoord <- c(xcoord, x[i], xrange[2L], x[i], xrange[4L])
            xr <- xrange
            xr[3L] <- x[i]
            ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i + 1L)
            xr <- xrange
            xr[1L] <- x[i]
            return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar, ll2$i + 1L))
        } else if(v[i] == tvar[2L]) {
            xcoord <- c(xcoord, xrange[1L], x[i], xrange[3L], x[i])
            xr <- xrange
            xr[4L] <- x[i]
            ll2 <- Recall(x, v, xr, xcoord, ycoord, tvar, i + 1L)
            xr <- xrange
            xr[2L] <- x[i]
            return(Recall(x, v, xr, ll2$xcoord, ll2$ycoord, tvar, ll2$i + 1L))
        }
        else stop("wrong variable numbers in tree.")
    }
    if(inherits(tree, "singlenode")) stop("cannot plot singlenode tree")
    if(!inherits(tree, "tree")) stop("not legitimate tree")
    frame <- tree$frame
    leaves <- frame$var == "<leaf>"
    var <- unique(as.character(frame$var[!leaves]))
    if(length(var) > 2L || length(var) < 1L)
        stop("tree can only have one or two predictors")
    nlevels <- sapply(attr(tree, "xlevels"), length)
    if(any(nlevels[var] > 0L))
        stop("tree can only have continuous predictors")
    x <- rep(NA, length(leaves))
    x[!leaves] <- as.double(substring(frame$splits[!leaves, "cutleft"],
                                      2L, 100L))
    m <- model.frame(tree)
    if(length(var) == 1L) {
        ## one x variable
        x <- sort(c(range(m[[var]]), x[!leaves]))
        if(is.null(attr(tree, "ylevels"))) y <- frame$yval[leaves]
        else y <- frame$yprob[, 1L]
        y <- c(y, y[length(y)])
        if(add) lines(x, y, type = "s", ...)
        else {
            a <- attributes(attr(m, "terms"))
            yvar <- as.character(a$variables[1+a$response])
            xo <- m[[yvar]]
            if(is.factor(xo)) ylim <- c(0,1) else ylim <- range(xo)
            plot(x, y, ylab = yvar, xlab = var, type = "s", ylim = ylim,
                 xaxs = "i", ...)
        }
        invisible(list(x = x, y = y))
    } else {
        ## two x variables
        if(!missing(ordvars)) {
            ind <- match(var, ordvars)
            if(any(is.na(ind))) stop("unmatched names in vars")
            var <- ordvars[sort(ind)]
        }
        lab <- frame$yval[leaves]
        if(is.null(frame$yprob)) lab <- format(signif(lab, 3L))
        else if(match(label, attr(tree, "ylevels"), nomatch = 0L))
            lab <- format(signif(frame$yprob[leaves, label], 3L))
        rx <- range(m[[var[1L]]])
        rx <- rx + c(-0.025, 0.025) * diff(rx)
        rz <- range(m[[var[2L]]])
        rz <- rz + c(-0.025, 0.025) * diff(rz)
        xrange <- c(rx, rz)[c(1, 3, 2, 4)]
        xcoord <- NULL                  # x1lo, x2lo, x1hi, x2hi
        ycoord <- NULL                  # y1, y2
        xy <- ptXlines(x, frame$var, xrange, xcoord, ycoord, var)
        xx <- matrix(xy$xcoord, nrow = 4L)
        yy <- matrix(xy$ycoord, nrow = 2L)
        if(!add)
            plot(rx, rz, xlab = var[1L], ylab = var[2L], type = "n",
                 xaxs = "i", yaxs = "i", ...)
        segments(xx[1L,  ], xx[2L,  ], xx[3L,  ], xx[4L,  ])
        text(yy[1L,  ], yy[2L,  ], as.character(lab), ...)
    }
}

plot.tree <- function (x, y = NULL,
                       type = c("proportional", "uniform"), ...)
{
    if(inherits(x, "singlenode")) stop("cannot plot singlenode tree")
    if(!inherits(x, "tree")) stop("not legitimate tree")
    type <- match.arg(type)
    uniform <- type == "uniform"
    dev <- dev.cur()
    if (dev == 1L) dev <- 2L # as device will be opened.
    assign(paste0("device", dev), uniform, envir = tree_env)
    invisible(treepl(treeco(x, uniform),
                     node = as.integer(row.names(x$frame)), ...))
}

print.tree <-
    function(x, pretty = 0, spaces = 2, digits = getOption("digits")-3, ...)
{
    if(!inherits(x, "tree")) stop("not legitimate tree")
    is.prob <- !is.null(ylevels <- attr(x, "ylevels"))
    if(is.prob) cat("node), split, n, deviance, yval, (yprob)\n")
    else cat("node), split, n, deviance, yval\n")
    cat("      * denotes terminal node\n\n")
    frame <- x$frame
    node <- as.integer(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = "")
                                     #32 is the maximal depth
    if(length(node) > 1L) {
        indent <- substring(indent, 1L, spaces * seq(depth))
        indent <- paste0(c("", indent[depth]), format(node), ")")
    } else
    indent <- paste0(format(node), ")")
    if(is.prob) {
        yval <- paste0(as.character(frame$yval), " (")
        yprob <- format(frame$yprob, digits = digits)
        for(i in seq(ylevels)) yval <- paste(yval, yprob[, i])
        yval <- paste(yval, ")")
    } else
    yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels.tree(x, pretty = pretty)
    z <- paste(indent, z, round(frame$n, 2L),
               format(signif(frame$dev, digits = digits)),
               yval, term)
    cat(z, sep = "\n")
    invisible(x)
}

residuals.tree <-
    function(object, type = c("usual", "pearson", "deviance"), ...)
{
    if(!inherits(object, "tree")) stop("not legitimate tree")
    if(is.null(y <- object$y))
        y <- model.extract(model.frame(object), "response")
    frame <- object$frame
    if(is.null(attr(object, "ylevels"))) return(y - frame$yval[object$where])
    type <- match.arg(type)
    if(type == "usual") yhat <- frame$yval[object$where]
    else yhat <- frame$yprob[object$where,  ][cbind(seq(y), unclass(y))]
    r <- switch(type,
                usual = as.integer(y != yhat),
                pearson = (1 - yhat)/yhat,
                deviance = -2 * log(yhat))
    names(r) <- names(y)
    r
}

snip.tree <-
    function(tree, nodes, xy.save = FALSE, digits = getOption("digits") - 3)
{
    .snip.tree <- function(tree, nodes)
    {
        where <- tree$where
        frame <- tree$frame
        node <- as.integer(row.names(frame))
        if(is.null(frame$which)) frame$which <- descendants(node, FALSE)
        i <- match(nodes, node)
        frame$var[i] <- "<leaf>"
        frame$splits[i,  ] <- ""
        ii <- frame$which[i,  , drop = FALSE]
        keep <- !apply(ii, 2L, any)
        frame <- frame[keep,  , drop = FALSE]
        frame$which <- frame$which[, keep, drop = FALSE]
        where <- match(node[where], node[keep])
        nn <- match(nodes, node[keep])
        for(n in seq(length(nodes))) where[ii[n, tree$where]] <- nn[n]
        tree$frame <- frame
        tree$where <- structure(where, names = names(tree$where))
        tree
    }

    if(inherits(tree, "singlenode")) stop("cannot snip singlenode tree")
    if(!inherits(tree, "tree")) stop("not legitimate tree")
    call <- match.call()
    node <- as.integer(row.names(tree$frame))
    if(missing(nodes)) {
        iprev <- i <- 0L
        nodes <- NULL
        totdev <- deviance(tree)
        xy <- treeco(tree)
        tree$frame$which <- descendants(node, include = FALSE)
        labs <- c("node number: ", "  tree deviance = ", "  subtree deviance = ")
        while(length(i <- identify(xy$x, xy$y, n=1, plot=FALSE)) > 0L) {
            if(tree$frame$var[i] == "<leaf>") {
                cat("Terminal node -- try again\n")
                next
            }
            ii <- tree$frame$which[i,  ]
            if(i != iprev) {
                frame <- tree$frame
                dev <- frame$dev
                node <- as.integer(row.names(frame))
                newdev <- totdev - sum(dev[ii & (frame$var == "<leaf>")]) + dev[i]
                stats <- c(node[i], format(signif(c(totdev, newdev), digits)))
                cat(paste(labs, stats, "\n"))
                iprev <- i
            } else {
                xy <- treepl(xy, node, erase = ii)
                nodes <- c(nodes, node[i])
                tree <- .snip.tree(tree, node[i])
                totdev <- newdev
            }
        }
        if(xy.save) attr(tree, ".xy") <- xy
    } else {
        i <- node.match(nodes, node, tree$frame$var == "<leaf>")
        if(length(i) == 0L) return(tree)
        nodes <- node[i]
        which <- descendants(node, include = TRUE)
        nodes <- nodes[colSums(which[i,  ] %*% which[, i]) == 1]
        diag(which) <- FALSE
        tree$frame$which <- which
        tree <- .snip.tree(tree, nodes)
    }
    tree$frame$which <- NULL
    call$nodes <- nodes
    tree$call <- call
    if(dim(tree$frame)[1L] == 1L) class(tree) <- "singlenode"
    tree
}

summary.tree <- function(object, ...)
{
    obj <- list(call = object$call)
    frame <- object$frame
    obj$type <- if(is.reg <- is.null(attr(object, "ylevels")))
        "\nRegression tree:\n" else "\nClassification tree:\n"
    leaves <- frame$var == "<leaf>"
    variables <- names(attr(object, "xlevels"))
    used <- unique(frame$var[!leaves])
    if(length(used) < length(variables)) obj$used <- used
    obj$size <- sum(leaves)
    obj$df <- frame$n[1L] - obj$size
    obj$dev <- deviance.tree(object)
    if(!is.reg) obj$misclass <- c(misclass.tree(object), frame$n[1L])
    else obj$residuals <- residuals(object)
    class(obj) <- "summary.tree"
    obj
}

print.summary.tree <-
    function(x, digits = max(getOption("digits") - 3, 3), ...)
{
    cat(x$type)
    dput(x$call)
    if(!is.null(x$used)) {
        cat("Variables actually used in tree construction:\n")
        print(as.character(x$used))
    }
    cat(paste("Number of terminal nodes: ", x$size, "\n"))
    if(!is.null(x$effect.size))
        cat("Effective number of terminal nodes: ",
            format(signif(x$ effect.size, digits)), "\n")
    cat(paste("Residual mean deviance: ",
              format(signif(x$dev/x$df, digits)), "=",
              format(signif(x$dev, digits)), "/",
              format(signif(x$df, digits)), "\n"))
    if(!is.null(x$misclass))
        cat("Misclassification error rate:",
            format(signif(x$misclass[1L]/x$misclass[2L], digits)),
            "=", x$misclass[1L], "/", x$ misclass[2L], "\n")
    else {
        cat("Distribution of residuals:\n")
        print(summary(x$residuals, digits = digits, ...))
    }
    invisible(x)
}

text.tree <-
    function(x, splits = TRUE, label = "yval", all = FALSE,
             pretty = NULL, digits = getOption("digits") - 3,
             adj = par("adj"), xpd = TRUE, ...)
{
    oldxpd <- par(xpd=xpd)
    on.exit(par(oldxpd))
    if(inherits(x, "singlenode")) stop("cannot plot singlenode tree")
    if(!inherits(x, "tree")) stop("not legitimate tree")
    frame <- x$frame
    column <- names(frame)
    if(!is.null(ylevels <- attr(x, "ylevels"))) column <- c(column, ylevels)
    if(!is.null(label) && is.na(match(label, column)))
        stop("label must be a column label of the frame component of the tree")
    charht <- par("cxy")[2L]
    if(!is.null(srt <- list(...)$srt) && srt == 90) {
        if(missing(adj)) adj <- 0
        ladj <- 1 - adj
    } else ladj <- adj
    xy <- treeco(x)
    if(splits) {
        node <- as.integer(row.names(frame))
        left.child <- match(2 * node, node)
        rows <- labels.tree(x, pretty = pretty)[left.child]
        ind <- !is.na(rows)
        text(xy$x[ind], xy$y[ind] + 0.5 * charht, rows[ind], adj=adj,  ...)
    }
    if(!is.null(label)) {
        leaves <- if(all) rep(TRUE, nrow(frame)) else frame$var == "<leaf>"
        if(label == "yval" & !is.null(ylevels))
            stat <- as.character(frame$yval[leaves])
        else if(!is.null(ylevels) && !is.na(lev <- match(label, ylevels)))
            stat <- format(signif(frame$yprob[leaves, lev], digits = digits))
        else stat <- format(signif(frame[leaves, label], digits = digits))
        if(!is.null(dim(stat)) && dim(stat)[2L] > 1) {
            if(length(dimnames(stat)[[2L]]))
                stat[1L,  ] <- paste(sep = ":", dimnames(stat)[[2L]], stat[1L,  ])
            stat <- do.call("paste",
                            c(list(sep = "\n"), split(stat, col(stat))))
        }
        text(xy$x[leaves], xy$y[leaves] - 0.5 * charht, labels = stat,
             adj = ladj, ...)
    }
    invisible()
}

tile.tree <- function(tree, var, screen.arg = ascr + 1, axes = TRUE)
{
    if(inherits(tree, "singlenode")) stop("cannot tile singlenode tree")
    if(!inherits(tree, "tree")) stop("not legitimate tree")
    where <- tree$where
    varname <- substitute(var)
    v <- match(deparse(varname), as.character(attr(tree$terms, "variables")),
               nomatch = 0L)
    if(v) var <- model.frame(tree)[, v]
    else {
        if(!is.null(data <- data.tree(tree))) {
            var <- eval(varname, data)
            names(var) <- row.names(data)
        } # else names(var) <- database.attr("row.names")
        var <- var[names(where)]
    }
    if(length(var) != length(where))
        stop("variable length must match that of data used in fit")
    if(any(is.na(var))) var <- na.tree.replace(list(var = var))$var
    if(!(levx <- length(levels(var)))) {
        var <- cut(var, quantile(var) - c(1, 0, 0, 0, 0))
        levx <- 4L
        names(var) <- names(where)
    }
    leaves <- tree$frame$var == "<leaf>"
    nl <- sum(leaves)
    x <- seq(nl)
    indices <- unlist(split(names(where), where))
    counts <- c(table(var[indices], rep(1L:nl, tree$frame$n[leaves])))
    xm <- max(x)
    y <- seq(0, 1, length = levx + 1L)[-1L]
    dy <- y[2L] - y[1L]
    y <- rep(y, nl)[as.logical(counts)]
    x <- rep(x, rep(levx, nl))[as.logical(counts)]
    dx <- counts[as.logical(counts)]/(max(counts) * 2)
    x <- rbind(x, x - dx, x - dx, x, x, NA)
    y <- rbind(y, y, y + dy, y + dy, y, NA)
    if(!(ascr <- screen())) stop("need to use split.screen first")
    screen(screen.arg)
    xpd <- par(xpd = TRUE, mar=rep(0,4))
    on.exit(par(xpd))
    plot(x, y, type = "l", axes = FALSE, xlab = "", ylab = "", xlim = c(1, xm))
    if(axes)
        abline(v = seq(xm))
    screen(ascr, new = FALSE)
    counts <- matrix(counts, levx, nl,
                     dimnames = list(levels(var), row.names(tree$frame)[leaves]))
    invisible(counts)
}

tree.control <- function(nobs, mincut = 5, minsize = 10, mindev = 0.01)
{
    mcut <- missing(mincut)
    msize <- missing(minsize)
    if(mincut > (minsize/2)) {
        if(mcut) mincut <- trunc(minsize/2)
        else if(msize && !mcut) minsize <- 2 * mincut
        else if(!msize && !mcut)
            stop("mincut cannot be greater than minsize/2")
    }
    mincut <- max(1, mincut)
    minsize <- max(2, minsize)
    nmax <- ceiling((4 * nobs)/(minsize - 1))
    list(mincut = mincut, minsize = minsize, mindev = mindev, nmax = nmax,
         nobs = nobs)
}

tree.depth <- function(nodes)
{
    depth <- floor(log(nodes, base = 2) + 1e-7)
    as.vector(depth - min(depth))
}

tree.matrix <- function(frame)
{
    if(!inherits(frame, "data.frame")) return(as.matrix(frame))
    frame$"(weights)" <- NULL
    terms <- attr(frame, "terms")
    if(is.null(terms)) predictors <- names(frame)
    else predictors <- as.character(attr(terms, "term.labels"))
    frame <- frame[predictors]
    term.levels <- lapply(frame, levels)
    factors <- sapply(term.levels, function(x) length(x) > 0L)
    if(any(factors)) {
        for(preds in predictors[factors])
            frame[[preds]] <- as.vector(unclass(frame[[preds]]))
        x <- as.matrix(frame)
        column.levels <- vector("list", ncol(x))
        names(column.levels) <- dimnames(x)[[2L]]
        TT <- term.levels[factors > 0L]
        column.levels[names(TT)] <- TT
        attr(x, "column.levels") <- column.levels
    }
    else x <- as.matrix(frame)
    class(x) <- "matrix"
    x
}

tree.screens <- function(figs, screen.arg = 0, ...)
{
    if(missing(figs))
        figs <- matrix(c(0, 0, 1, 1, 0.25, 0.05, 1, 0.25), 2L, 4L)
    which <- split.screen(figs, screen.arg)
    par(mar = rep(0, 4), cex = 1, ...)
    which
}

treeco <-
    function(tree, uniform)
{
    if(missing(uniform)) {
        pn <- paste0("device", dev.cur())
        uniform <- if(exists(pn, envir = tree_env, inherits = FALSE))
            get(pn, envir = tree_env, inherits = FALSE)
        else FALSE
    }

    frame <- tree$frame
    node <- as.integer(row.names(frame))
    depth <- tree.depth(node)
    x <- -depth
    if(uniform) y <- x
    else {
        y <- dev <- frame$dev
        depth <- split(seq(node), depth)
        parent <- match(node %/% 2L, node)
        sibling <- match(ifelse(node %% 2L, node - 1L, node + 1L), node)
        for(i in depth[-1L])
            y[i] <- y[parent[i]] - dev[parent[i]] + dev[i] + dev[sibling[i]]
    }
    depth <- -x # seems unneeded
    leaves <- frame$var == "<leaf>"
    x[leaves] <- seq(sum(leaves))
    depth <- split(seq(node)[!leaves], depth[!leaves])
    left.child <- match(node * 2L, node)
    right.child <- match(node * 2 + 1L, node)
    for(i in rev(depth)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
    list(x = x, y = y)
}

treepl <- function(xy, node, erase = FALSE, ...)
{
    x <- xy$x; y <- xy$y
    parent <- match((node %/% 2L), node)
    sibling <- match(ifelse(node %% 2L, node - 1L, node + 1L), node)
    xx <- rbind(x, x, x[sibling], x[sibling], NA)
    yy <- rbind(y, y[parent], y[parent], y[sibling], NA)
    if(any(erase)) {
        lines(c(xx[, erase]), c(yy[, erase]), col = par("bg"))
        return(list(x = x[!erase], y = y[!erase]))
    }
    plot(range(x), range(y), type = "n", axes = FALSE, xlab = "", ylab = "")
    text(x[1L], y[1L], "|", ...)
    lines(c(xx[, -1L]), c(yy[, -1L]), ...)
    list(x = x, y = y)
}
