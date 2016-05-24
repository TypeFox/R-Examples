### File R/plot.sisal.R
### This file is part of the sisal package for R.
###
### Copyright (C) 2015 Aalto University
###
### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### A copy of the GNU General Public License is available at
### http://www.r-project.org/Licenses/

graphvizInstructions <- function() {
    cat(gettext("Installation instructions:",
                "1. Source the biocLite script", domain="R-sisal"),
        '  source("http://www.bioconductor.org/biocLite.R")',
        gettext("2. Check URL of current Rgraphviz README",
                domain="R-sisal"),
        paste('  readme <- file.path(biocinstallRepos()["BioCsoft"],',
              '"readmes", "Rgraphviz", "README", fsep="/")'),
        gettext("3. View README (pressing 'q' is a likely exit route)",
                domain="R-sisal"),
        "  url.show(readme)",
        gettext("4a. If Rgraphviz comes bundled with Graphviz",
                domain="R-sisal"),
        '  biocLite("Rgraphviz")',
        gettext("4b. Otherwise, follow instructions in README",
                domain="R-sisal"), sep="\n")
}

plot.sisal <- function(x, which = 1, standardize = "inherit", ...,
                       plotArgs = list(list(), list(mai = rep(0.1, 4))),
                       xlim = c(x[["d"]], 0), ylim = NULL,
                       ask = TRUE, dev.set = !ask,
                       draw.node.labels = TRUE, draw.edge.labels = TRUE,
                       draw.selected.labels = TRUE,
                       rankdir = c("TB", "LR", "BT", "RL"),
                       fillcolor.normal = "deepskyblue",
                       fillcolor.pruned = "deeppink",
                       fillcolor.selected = "chartreuse",
                       fillcolor.levelbest = "gold",
                       fillcolor.small = "moccasin",
                       fillcolor.large = "black",
                       fillcolor.NA = "white",
                       bordercolor.normal = "black",
                       bordercolor.special.levelbest = fillcolor.levelbest,
                       bordercolor.special.selected = fillcolor.selected,
                       color.by.error = FALSE,
                       ramp.space = c("Lab", "rgb"),
                       ramp.size = 128,
                       error.limits = c(NA_real_, NA_real_),
                       category.labels = c(normal=gettext("Other",
                                           domain="R-sisal"),
                       pruned=gettext("Pruned", domain="R-sisal"),
                       levelbest=gettext("Best\nin class", domain="R-sisal"),
                       selected=gettext("Selected", domain="R-sisal"),
                       special.levelbest=gettext("Best\n(no branching)",
                       domain="R-sisal"),
                       special.selected=gettext("Selected\n(no branching)",
                       domain="R-sisal"),
                       shape.normal=gettext("Other", domain="R-sisal"),
                       shape.highlighted=gettext("Highlighted",
                       domain="R-sisal")),
                       integrate.colorkey = TRUE, colorkey.gap = 0.1,
                       colorkey.space = c("right", "bottom", "left", "top"),
                       colorkey.title.gp = gpar(fontface = "bold"),
                       nodesep = 0.25, ranksep = 0.5,
                       graph.attributes = character(0),
                       node.attributes = character(0),
                       edge.attributes = character(0)) {
    if (!inherits(x, "sisal")) {
        stop('use only with "sisal" objects')
    }
    MAX.WHICH <- 3

    ## Conversion from RGB to brightness / luminance / luma...
    ## Sources referenced on 2015-10-09.
    ## http://stackoverflow.com/q/596216
    ## https://en.wikipedia.org/wiki/Luma_%28video%29
    LUMA.CONVERSION <- c(0.2126, 0.7152, 0.0722) # R, G, B

### Checks
    rankdir2 <- match.arg(rankdir)
    colorkey.space2 <- match.arg(colorkey.space)

    if (color.by.error) {
        ## Check color ramp settings
        ramp.space2 <- match.arg(ramp.space)
        stopifnot(is.numeric(ramp.size), length(ramp.size) == 1,
                  round(ramp.size) == ramp.size, ramp.size >= 2)
    } else {
        ## Check category labels
        stopifnot(is.character(category.labels))
        cl.names <- names(category.labels)
        required.names <- c("normal", "pruned", "selected", "levelbest",
                            "special.levelbest", "special.selected",
                            "shape.normal", "shape.highlighted")
        if (!all(required.names %in% cl.names)) {
            stop(gettextf("names required in 'category.labels'",
                          paste0(required.names, collapse=", "),
                          domain="R-sisal"), domain=NA)
        }
    }

    ## Check which
    if (!(is.numeric(which) || length(which) == 0) ||
        !all(which %in% 1:MAX.WHICH)) {
        stop(gettextf("'%s' must be a numeric vector with all elements in %.0f:%.0f",
                      "which", 1, MAX.WHICH, domain="R-sisal"),
             domain=NA)
    }

    ## Check color parameters (must be valid R colors)
    color.pars <- ls(pattern="^fillcolor[.]")
    for (the.name in color.pars) {
        the.value <- get(the.name, inherits=FALSE)
        if (inherits(try(col2rgb(the.value), silent=FALSE), "try-error")) {
            stop(gettextf("Unknown R color in '%s'. See ?col2rgb.",
                          the.name, domain = "R-sisal"),
                 domain=NA)
        }
    }

    ## Check logical parameters
    logical.pars <- c("ask", "dev.set", "color.by.error",
                      "draw.node.labels", "draw.edge.labels",
                      "draw.selected.labels", "integrate.colorkey")
    na.forbidden <- logical.pars
    for (the.name in logical.pars) {
        the.value <- get(the.name, inherits=FALSE)
        if (!is.logical(the.value) || length(the.value) != 1) {
            stop(gettextf("'%s' must have a logical value", the.name,
                          domain="R-sisal"),
                 domain=NA)
        }
    }
    for (the.name in na.forbidden) {
        the.value <- get(the.name, inherits=FALSE)
        if (is.na(the.value)) {
            stop(gettextf("'%s' must not be NA", the.name, domain="R-sisal"),
                 domain=NA)
        }
    }

    if (!identical(standardize, "inherit") &&
        !identical(standardize, TRUE) && !identical(standardize, FALSE)) {
        stop("'standardize' must be \"inherit\", TRUE or FALSE")
    }
    orig.standardize <- x[["params"]][["standardize"]]
    if (is.character(standardize)) {
        standardize2 <- orig.standardize
    } else {
        standardize2 <- standardize
    }

    ## Check that required packages are available
    show <- rep.int(FALSE, MAX.WHICH)
    show[which] <- TRUE
    if (show[2]) {
        req.gr <- FALSE
        req.Rgv <- FALSE
        if (inherits(try(suppressWarnings(req.Rgv <-
                                          requireNamespace("Rgraphviz")),
                         silent = TRUE),
                     "try-error") || !req.Rgv) {
            warning(gettextf("package '%s' is needed but not installed",
                             "Rgraphviz", domain="R-sisal"),
                    domain=NA)
            graphvizInstructions()
        } else if (inherits(try(suppressWarnings(req.gr <-
                                                 requireNamespace("graph")),
                                silent = TRUE),
                            "try-error") || !req.gr) {
            warning(gettextf("package '%s' is needed but not installed",
                             "graph", domain="R-sisal"),
                    domain=NA)
        }
        if (!req.Rgv || !req.gr) {
            warning("not doing plot number 2")
            show[2] <- FALSE
        }
    }

    ## Check xlim
    stopifnot(is.numeric(xlim), length(xlim) == 2)

    ## Check ylim
    if (!is.null(ylim) && (!is.numeric(ylim) || length(ylim) != 2)) {
        stop("'ylim' must be NULL or a numeric vector with 2 elements")
    }

    ## Check error.limits
    stopifnot(is.numeric(error.limits), length(error.limits) == 2)

    ## Check if the plot method of Ragraph in Rgraphviz supports the 'mai'
    ## argument required by (this implementation of) 'integrate.colorkey'.
    ## Rgraphviz 2.5 or later required.
    if (all(show[c(2, 3)]) && integrate.colorkey &&
        packageVersion("Rgraphviz") < "2.5") {
        warning("'integrate.colorkey' requires Rgraphviz version >= \"2.5\"")
        integrate.colorkey2 <- FALSE
    } else {
        integrate.colorkey2 <- integrate.colorkey
    }

    make.plots <- any(show)
    res <- NULL
    d <- x[["d"]]
    D <- d + 1

### Prepare for plotting
    if (show[2]) {
        ## Make graph
        g <- graph::graphNEL(nodes=x[["vertices"]], edgeL=x[["edges"]],
                             edgemode="directed")
        res <- g
        ## Check numeric parameters (they must all be positive)
        numeric.pars <- c("nodesep", "ranksep")
        for (the.name in numeric.pars) {
            the.value <- get(the.name, inherits=FALSE)
            if (!is.numeric(the.value) || length(the.value) != 1 ||
                the.value <= 0) {
                stop(gettextf("'%s' must have a positive numeric value",
                              the.name, domain="R-sisal"), domain=NA)
            }
        }

        ## Node shape (circle, box), colors (fill, border) and border width
        oD <- graph::degree(g)[["outDegree"]]
        nN <- names(oD)
        n.nodes <- length(nN)
        n.vars.in.node <- sapply(strsplit(nN, ".", fixed=TRUE), length)
        n.vars.in.node[nN == vars.to.name(empty = TRUE)] <- 0
        L.v <- x[["L.v"]]
        L.f <- x[["L.f"]]
        L.v.name <- paste0(L.v, collapse=".")
        L.f.name <- paste0(L.f, collapse=".")
        L.v.loc <- which(nN == L.v.name)
        L.f.loc <- which(nN == L.f.name)
        vertex.data <- x[["vertex.data"]][nN, ]
        bordercolor <- rep.int(rgb(t(col2rgb(bordercolor.normal)),
                                   maxColorValue=255), n.nodes)
        names(bordercolor) <- nN
        nodeshape <- rep.int("box", n.nodes)
        names(nodeshape) <- nN
        which.pruned <- which(oD == 0 & n.vars.in.node > 0)
        n.pruned <- length(which.pruned)
        levelbest.loc <- which(vertex.data[["E.v.level.rank"]] == 1)
        c.good <- c(levelbest.loc, L.v.loc, L.f.loc)
        nodeshape[c.good] <- "circle"
        if (color.by.error) {
            E.v <- vertex.data[["E.v"]]
            min.E.v <- min(E.v, na.rm=TRUE)
            max.E.v <- max(E.v, na.rm=TRUE)
            if (!anyNA(error.limits) &&
                error.limits[1] > error.limits[2]) {
                stop("error.limits[1] > error.limits[2]")
            }
            error.limits2 <- error.limits
            if (is.na(error.limits[1])) {
                error.limits2[1] <- min.E.v
                if (!is.na(error.limits[2]) &&
                    error.limits2[1] > error.limits[2]) {
                    warning("'error.limits[2]' too small, adjusting")
                    error.limits2[2] <- max.E.v
                }
            }
            if (is.na(error.limits2[2])) {
                error.limits2[2] <- max.E.v
                if (!is.na(error.limits[1]) &&
                    error.limits[1] > error.limits2[2]) {
                    warning("'error.limits[1]' too large, adjusting")
                    error.limits2[1] <- min.E.v
                }
            }
            E.v.norm <- E.v - error.limits2[1]
            E.v.norm <- E.v.norm / (error.limits2[2] - error.limits2[1])
            E.v.norm <- pmin(pmax(E.v.norm, 0), 1)
            not.na.idx <- which(!is.na(E.v.norm))
            color.mat <- matrix(rep(col2rgb(fillcolor.NA), each=n.nodes),
                                n.nodes, 3)
            color.mat[not.na.idx, ] <-
                colorRamp(c(fillcolor.small, fillcolor.large),
                          space = ramp.space2)(E.v.norm[not.na.idx])

            bordercolor[which.pruned] <-
                rgb(t(col2rgb(fillcolor.pruned)), maxColorValue=255)
            bordercolor[levelbest.loc] <-
                rgb(t(col2rgb(fillcolor.levelbest)), maxColorValue=255)
            bordercolor[c(L.v.loc, L.f.loc)] <-
                rgb(t(col2rgb(fillcolor.selected)), maxColorValue=255)
        } else {
            color.mat <- matrix(rep(col2rgb(fillcolor.normal), each=n.nodes),
                                n.nodes, 3)
            color.mat[which.pruned, ] <-
                matrix(rep(col2rgb(fillcolor.pruned), each=n.pruned),
                       n.pruned, 3)
            color.mat[levelbest.loc, ] <-
                matrix(rep(col2rgb(fillcolor.levelbest), each=D), D, 3)
            color.mat[c(L.v.loc, L.f.loc), ] <-
                matrix(rep(col2rgb(fillcolor.selected), each=2), 2, 3)
        }
        color.vec <- rgb(color.mat, maxColorValue=255)
        names(color.vec) <- nN
        ## Font color is black or white, depending on background color
        luma.vec <- as.numeric(color.mat %*% LUMA.CONVERSION / 255)
        fontcolor.vec <- rep.int("black", n.nodes)
        fontcolor.vec[luma.vec < 0.5] <- "white"
        names(fontcolor.vec) <- nN
        ## Check if "best nodes with a certain number of inputs"
        ## differ between "branching" and "no branching" solutions
        ## (when available)
        n.inputs <- vertex.data[["n.inputs"]]
        levelbest.n.inputs <- n.inputs[levelbest.loc]
        nobranch.loc <- which(vertex.data[["min.branches"]] == 1)
        nobranch.n.inputs <- n.inputs[nobranch.loc]
        dSeq <- seq_len(d - 1)
        loc1 <- levelbest.loc[match(dSeq, levelbest.n.inputs)]
        ## NA values possible in loc2 and diff.loc but it is not a problem
        loc2 <- nobranch.loc[match(dSeq, nobranch.n.inputs)]
        diff.loc <- loc2[loc1 != loc2]
        bordercolor[diff.loc] <-
            rgb(t(col2rgb(bordercolor.special.levelbest)),
                maxColorValue=255)
        nodeshape[diff.loc] <- "circle"
        ## Check if L.v and/or L.f sets differ between "branching" and
        ## "no branching" solutions (when available)
        L.v.nobranch <- x[["L.v.nobranch"]]
        if (!is.null(L.v.nobranch)) {
            L.f.nobranch <- x[["L.f.nobranch"]]
            L.v.nobranch.name <- paste0(L.v.nobranch, collapse=".")
            L.f.nobranch.name <- paste0(L.f.nobranch, collapse=".")
            L.v.nobranch.loc <- which(nN == L.v.nobranch.name)
            L.f.nobranch.loc <- which(nN == L.f.nobranch.name)
            L.v.diff <- L.v.nobranch.loc != L.v.loc
            L.f.diff <- L.f.nobranch.loc != L.f.loc
            if (L.v.diff || L.f.diff) {
                selected.rgb <- rgb(t(col2rgb(bordercolor.special.selected)),
                                    maxColorValue=255)
            }
            border.selected <-
                c(L.v.nobranch.loc, L.f.nobranch.loc)[c(L.v.diff, L.f.diff)]
            for (loc in border.selected) {
                bordercolor[loc] <- selected.rgb
                nodeshape[loc] <- "circle"
            }
            if (L.v.diff && L.f.diff && L.v.nobranch.loc == L.f.nobranch.loc) {
                node.labels[L.v.nobranch.loc] <- "v,f"
            } else {
                if (L.f.diff) {
                    node.labels[L.f.nobranch.loc] <- "f"
                }
                if (L.v.diff && L.v.nobranch.loc !=  L.f.loc) {
                    node.labels[L.v.nobranch.loc] <- "v"
                }
            }
        } else {
            border.selected <- numeric(0)
        }

        ## Number of occurrences of different border and fill colors
        n.border.selected <- length(border.selected)
        border.levelbest <-
            as.vector(na.omit(setdiff(diff.loc, border.selected)))
        n.border.levelbest <- length(border.levelbest)
        n.border.normal <- n.nodes - (n.border.selected + n.border.levelbest)
        if (color.by.error) {
            loc.selected <- setdiff(union(L.v.loc, L.f.loc), border.selected)
            n.selected <- length(loc.selected)
            n.levelbest <- length(setdiff(levelbest.loc,
                                          c(loc.selected, border.selected)))
            n.border.normal <- n.border.normal -
                (n.selected + n.levelbest + n.pruned)
        } else {
            loc.selected <- union(L.v.loc, L.f.loc)
            n.selected <- length(loc.selected)
            n.levelbest <- length(setdiff(levelbest.loc, loc.selected))
            n.normal <- n.nodes - (n.selected + n.levelbest + n.pruned)
        }

        ## Node labels
        node.labels <- rep.int("", n.nodes)
        names(node.labels) <- nN
        if (draw.node.labels) {
            which.onevar <- which(n.vars.in.node == 1)
            node.labels[which.onevar] <- nN[which.onevar]
        }
        if (draw.selected.labels) {
            node.labels[L.v.loc] <- "L.v"
            node.labels[L.f.loc] <- "L.f"
        }

        ## Edge labels
        eN <- graph::edgeNames(g)
        if (draw.edge.labels) {
            n.edges <- length(eN)
            eN.split <- strsplit(eN, "~", fixed=TRUE)
            eN.fromto <- unlist(eN.split)
            eN.from <- eN.fromto[seq(from=1, by=2, length.out=n.edges)]
            eN.to <- eN.fromto[seq(from=2, by=2, length.out=n.edges)]
            eN.dropped <- character(n.edges)
            for (k in seq_len(n.edges)) {
                vars.in.from <- name.to.vars.char(eN.from[k])
                vars.in.to <- name.to.vars.char(eN.to[k])
                eN.dropped[k] <- setdiff(vars.in.from, vars.in.to)
            }
            names(eN.dropped) <- eN
            edge.labels <- eN.dropped
        } else {
            edge.labels <- character(0)
            names(edge.labels) <- character(0)
        }
    }

### Plot
    ## Setting the attributes for Graphviz is rather tricky.  In
    ## addition to documentation, the approach below is based on trial
    ## and error. It seems to work but may contain some unnecessary
    ## steps.
    if (show[2]) {

        ## Handler for Graphviz attributes given by the user. Returns
        ## a character vector.  If the input x is a _named_ character,
        ## numeric or logical vector, the output is the character
        ## representation of x, with 'forbidden' attributes removed.
        ## The names are preserved.
        fix.extra.attrs <- function(x, forbidden=character(0)) {
            y <- character(0)
            names(y) <- character(0)
            x2 <- x
            if (is.numeric(x) || is.logical(x)) {
                x2 <- as.character(x)
                names(x2) <- names(x)
            }
            if (is.character(x2) && length(x2) > 0) {
                extra.names <- names(x2)
                if (!is.null(extra.names)) {
                    good.names <- unique(extra.names[!is.na(extra.names) &
                                                     nchar(extra.names) > 0])
                    good.names <- setdiff(good.names, forbidden)
                    y <- x2[good.names]
                }
            }
            y
        }

        ## From graph attributes, remove those that are controlled with
        ## separate parameters
        extra.graph.attrs <-
            fix.extra.attrs(graph.attributes,
                            c("nodesep", "ranksep", "rankdir"))
        ## From node and edge attributes, remove those that are set
        ## separately for each node and edge
        extra.node.attrs <-
            fix.extra.attrs(node.attributes,
                            c("fillcolor", "fontcolor", "label"))
        extra.edge.attrs <- fix.extra.attrs(edge.attributes, "label")

        ## Default "style" of nodes is "filled"
        if (is.na(extra.node.attrs["style"])) {
            extra.node.attrs["style"] <- "filled"
        }

        attrs <- list(node=as.list(extra.node.attrs),
                      edge=as.list(extra.edge.attrs),
                      graph=c(list(nodesep=nodesep, ranksep=ranksep,
                      rankdir=rankdir2), extra.graph.attrs))

        ## We create attribute lists where named vectors contain a
        ## copy of each global node and edge attribute for every node
        ## and edge.  The attribute lists are later combined with
        ## truly node- and edge-specific attributes and then used in
        ## buildNodeList() and buildEdgeList().
        n.node.attrs <- length(extra.node.attrs)
        node.attr.list <- vector("list", n.node.attrs)
        for (k in seq_len(n.node.attrs)) {
            temp.vec <- rep.int(extra.node.attrs[k], n.nodes)
            names(temp.vec) <- nN
            node.attr.list[[k]] <- temp.vec

        }
        names(node.attr.list) <- names(extra.node.attrs)
        n.edge.attrs <- length(extra.edge.attrs)
        edge.attr.list <- vector("list", n.edge.attrs)
        for (k in seq_len(n.edge.attrs)) {
            temp.vec <- rep.int(extra.edge.attrs[k], n.edges)
            names(temp.vec) <- eN
            edge.attr.list[[k]] <- temp.vec

        }
        names(edge.attr.list) <- names(extra.edge.attrs)
        ## It seems the attribute "penwidth" (for thickness of node
        ## borders) is ignored. It should work through the alternative
        ## (newer) layoutGraph() + renderGraph() interface, though.
        ## There the attribute is named "lwd".  However, renderGraph()
        ## appears not to support setting the margin of the plot (with
        ## "mai"), which is a new feature in the Ragraph method of
        ## plot().  Also, toFile() is convenient and works on Ragraphs
        ## (only available on the old interface).  It is confusing
        ## that there are (at least) two interfaces for laying out and
        ## plotting graphs, and each has its own advantages.
        nAttrs <- c(list(label=node.labels, fillcolor=color.vec,
                         fontcolor=fontcolor.vec, color=bordercolor,
                         shape=nodeshape),
                    node.attr.list)
        eAttrs <- c(list(label=edge.labels), edge.attr.list)
        nodes <- Rgraphviz::buildNodeList(g, nodeAttrs = nAttrs)
        edges <- Rgraphviz::buildEdgeList(g, edgeAttrs = eAttrs)
        ## Creates a 'Ragraph' object that can be plot()ted and
        ## written toFile()
        if (dev.cur() == 1) {
            vv <- Rgraphviz::agopen(name="foo", nodes=nodes, edges=edges,
                                    attrs=attrs, edgeMode="directed")
            ## Due to a bug in early Rgraphviz (version 1.35.2,
            ## others?), it's best to do the layout before opening a
            ## graphics device (if the default device size is OK). If
            ## one or several devices were already open, we postpone
            ## the layout procedure until the correct device has been
            ## selected, just before plotting. That way, the size of
            ## some graphical elements can be properly adjusted
            ## according to a possibly non-default device
            ## size. vv.done == FALSE indicates that the layout must
            ## be (re)computed just before plotting.
            vv.done <- dev.cur() == 1
        } else {
            vv.done <- FALSE
        }
    }

    ## Interactive plot
    if (make.plots) {
        plots.on.page <- 0
        max.plots.on.page <- prod(par("mfcol"))
        is.interactive <- dev.interactive()

        ## Call this before starting each new plot.  If the return
        ## value is TRUE, the caller should run devAskNewPage(TRUE)
        next.plot <- function() {
            plots.on.page <<- plots.on.page + 1
            if (is.interactive && plots.on.page == max.plots.on.page + 1) {
                if (dev.set) {
                    dev.set()
                    plots.on.page <<- 1
                    max.plots.on.page <<- prod(par("mfcol"))
                    is.interactive <<- dev.interactive()
                } else if (ask) {
                    return(TRUE)
                }
            }
            FALSE
        }

        onexit <- FALSE
        ## Plot 1: Error as a function of the number of variables
        ## Training error: Solid grey line with circles
        if (show[1]) {
            if (next.plot()) {
                originalAsk <- devAskNewPage(TRUE)
                on.exit(devAskNewPage(originalAsk))
                onexit <- TRUE
            }
            E.tr <- x[["E.tr"]]
            E.v <- x[["E.v"]]
            n.L.v <- length(x[["L.v"]])
            idx.L.v <- n.L.v + 1
            s.tr.temp <- x[["s.tr"]][idx.L.v]
            sd.y <- x[["sd.y"]]
            sd.zero.y <- isTRUE(all.equal(0, sd.y, check.attributes = FALSE))
            if (!sd.zero.y && standardize2 != orig.standardize) {
                var.y <- sd.y * sd.y
                if (standardize2) {
                    E.tr <- E.tr / var.y
                    E.v <- E.v / var.y
                    s.tr.temp <- s.tr.temp / var.y
                } else {
                    E.tr <- E.tr * var.y
                    E.v <- E.v * var.y
                    s.tr.temp <- s.tr.temp * var.y
                }
            }
            min.E.v <- E.v[idx.L.v]
            ## * Validation error: Longdash black line with bullets
            ## * Loc. of min. validation error and SD of training error
            ##   there: Solid black line
            ## * Threshold: Dotdash black line
            plotX <- matrix(NA_real_, max(D, 2), 4)
            plotY <- plotX
            tmpSeq <- seq_len(D)
            plotX[tmpSeq, 1:2] <- 0:d
            plotY[tmpSeq, 1] <- E.tr
            plotY[tmpSeq, 2] <- E.v
            plotX[1:2, 3] <- n.L.v
            plotY[1:2, 3] <- c(max(0, min.E.v - s.tr.temp),
                               min.E.v + s.tr.temp)
            plotX[1:2, 4] <- c(0, d)
            plotY[1:2, 4] <- min.E.v + s.tr.temp
            matplotArgs <- as.list(match.call(matplot,
                                        as.call(c(as.name("matplot"),
                                                  list(...)))))[-1L]
            if (is.list(plotArgs) && length(plotArgs) >= 1 &&
                is.list(plotArgs[[1]]) && length(plotArgs[[1]]) >= 1) {
                matplotArgs2 <- as.list(match.call(matplot,
                                             as.call(c(as.name("matplot"),
                                                       plotArgs[[1]]))))[-1L]
                matplotArgs[names(matplotArgs2)] <- matplotArgs2
            }
            ylim2 <- ylim
            if (is.null(ylim)) {
                xTmp <- max(1, floor(min(xlim)) + 1) :
                    min(D, ceiling(max(xlim)) + 1)
                E.range <- c(E.tr[xTmp], E.v[xTmp])
                ylim2 <- c(min(c(E.range, max(0, min.E.v - s.tr.temp))),
                           max(c(E.range, min.E.v + s.tr.temp)))
            }
            matplotArgs[["x"]] <- plotX
            matplotArgs[["y"]] <- plotY
            matplotArgs[["xlim"]] <- xlim
            matplotArgs[["ylim"]] <- ylim2
            argNames <- names(matplotArgs)
            if (!"type" %in% argNames) {
                matplotArgs[["type"]] <- c("b", "b", "l", "l")
            }
            if (!"pch" %in% argNames) {
                matplotArgs[["pch"]] <- c(1, 20)
            }
            if (!"col" %in% argNames) {
                matplotArgs[["col"]] <- c("grey", "black", "black", "black")
            }
            if (!"lty" %in% argNames) {
                matplotArgs[["lty"]] <- c(1, 5, 1, 4)
            }
            if (!"xlab" %in% argNames) {
                matplotArgs[["xlab"]] <- "Number of inputs"
            }
            if (!"ylab" %in% argNames) {
                matplotArgs[["ylab"]] <- "MSE"
            }
            do.call(matplot, matplotArgs, quote = TRUE)
        }

        ## Plot 2: Search graph
        if (show[2]) {
            if (next.plot()) {
                originalAsk <- devAskNewPage(TRUE)
                on.exit(devAskNewPage(originalAsk))
                onexit <- TRUE
            }
            if (!vv.done) {
                vv <- Rgraphviz::agopen(name="foo", nodes=nodes, edges=edges,
                                        attrs=attrs, edgeMode="directed")
            }
            graphplotArgs <- list(...)
            if (is.list(plotArgs) && length(plotArgs) >= 2 &&
                is.list(plotArgs[[2]]) && length(plotArgs[[2]]) >= 1) {
                argNames <- names(plotArgs[[2]])
                if (!is.null(argNames)) {
                    graphplotArgs2 <- plotArgs[[2]][!is.na(argNames) &
                                                    nchar(argNames) > 0]
                    graphplotArgs[names(graphplotArgs2)] <- graphplotArgs2
                }
            }
            graphplotArgs[["x"]] <- vv
            plotGraph <- function(fArgs) {
                if (isS4(vv)) {
                    ## Since plot is imported from the graphics package,
                    ## calling it the normal way would not invoke the S4
                    ## method for plot(Ragraph, ...) (defined in
                    ## Rgraphviz).  Plot needs to be imported, because we
                    ## define an S3 method for plot in this file and
                    ## register it in the NAMESPACE file. Also, the plot
                    ## function from graphics is actually used when
                    ## show[1] is TRUE.
                    do.call(selectMethod("plot", class(fArgs[["x"]])), fArgs,
                            quote = TRUE)
                } else {
                    do.call("plot", fArgs, quote = TRUE)
                }
            }

            ## Create plot 3: Color and shape key for plot 2.
            ## Drawing will be done later.
            if (show[3]) {
                colorkey.right    <- colorkey.space2 == "right"
                colorkey.bottom   <- colorkey.space2 == "bottom"
                colorkey.left     <- colorkey.space2 == "left"
                colorkey.top      <- colorkey.space2 == "top"
                colorkey.vertical <- colorkey.left || colorkey.right
                ## Width of a grob or a list of grobs in inches
                widthInches <- function(x) {
                    if (inherits(x, "grob")) {
                        convertWidth(grobWidth(x), unitTo="inches",
                                     valueOnly=TRUE)
                    } else if (is.vector(x)) {
                        vapply(x, widthInches, numeric(1))
                    } else {
                        NA_real_
                    }
                }
                ## Height of a grob or a list of grobs in inches
                heightInches <- function(x) {
                    if (inherits(x, "grob")) {
                        convertHeight(grobHeight(x), unitTo="inches",
                                      valueOnly=TRUE)
                    } else if (is.vector(x)) {
                        vapply(x, heightInches, numeric(1))
                    } else {
                        NA_real_
                    }
                }
                ## Finds size of gap between the "color" and "labels"
                ## parts of a colorkey object.  Assumes things about
                ## the structure of the object but falls back to a
                ## default value if necessary.
                colorkeyGap <- function(x, vertical,
                                        default = unit(0.6, "lines")) {
                    if (inherits(x, "frame") && inherits(x, "grob")) {
                        if (vertical) {
                            measures <- x[["framevp"]][["layout"]][["widths"]]
                        } else {
                            measures <- x[["framevp"]][["layout"]][["heights"]]
                        }
                        if (inherits(measures, "unit") &&
                            length(measures) == 3) {
                            measures[2]
                        } else {
                            default
                        }
                    } else {
                        default
                    }
                }
                ## Accepts a text grob x possibly representing many
                ## pieces of text, e.g. length(x) > 1.
                ## Alternatively accepts a gTree x, in which case the
                ## function recursively searches for text grobs.
                ## Returns a list with one grob for each piece of text.
                harvestText <- function(x) {
                    if (inherits(x, "gTree")) {
                        cNames <- childNames(x)
                        res <- list()
                        for (name in cNames) {
                            res <- c(res, harvestText(getGrob(x, name)))
                        }
                        res
                    } else if (inherits(x, "text") && inherits(x, "grob")) {
                        X <- x[["x"]]
                        Y <- x[["y"]]
                        n <- max(length(X), length(Y))
                        X <- rep(X, length.out = n)
                        Y <- rep(Y, length.out = n)
                        label <- rep_len(x[["label"]], n)
                        rot <- rep_len(x[["rot"]], n)
                        res <- vector(mode = "list", length = n)
                        name <- paste(x[["name"]], 1:n, sep=".")
                        for (k in seq_len(n)) {
                            res[[k]] <- editGrob(x, x = X[k], y = Y[k],
                                                 label = label[k], rot = rot[k],
                                                 name = name[k])
                        }
                        res
                    }
                }
                ## Create keys explaining the colors and shapes used in plot 2
                createKeys <- function() {
                    if (colorkey.right) {
                        titleX <- 0
                        titleJust <- "left"
                    } else if (colorkey.bottom) {
                        titleX <- 0.5
                        titleJust <- "centre"
                    } else if (colorkey.left) {
                        titleX <- 1
                        titleJust <- "right"
                    } else {
                        titleX <- 0.5
                        titleJust <- "centre"
                    }
                    if (is.list(plotArgs) && length(plotArgs) >= 3 &&
                        is.list(plotArgs[[3]])) {
                        ck.key <- plotArgs[[3]]
                    } else {
                        ck.key <- list()
                    }
                    ck.key[["space"]] <- colorkey.space2
                    border.ck.key <- ck.key
                    have.label.list <-
                        ("labels" %in% names(ck.key) &&
                         is.list(ck.key[["labels"]]))

                    if (color.by.error) {
                        key.colors <-
                            colorRampPalette(c(fillcolor.small,
                                               fillcolor.large),
                                             space = ramp.space2)(ramp.size)
                        col.at <- seq(from = error.limits2[1],
                                      to = error.limits2[2],
                                      length.out = ramp.size)
                        col.at <- c(col.at,
                                    col.at[ramp.size] + col.at[2] - col.at[1])
                        if (min.E.v < col.at[1]) {
                            col.at[1] <- min.E.v
                        }
                        if (max.E.v > col.at[ramp.size + 1]) {
                            col.at[ramp.size + 1] <- max.E.v
                        }
                        n.fill.styles <- 4 # value affects layout proportions
                        border.used <- c(n.border.normal > 0, n.pruned > 0,
                                         n.levelbest > 0, n.selected > 0,
                                         n.border.levelbest > 0,
                                         n.border.selected > 0)
                    } else {
                        fill.used <- c(n.pruned > 0,
                                       n.levelbest > 0, n.selected > 0,
                                       n.normal > 0)
                        key.colors <- rbind(t(col2rgb(fillcolor.pruned)),
                                            t(col2rgb(fillcolor.levelbest)),
                                            t(col2rgb(fillcolor.selected)),
                                            t(col2rgb(fillcolor.normal)))
                        key.colors <- rgb(key.colors[fill.used, , drop=FALSE],
                                          maxColorValue=255)
                        n.fill.styles <- length(key.colors)
                        col.at <- seq_len(n.fill.styles + 1)
                        key.labels <-
                            category.labels[c("pruned", "levelbest",
                                              "selected", "normal")[fill.used]]
                        if (colorkey.vertical) {
                            key.colors <- rev(key.colors)
                            key.labels <- rev(key.labels)
                        }
                        if (have.label.list) {
                            ck.key[["labels"]][["labels"]] <-
                                as.character(key.labels)
                            ck.key[["labels"]][["at"]] <-
                                seq_len(n.fill.styles) + 0.5
                        } else {
                            ck.key[["labels"]] <-
                                list(labels = as.character(key.labels),
                                     at = seq_len(n.fill.styles) + 0.5)
                        }
                        border.used <-
                            c(FALSE, FALSE, FALSE, n.border.levelbest > 0,
                              n.border.selected > 0, n.border.normal > 0)
                    }

                    ck.key[["at"]] <- col.at
                    ck.key[["col"]] <- key.colors
                    ck <- draw.colorkey(key = ck.key)
                    width.fill <- widthInches(ck)
                    height.fill <- heightInches(ck)
                    ppMar <- convertWidth(colorkeyGap(ck, colorkey.vertical),
                                          unitTo = "inches",
                                          valueOnly = TRUE)
                    if (colorkey.vertical) {
                        if (colorkey.right) {
                            ppW <- unit.c(unit(width.fill, "inches"),
                                          unit(1, "null"))
                            contentCol <- 1
                        } else {
                            ppW <- unit.c(unit(1, "null"),
                                          unit(width.fill, "inches"))
                            contentCol <- 2
                        }
                        ppLayout <- grid.layout(nrow = 1, ncol = length(ppW),
                                                widths = ppW)
                        childPorts <-
                            vpList(viewport(layout.pos.row = 1,
                                            layout.pos.col = contentCol,
                                            name = "fill.content"))
                    } else {
                        if (colorkey.top) {
                            ppH <- unit.c(unit(1, "null"),
                                          unit(height.fill, "inches"))
                            contentRow <- 2
                        } else {
                            ppH <- unit.c(unit(height.fill, "inches"),
                                          unit(1, "null"))
                            contentRow <- 1
                        }
                        ppLayout <- grid.layout(nrow = length(ppH), ncol = 1,
                                                heights = ppH)
                        childPorts <-
                            vpList(viewport(layout.pos.row = contentRow,
                                            layout.pos.col = 1,
                                            name = "fill.content"))
                    }
                    parentPort <- viewport(layout = ppLayout, name = "fill")
                    fillPorts <- vpTree(parent = parentPort,
                                        children = childPorts)
                    fillGrobs <-
                        gList(editGrob(ck, vp = vpPath("fill", "fill.content")))
                    ck.tree <- gTree(children = fillGrobs, name = "fillTree",
                                     childrenvp = fillPorts)

                    fillTextGrobs <- harvestText(ck)
                    width.fill.text <- max(widthInches(fillTextGrobs))
                    height.fill.text <- max(heightInches(fillTextGrobs))
                    if (colorkey.vertical) {
                        maxSize <- height.fill.text
                    } else {
                        maxSize <- width.fill.text
                    }
                    n.border.styles <- sum(border.used)
                    do.border.ck <- n.border.styles > 1
                    n.shapes <- 1 + ("circle" %in% nodeshape)
                    do.shape.key <- n.shapes > 1
                    key.titles <-
                        list(fill = textGrob(gettext("Fill", domain="R-sisal"),
                             x = titleX, just = titleJust,
                             gp = colorkey.title.gp))
                    if (do.border.ck) {
                        key.titles[["border"]] <-
                            textGrob(gettext("Border", domain="R-sisal"),
                                     x = titleX, just = titleJust,
                                     gp = colorkey.title.gp)
                    }
                    if (do.shape.key) {
                        key.titles[["shape"]] <-
                            textGrob(gettext("Shape", domain="R-sisal"),
                                     x = titleX, just = titleJust,
                                     gp = colorkey.title.gp)
                    }
                    marginNum <- 0.5
                    margin <- unit(marginNum, "null")
                    titleH <- heightInches(key.titles)
                    titleW <- widthInches(key.titles)
                    titleD <- vapply(key.titles, function (x) {
                        convertHeight(grobDescent(x), unitTo = "inches",
                                      valueOnly = TRUE)
                    }, numeric(1))
                    titleHeight <- titleH + 2 * titleD +
                        convertHeight(unit(1, "lines"), unitTo = "inches",
                                      valueOnly = TRUE)
                    if (colorkey.vertical) {
                        colorkey.size <- max(width.fill, titleW)
                        widths <- 1
                        title.row <- c(fill = 1)
                        title.col <- c(fill = 1)
                        content.row <- c(fill = 2)
                        content.col <- c(fill = 1)
                        heights <- unit.c(unit(titleHeight["fill"], "inches"),
                                          unit(n.fill.styles, "null"))
                        if (do.border.ck) {
                            heights <-
                                unit.c(heights, margin,
                                       unit(titleHeight["border"], "inches"),
                                       unit(n.border.styles, "null"))
                            title.row <- c(title.row, border = 4)
                            title.col <- c(title.col, border = 1)
                            content.row <- c(content.row, border = 5)
                            content.col <- c(content.col, border = 1)
                        }
                        if (do.shape.key) {
                            heights <-
                                unit.c(unit(titleHeight["shape"], "inches"),
                                       unit(n.shapes, "null"), margin, heights)
                            title.row <- c(title.row + 3, shape = 1)
                            title.col <- c(title.col, shape = 1)
                            content.row <- c(content.row + 3, shape = 2)
                            content.col <- c(content.col, shape = 1)
                        }
                    } else { # horizontal
                        maxContentHeight <- height.fill
                        if (colorkey.top) {
                            heights <- unit.c(unit(max(titleHeight), "inches"),
                                              unit(1, "null"))
                            title.row <- c(fill = 1)
                            content.row <- c(fill = 2)
                        } else {
                            heights <- unit.c(unit(1, "null"),
                                              unit(max(titleHeight), "inches"))
                            title.row <- c(fill = 2)
                            content.row <- c(fill = 1)
                        }
                        title.col <- c(fill = 1)
                        content.col <- c(fill = 1)
                        widths <- unit.c(unit(n.fill.styles, "null"))
                        if (do.border.ck) {
                            widths <- unit.c(widths, margin,
                                             unit(n.border.styles, "null"))
                            title.row <- c(title.row,
                                           border = as.vector(title.row[1]))
                            title.col <- c(title.col, border = 3)
                            content.row <- c(content.row,
                                             border = as.vector(content.row[1]))
                            content.col <- c(content.col, border = 3)
                        }
                        if (do.shape.key) {
                            widths <- unit.c(unit(n.shapes, "null"),
                                             margin, widths)
                            title.row <- c(title.row,
                                           shape = as.vector(title.row[1]))
                            title.col <- c(title.col + 2, shape = 1)
                            content.row <- c(content.row,
                                             shape = as.vector(content.row[1]))
                            content.col <- c(content.col + 2, shape = 1)
                        }
                    }
                    nRow <- length(heights)
                    nCol <- length(widths)
                    width.border <- 0
                    if (do.border.ck) {
                        border.key.colors <-
                            rbind(t(col2rgb(fillcolor.pruned)),
                                  t(col2rgb(fillcolor.levelbest)),
                                  t(col2rgb(fillcolor.selected)),
                                  t(col2rgb(bordercolor.special.levelbest)),
                                  t(col2rgb(bordercolor.special.selected)),
                                  t(col2rgb(bordercolor.normal)))
                        border.key.colors <-
                            rgb(border.key.colors[border.used, , drop=FALSE],
                                maxColorValue=255)
                        border.key.colors <- rep(border.key.colors, each=2)
                        border.key.labels <-
                            category.labels[c("pruned", "levelbest",
                                              "selected", "special.levelbest",
                                              "special.selected", "normal")]
                        border.key.labels <- border.key.labels[border.used]
                        if (colorkey.vertical) {
                            border.key.colors <- rev(border.key.colors)
                            border.key.labels <- rev(border.key.labels)
                        }
                        if (have.label.list) {
                            border.ck.key[["labels"]][["labels"]] <-
                                as.character(border.key.labels)
                            border.ck.key[["labels"]][["at"]] <-
                                seq_len(n.border.styles) + 0.5
                        } else {
                            border.ck.key[["labels"]] <-
                                list(labels = as.character(border.key.labels),
                                     at = seq_len(n.border.styles) + 0.5)
                        }
                        border.ck.key[["at"]] <- seq_len(n.border.styles + 1)
                        border.ck.key[["col"]] <- border.key.colors
                        border.ck <- draw.colorkey(key = border.ck.key)
                        width.border <- widthInches(border.ck)
                        height.border <- heightInches(border.ck)
                        if (colorkey.vertical) {
                            if (colorkey.right) {
                                ppW <- unit.c(unit(width.border, "inches"),
                                              unit(1, "null"))
                                contentCol <- 1
                            } else {
                                ppW <- unit.c(unit(1, "null"),
                                              unit(width.border, "inches"))
                                contentCol <- 2
                            }
                            ppLayout <- grid.layout(nrow = 1,
                                                    ncol = length(ppW),
                                                    widths = ppW)
                            childPorts <-
                                vpList(viewport(layout.pos.row = 1,
                                                layout.pos.col = contentCol,
                                                name = "border.content"))
                        } else {
                            if (colorkey.top) {
                                ppH <- unit.c(unit(1, "null"),
                                              unit(height.border, "inches"))
                                contentRow <- 2
                            } else {
                                ppH <- unit.c(unit(height.border, "inches"),
                                              unit(1, "null"))
                                contentRow <- 1
                            }
                            ppLayout <- grid.layout(nrow = length(ppH),
                                                    ncol = 1,
                                                    heights = ppH)
                            childPorts <-
                                vpList(viewport(layout.pos.row = contentRow,
                                                layout.pos.col = 1,
                                                name = "border.content"))
                        }
                        parentPort <- viewport(layout = ppLayout,
                                               name = "border")
                        borderPorts <- vpTree(parent = parentPort,
                                            children = childPorts)
                        borderGrobs <-
                            gList(editGrob(border.ck,
                                           vp = vpPath("border",
                                           "border.content")))
                        border.ck.tree <- gTree(children = borderGrobs,
                                                name = "borderTree",
                                                childrenvp = borderPorts)

                        borderTextGrobs <- harvestText(border.ck)
                        height.border.text <-
                            max(heightInches(borderTextGrobs))
                        if (colorkey.vertical) {
                            maxSize <- max(maxSize,
                                           height.border.text)
                            colorkey.size <- max(colorkey.size, width.border)
                        } else { # horizontal
                            maxSize <- max(maxSize,
                                           widthInches(borderTextGrobs))
                            maxContentHeight <-
                                max(maxContentHeight, height.border)
                        }
                    }
                    if (do.shape.key) {
                        textGp <- fillTextGrobs[[1]][["gp"]]
                        labelNormal <- category.labels["shape.normal"]
                        labelSpecial <- category.labels["shape.highlighted"]
                        if (colorkey.vertical) {
                            textY <- unit(0.5, "npc")
                            if (colorkey.right) {
                                textJust <- c("left", "centre")
                                textX <- unit(0, "npc")
                            } else {
                                textJust <- c("right", "centre")
                                textX <- unit(1, "npc")
                            }
                        } else {
                            textX <- unit(0.5, "npc")
                            if (colorkey.top) {
                                textJust <- c("centre", "bottom")
                                textY <- unit(0, "npc")
                            } else {
                                textJust <- c("centre", "top")
                                textY <- unit(1, "npc")
                            }
                        }
                        normalTextGrob <- textGrob(labelNormal,
                                                   x = textX, y = textY,
                                                   just = textJust,
                                                   gp = textGp,
                                                   name = "normalText")
                        specialTextGrob <- editGrob(normalTextGrob,
                                                    label = labelSpecial,
                                                    name = "specialText")
                        width.shape.text <- max(widthInches(normalTextGrob),
                                                widthInches(specialTextGrob))
                        height.shape.text <- max(heightInches(normalTextGrob),
                                                 heightInches(specialTextGrob))
                        if (colorkey.vertical) {
                            maxSize <- max(maxSize, height.shape.text)
                            ppSymbol <- width.fill - ppMar - width.fill.text
                            if (colorkey.right) {
                                ppW <-
                                    unit.c(unit(c(ppSymbol, ppMar), "inches"),
                                           unit(1, "null"))
                                symbolCol <- 1
                                textCol <- 3
                            } else {
                                ppW <-
                                    unit.c(unit(1, "null"),
                                           unit(c(ppMar, ppSymbol), "inches"))
                                symbolCol <- 3
                                textCol <- 1
                            }
                            ppLayout <- grid.layout(nrow = 2, widths = ppW,
                                                    ncol = length(ppW))
                            childPorts <-
                                vpList(viewport(layout.pos.row = 2,
                                                layout.pos.col = symbolCol,
                                                name = "normal.symbol"),
                                       viewport(layout.pos.row = 2,
                                                layout.pos.col = textCol,
                                                name = "normal.text"),
                                       viewport(layout.pos.row = 1,
                                                layout.pos.col = symbolCol,
                                                name = "special.symbol"),
                                       viewport(layout.pos.row = 1,
                                                layout.pos.col = textCol,
                                                name = "special.text"))
                            colorkey.size <-
                                max(colorkey.size,
                                    ppSymbol + ppMar + width.shape.text)
                        } else { # horizontal
                            maxSize <- max(maxSize, width.shape.text)
                            ppSymbol <- height.fill - ppMar - height.fill.text
                            if (colorkey.top) {
                                ppH <-
                                    unit.c(unit(1, "null"),
                                           unit(c(ppMar, ppSymbol), "inches"))
                                symbolRow <- 3
                                textRow <- 1
                            } else {
                                ppH <-
                                    unit.c(unit(c(ppSymbol, ppMar), "inches"),
                                           unit(1, "null"))
                                symbolRow <- 1
                                textRow <- 3
                            }
                            ppLayout <- grid.layout(nrow = length(ppH),
                                                    ncol = 2, heights = ppH)
                            childPorts <-
                                vpList(viewport(layout.pos.row = symbolRow,
                                                layout.pos.col = 2,
                                                name = "normal.symbol"),
                                       viewport(layout.pos.row = textRow,
                                                layout.pos.col = 2,
                                                name = "normal.text"),
                                       viewport(layout.pos.row = symbolRow,
                                                layout.pos.col = 1,
                                                name = "special.symbol"),
                                       viewport(layout.pos.row = textRow,
                                                layout.pos.col = 1,
                                                name = "special.text"))
                            maxContentHeight <-
                                max(maxContentHeight,
                                    ppSymbol + ppMar + height.shape.text)
                        }
                        parentPort <- viewport(layout = ppLayout,
                                               name = "shape")
                        shapePorts <- vpTree(parent = parentPort,
                                             children = childPorts)
                        normalTextGrob <-
                            editGrob(normalTextGrob,
                                     vp = vpPath("shape", "normal.text"))
                        specialTextGrob <-
                            editGrob(specialTextGrob,
                                     vp = vpPath("shape", "special.text"))
                        normalSymbolGrob <-
                            rectGrob(width = 1, height = 1,
                                     default.units = "snpc",
                                     vp = vpPath("shape", "normal.symbol"),
                                     name = "normalSymbol")
                        specialSymbolGrob <-
                            circleGrob(vp = vpPath("shape", "special.symbol"),
                                       name = "specialSymbol")
                        shapeGrobs <-
                            gList(normalSymbolGrob, normalTextGrob,
                                  specialSymbolGrob, specialTextGrob)
                        shape.key <- gTree(children = shapeGrobs,
                                           childrenvp = shapePorts,
                                           name = "shapeTree")
                    }
                    if (!colorkey.vertical) {
                        colorkey.size <- max(titleHeight) + maxContentHeight
                    }
                    drawKeys <- function() {
                        totalNulls <- n.fill.styles +
                            do.border.ck * (n.border.styles + marginNum) +
                                do.shape.key * (n.shapes + marginNum)
                        if (colorkey.vertical) {
                            totalInches <- titleHeight["fill"]
                            if (do.border.ck) {
                                totalInches <-
                                    totalInches + titleHeight["border"]
                            }
                            if (do.shape.key) {
                                totalInches <-
                                    totalInches + titleHeight["shape"]
                            }
                            keys.port.width <- unit(colorkey.size, "inches")
                            ## magic constant 2.5
                            keys.port.height <-
                                min(unit(1, "npc"),
                                    unit(totalInches + 3 * totalNulls * maxSize,
                                         "inches"))
                        } else  {
                            ## magic constant 1.25
                            keys.port.width <-
                                min(unit(1, "npc"),
                                    unit(max(titleW["fill"] / n.fill.styles,
                                             if (do.border.ck) {
                                                 titleW["border"] /
                                                     n.border.styles
                                             },
                                             if (do.shape.key) {
                                                 titleW["shape"] /
                                                     n.shapes
                                             },
                                             1.25 * maxSize) * totalNulls,
                                         "inches"))
                            keys.port.height <- unit(colorkey.size, "inches")
                        }
                        pushViewport(viewport(layout = grid.layout(nrow = nRow,
                                              ncol = nCol, heights = heights,
                                              widths = widths),
                                              width = keys.port.width,
                                              height = keys.port.height,
                                              name = "keys"))
                        drawOneKey <- function(x, name) {
                            tmp <-
                                viewport(layout.pos.row = content.row[name],
                                         layout.pos.col = content.col[name],
                                         name = paste("key", name, sep="."))
                            pushViewport(tmp)
                            grid.draw(x)
                            popViewport(1)
                            tmp <-
                                viewport(layout.pos.row = title.row[name],
                                         layout.pos.col = title.col[name],
                                         name = paste("title", name, sep="."))
                            pushViewport(tmp)
                            grid.draw(key.titles[[name]])
                            popViewport(1)
                        }
                        ## Draw border key
                        if (do.border.ck) {
                            drawOneKey(border.ck.tree, "border")
                        }
                        ## Draw shape key
                        if (do.shape.key) {
                            drawOneKey(shape.key, "shape")
                        }
                        ## Draw fill key
                        drawOneKey(ck.tree, "fill")
                        popViewport(1)
                    }
                    list(size = colorkey.size, draw = drawKeys)
                }
            }

            op <- par(no.readonly = TRUE)
            nochange <- c("ask", "fin", "new", "mar", "mfcol", "mfg",
                          "mfrow", "oma", "omd", "pin")
            op <- op[!(names(op) %in% nochange)]
            on.exit(par(op), add = TRUE)
            if (show[3] && integrate.colorkey2) {
                grid.newpage()
                plot(1:5, type="n", axes=FALSE, ann=FALSE)
                ## Create layout (mixed base and grid graphics...)
                mai <- graphplotArgs[["mai"]]
                if (is.null(mai)) {
                    mai <- rep.int(0.1, 4)
                }
                keyList <- createKeys()
                ck.size <- keyList[["size"]]
                bottomMargin <- unit(mai[1], "inches")
                leftMargin <- unit(mai[2], "inches")
                topMargin <- unit(mai[3], "inches")
                rightMargin <- unit(mai[4], "inches")
                nullWidths <- unit.c(leftMargin, unit(1, "null"), rightMargin)
                nullHeights <- unit.c(topMargin, unit(1, "null"), bottomMargin)
                if (colorkey.right) {
                    keys.pos.row <- 2
                    keys.pos.col <- 3
                    pushViewport(viewport(layout = grid.layout(nrow = 3,
                                          ncol = 4,
                                          widths = unit.c(leftMargin,
                                          unit(1, "null"),
                                          unit(ck.size, "inches"),
                                          rightMargin),
                                          heights = nullHeights),
                                          name = "top"))
                    pushViewport(viewport(layout.pos.row = 2,
                                          layout.pos.col = 2))
                    mai[4] <- mai[4] + ck.size + colorkey.gap
                } else if (colorkey.bottom) {
                    keys.pos.row <- 3
                    keys.pos.col <- 2
                    pushViewport(viewport(layout = grid.layout(nrow = 4,
                                          ncol = 3,
                                          widths = nullWidths,
                                          heights = unit.c(topMargin,
                                          unit(1, "null"),
                                          unit(ck.size, "inches"),
                                          bottomMargin)),
                                          name = "top"))
                    pushViewport(viewport(layout.pos.row = 2,
                                          layout.pos.col = 2))
                    mai[1] <- mai[1] + ck.size + colorkey.gap
                } else if (colorkey.left) {
                    keys.pos.row <- 2
                    keys.pos.col <- 2
                    pushViewport(viewport(layout = grid.layout(nrow = 3,
                                          ncol = 4,
                                          widths = unit.c(leftMargin,
                                          unit(ck.size, "inches"),
                                          unit(1, "null"),
                                          rightMargin),
                                          heights = nullHeights),
                                          name = "top"))
                    pushViewport(viewport(layout.pos.row = 2,
                                          layout.pos.col = 3))
                    mai[2] <- mai[2] + ck.size + colorkey.gap
                } else {
                    keys.pos.row <- 2
                    keys.pos.col <- 2
                    pushViewport(viewport(layout = grid.layout(nrow = 4,
                                          ncol = 3,
                                          widths = nullWidths,
                                          heights = unit.c(topMargin,
                                          unit(ck.size, "inches"),
                                          unit(1, "null"), bottomMargin)),
                                          name = "top"))
                    pushViewport(viewport(layout.pos.row = 3,
                                          layout.pos.col = 2))
                    mai[3] <- mai[3] + ck.size + colorkey.gap
                }
                graphplotArgs[["mai"]] <- mai
                par(new=TRUE)
                plotGraph(graphplotArgs)
                popViewport(1)
            } else {
                plotGraph(graphplotArgs)
            }
            par(op)
            if (onexit) {
                on.exit(devAskNewPage(originalAsk))
            } else {
                on.exit()
            }
        }

        ## Plot 3: Color and shape key for plot 2.
        if (show[2] && show[3]) {
            if (!integrate.colorkey2) {
                if (next.plot()) {
                    originalAsk <- devAskNewPage(TRUE)
                    on.exit(devAskNewPage(originalAsk))
                }
                grid.newpage()
                keyList <- createKeys()
                keys.pos.row <- NULL
                keys.pos.col <- NULL
            }
            pushViewport(viewport(layout.pos.row = keys.pos.row,
                                  layout.pos.col = keys.pos.col,
                                  name = "keys.top"))
            keyList[["draw"]]()
            popViewport(1)
        }
    }
    invisible(res)
}
