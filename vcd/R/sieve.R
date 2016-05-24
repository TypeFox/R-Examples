###########################################################
## sieveplot

sieve <- function(x, ...)
    UseMethod("sieve")

sieve.formula <-
    function(formula, data = NULL, ..., main = NULL, sub = NULL, subset = NULL)
{
    if (is.logical(main) && main)
        main <- deparse(substitute(data))
    else if (is.logical(sub) && sub)
        sub <- deparse(substitute(data))

    m <- match.call(expand.dots = FALSE)
    edata <- eval(m$data, parent.frame())

    fstr <- strsplit(paste(deparse(formula), collapse = ""), "~")
    vars <- strsplit(strsplit(gsub(" ", "", fstr[[1]][2]), "\\|")[[1]], "\\+")
    varnames <- vars[[1]]
    condnames <- if (length(vars) > 1) vars[[2]] else NULL

    if (inherits(edata, "ftable") || inherits(edata, "table") || length(dim(edata)) > 2) {
        condind <- NULL
        dat <- as.table(data)
        if(all(varnames != ".")) {
            ind <- match(varnames, names(dimnames(dat)))
            if (any(is.na(ind)))
                stop(paste("Can't find", paste(varnames[is.na(ind)], collapse=" / "), "in", deparse(substitute(data))))

            if (!is.null(condnames)) {
                condind <- match(condnames, names(dimnames(dat)))
                if (any(is.na(condind)))
                    stop(paste("Can't find", paste(condnames[is.na(condind)], collapse=" / "), "in", deparse(substitute(data))))
                ind <- c(condind, ind)
            }
            dat <- margin.table(dat, ind)
        }
        sieve.default(dat, main = main, sub = sub,
                      condvars = if (is.null(condind)) NULL else match(condnames, names(dimnames(dat))), ...)
    } else {
        tab <- if ("Freq" %in% colnames(data))
            xtabs(formula(paste("Freq~", paste(c(condnames, varnames), collapse = "+"))),
                  data = data, subset = subset)
        else
            xtabs(formula(paste("~", paste(c(condnames, varnames), collapse = "+"))),
                  data = data, subset = subset)

        sieve.default(tab, main = main, sub = sub, ...)
    }
}

sieve.default <- function(x, condvars = NULL, gp = NULL,
                          shade = NULL, legend = FALSE,
                          split_vertical = NULL, direction = NULL,
                          spacing = NULL, spacing_args = list(),
                          sievetype = c("observed","expected"),
                          gp_tile = gpar(),
                          scale = 1,
                          main = NULL, sub = NULL, ...) {
    if (is.logical(main) && main)
        main <- deparse(substitute(x))
    else if (is.logical(sub) && sub)
        sub <- deparse(substitute(x))
    sievetype = match.arg(sievetype)
    if (is.logical(shade) && shade && is.null(gp))
        gp <- if (sievetype == "observed")
#            shading_sieve(interpolate = 0, lty = c("longdash", "solid"))
            shading_sieve(interpolate = 0, lty = c("solid", "longdash"))
	else
            shading_sieve(interpolate = 0, line_col = "darkgray", eps = Inf, lty = "dotted")

    if (is.structable(x)) {
        if (is.null(direction) && is.null(split_vertical))
            split_vertical <- attr(x, "split_vertical")
        x <- as.table(x)
    }
    if (is.null(split_vertical))
        split_vertical <- FALSE

    dl <- length(dim(x))

    ## splitting argument
    if (!is.null(direction))
        split_vertical <- direction == "v"
    if (length(split_vertical) == 1)
        split_vertical <- rep(c(split_vertical, !split_vertical), length.out = dl)
    if (length(split_vertical) < dl)
        split_vertical <- rep(split_vertical, length.out = dl)

    ## condvars
    if (!is.null(condvars)) {
        if (is.character(condvars))
            condvars <- match(condvars, names(dimnames(x)))
        x <- aperm(x, c(condvars, seq(dl)[-condvars]))
        if (is.null(spacing))
            spacing <- spacing_conditional
    }

    ## spacing argument
    if (is.null(spacing))
        spacing <- if (dl < 3) spacing_equal(sp = 0) else spacing_increase

    strucplot(x,
              condvars = if (is.null(condvars)) NULL else length(condvars),
              core = struc_sieve(sievetype = sievetype, gp_tile = gp_tile,
                                 scale = scale),
              split_vertical = split_vertical,
              spacing = spacing,
              spacing_args = spacing_args,
              main = main,
              sub = sub,
              shade = shade,
              legend = legend,
              gp = gp,
              ...)
}

## old version (not performant enough)
##
## struc_sieve <- function(sievetype = c("observed", "expected")) {
##   sievetype = match.arg(sievetype)
##   function(residuals, observed, expected, spacing, gp, split_vertical, prefix = "") {
##     dn <- dimnames(expected)
##     dnn <- names(dn)
##     dx <- dim(expected)
##     dl <- length(dx)
##     n <- sum(expected)

##     ## split workhorse
##     split <- function(x, i, name, row, col, rowmargin, colmargin) {
##       cotab <- co_table(x, 1)
##       margin <- sapply(cotab, sum)
##       v <- split_vertical[i]
##       d <- dx[i]

##       ## compute total cols/rows and build split layout
##       dist <- unit.c(unit(margin, "null"), spacing[[i]])
##       idx <- matrix(1:(2 * d), nrow = 2, byrow = TRUE)[-2 * d]
##       layout <- if (v)
##         grid.layout(ncol = 2 * d - 1, widths = dist[idx])
##       else
##         grid.layout(nrow = 2 * d - 1, heights = dist[idx])
##       vproot <- viewport(layout.pos.col = col, layout.pos.row = row,
##                          layout = layout, name = remove_trailing_comma(name))

##       ## next level: either create further splits, or final viewports
##       name <- paste(name, dnn[i], "=", dn[[i]], ",", sep = "")
##       row <- col <- rep.int(1, d)
##       if (v) col <- 2 * 1:d - 1 else row <- 2 * 1:d - 1
##       proptab <- function(x) x / max(sum(x), 1)
##       f <- if (i < dl) {
##         if (v)
##           function(m) split(cotab[[m]], i + 1, name[m], row[m], col[m],
##                             colmargin = colmargin * proptab(margin)[m],
##                             rowmargin = rowmargin)
##         else
##           function(m) split(cotab[[m]], i + 1, name[m], row[m], col[m],
##                             colmargin = colmargin,
##                             rowmargin = rowmargin * proptab(margin)[m])
##       } else {
##         if (v)
##           function(m) viewport(layout.pos.col = col[m], layout.pos.row = row[m],
##                                name = remove_trailing_comma(name[m]),
##                                yscale = c(0, rowmargin),
##                                xscale = c(0, colmargin * proptab(margin)[m]))
##         else
##           function(m) viewport(layout.pos.col = col[m], layout.pos.row = row[m],
##                                name = remove_trailing_comma(name[m]),
##                                yscale = c(0, rowmargin * proptab(margin)[m]),
##                                xscale = c(0, colmargin))
##       }
##       vpleaves <- structure(lapply(1:d, f), class = c("vpList", "viewport"))

##       vpTree(vproot, vpleaves)
##     }

##     ## start splitting on top, creates viewport-tree
##     pushViewport(split(expected + .Machine$double.eps,
##                        i = 1, name = paste(prefix, "cell:", sep = ""), row = 1, col = 1,
##                        rowmargin = n, colmargin = n))

##     ## draw rectangles
##     mnames <- apply(expand.grid(dn), 1,
##                     function(i) paste(dnn, i, collapse=",", sep = "=")
##                     )

##     for (i in seq_along(mnames)) {
##       seekViewport(paste(prefix, "cell:", mnames[i], sep = ""))
##       vp <- current.viewport()
##       gpobj <- structure(lapply(gp, function(x) x[i]), class = "gpar")

##       div <- if (sievetype == "observed") observed[i] else expected[i]
##       if (div > 0) {
##         square.side <- sqrt(vp$yscale[2] * vp$xscale[2] / div)

##         ii <- seq(0, vp$xscale[2], by = square.side)
##         jj <- seq(0, vp$yscale[2], by = square.side)

##         grid.segments(x0 = ii, x1 = ii, y0 = 0, y1 = vp$yscale[2],
##                       default.units = "native", gp = gpobj)
##         grid.segments(x0 = 0, x1 = vp$xscale[2], y0 = jj, y1 = jj,
##                       default.units = "native", gp = gpobj)
##       }
##       grid.rect(name = paste(prefix, "rect:", mnames[i], sep = ""),
##                 gp = gpar(fill = "transparent"))
##     }
##   }
## }
##class(struc_sieve) <- "grapcon_generator"

struc_sieve <- function(sievetype = c("observed", "expected"),
                        gp_tile = gpar(), scale = 1) {
    sievetype = match.arg(sievetype)
    function(residuals, observed, expected, spacing,
             gp, split_vertical, prefix = "") {
        observed <- scale * observed
        expected <- scale * expected
        if (is.null(expected)) stop("Need expected values.")
        dn <- dimnames(expected)
        dnn <- names(dn)
        dx <- dim(expected)
        dl <- length(dx)
        n <- sum(expected)

        ## split workhorse
        split <- function(x, i, name, row, col, rowmargin, colmargin, index) {
            cotab <- co_table(x, 1)
            margin <- sapply(cotab, sum)
            v <- split_vertical[i]
            d <- dx[i]

            ## compute total cols/rows and build split layout
            dist <- if (d > 1)
                unit.c(unit(margin, "null"), spacing[[i]])
            else
                unit(margin, "null")
            idx <- matrix(1:(2 * d), nrow = 2, byrow = TRUE)[-2 * d]
            layout <- if (v)
                grid.layout(ncol = 2 * d - 1, widths = dist[idx])
            else
                grid.layout(nrow = 2 * d - 1, heights = dist[idx])
            pushViewport(viewport(layout.pos.col = col, layout.pos.row = row,
                                  layout = layout, name = paste(prefix, "cell:",
                                                   remove_trailing_comma(name), sep = "")))

            ## next level: either create further splits, or final viewports
            row <- col <- rep.int(1, d)
            if (v) col <- 2 * 1:d - 1 else row <- 2 * 1:d - 1
            proptab <- function(x) x / max(sum(x), 1)
            for (m in 1:d) {
                nametmp <- paste(name, dnn[i], "=", dn[[i]][m], ",", sep = "")
                if (v) {
                    colmargintmp <- colmargin * proptab(margin)[m]
                    rowmargintmp <- rowmargin
                } else {
                    rowmargintmp <- rowmargin * proptab(margin)[m]
                    colmargintmp <- colmargin
                }
                if (i < dl)
                    split(cotab[[m]], i + 1, nametmp, row[m], col[m],
                          colmargin = colmargintmp,
                          rowmargin = rowmargintmp, index = cbind(index, m))
                else {
                    pushViewport(viewport(layout.pos.col = col[m],
                                          layout.pos.row = row[m],
                                          name = paste(prefix, "cell:",
                                              remove_trailing_comma(nametmp),
                                              sep = ""),
                                          yscale = c(0, rowmargintmp),
                                          xscale = c(0, colmargintmp)))

                    gpobj <- structure(lapply(gp,
                                              function(x) x[cbind(index, m)]),
                                       class = "gpar")

                    ## draw sieve
                    div <- if (sievetype == "observed")
                        observed[cbind(index, m)]
                    else
                        expected[cbind(index, m)]
                    gptmp <- gp_tile
                    gptmp$col <- "transparent"
                    grid.rect(name = paste(prefix, "rect:", nametmp, sep = ""),
                              gp = gptmp)
                    if (div > 0) {
                        square.side <- sqrt(colmargintmp * rowmargintmp / div)

                        ii <- seq(0, colmargintmp, by = square.side)
                        jj <- seq(0, rowmargintmp, by = square.side)

                        grid.segments(x0 = ii, x1 = ii, y0 = 0, y1 = rowmargintmp,
                                      default.units = "native", gp = gpobj)
                        grid.segments(x0 = 0, x1 = colmargintmp, y0 = jj, y1 = jj,
                                      default.units = "native", gp = gpobj)
                    }
                    gptmp <- gp_tile
                    gptmp$fill <- "transparent"
                    grid.rect(name = paste(prefix, "border:", nametmp, sep = ""),
                              gp = gptmp)

                }
                upViewport(1)
            }
        }

        ## start splitting on top, creates viewport-tree
        split(expected + .Machine$double.eps,
              i = 1, name = "", row = 1, col = 1,
              rowmargin = n, colmargin = n, index = cbind())
    }
}
class(struc_sieve) <- "grapcon_generator"
