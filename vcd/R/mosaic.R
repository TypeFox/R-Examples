###########################################################
## mosaicplot

mosaic <- function(x, ...)
    UseMethod("mosaic")

mosaic.formula <-
    function(formula, data = NULL, highlighting = NULL,
             ..., main = NULL, sub = NULL, subset = NULL, na.action = NULL)
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

    dep <- gsub(" ", "", fstr[[1]][1])
    if (is.null(highlighting) && (!dep %in% c("","Freq"))) {
        if (all(varnames == ".")) {
            varnames <- if (is.data.frame(data))
                            colnames(data)
                        else
                            names(dimnames(as.table(data)))
            varnames <- varnames[-which(varnames %in% dep)]
        }

        varnames <- c(varnames, dep)
        highlighting <- length(varnames) + length(condnames)
    }


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
        mosaic.default(dat, main = main, sub = sub, highlighting = highlighting,
                       condvars = if (is.null(condind)) NULL else match(condnames, names(dimnames(dat))), ...)
    } else {
        m <- m[c(1, match(c("formula", "data", "subset", "na.action"), names(m), 0))]
        m[[1]] <- as.name("xtabs")
        m$formula <-
            formula(paste(if("Freq" %in% colnames(data)) "Freq",
                          "~",
                          paste(c(condnames, varnames), collapse = "+")))
        tab <- eval(m, parent.frame())
        mosaic.default(tab, main = main, sub = sub, highlighting = highlighting, ...)
    }
}

mosaic.default <- function(x, condvars = NULL,
                           split_vertical = NULL, direction = NULL,
                           spacing = NULL, spacing_args = list(),
                           gp = NULL, expected = NULL, shade = NULL,
                           highlighting = NULL,
                           highlighting_fill = grey.colors,
                           highlighting_direction = NULL,
                           zero_size = 0.5,
                           zero_split = FALSE,
                           zero_shade = NULL,
                           zero_gp = gpar(col = 0),
                           panel = NULL,
                           main = NULL, sub = NULL, ...) {
    zero_shade <- !is.null(shade) && shade || !is.null(expected) || !is.null(gp)
    if (!is.null(shade) && !shade) zero_shade = FALSE

    if (is.logical(main) && main)
        main <- deparse(substitute(x))
    else if (is.logical(sub) && sub)
        sub <- deparse(substitute(x))

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

    ## highlighting
    if (!is.null(highlighting)) {
        if (is.character(highlighting))
            highlighting <- match(highlighting, names(dimnames(x)))
        if (length(highlighting) > 0) {
            if (is.character(condvars))
                condvars <- match(condvars, names(dimnames(x)))
            perm <- if (length(condvars) > 0)
                     c(condvars, seq(dl)[-c(condvars,highlighting)], highlighting)
                 else
                     c(seq(dl)[-highlighting], highlighting)
            x <- aperm(x, perm)
            split_vertical <- split_vertical[perm]
            if (is.null(spacing))
                spacing <- spacing_highlighting
            if (is.function(highlighting_fill))
                highlighting_fill <- rev(highlighting_fill(dim(x)[dl]))
            if (is.null(gp))
                gp <- gpar(fill = highlighting_fill)
            if (!is.null(highlighting_direction)) {
                split_vertical[dl] <- highlighting_direction %in% c("left", "right")
                if (highlighting_direction %in% c("left", "top")) {
                    ## ugly:
                    tmp <- as.data.frame.table(x)
                    tmp[,dl] <- factor(tmp[,dl], rev(levels(tmp[,dl])))
                    x <- xtabs(Freq ~ ., data = tmp)
                    gp <- gpar(fill = rev(highlighting_fill))
                }
            }
        }
    } else if (!is.null(condvars)) { # Conditioning only
        if (is.character(condvars))
            condvars <- match(condvars, names(dimnames(x)))
        if (length(condvars) > 0) {
            perm <- c(condvars, seq(dl)[-condvars])
            x <- aperm(x, perm)
            split_vertical <- split_vertical[perm]
        }
        if (is.null(spacing))
            spacing <- spacing_conditional
    }

    ## spacing argument
    if (is.null(spacing))
        spacing <- if (dl < 3) spacing_equal else spacing_increase

    strucplot(x,
              condvars = if (is.null(condvars)) NULL else length(condvars),
              core = struc_mosaic(zero_size = zero_size, zero_split = zero_split,
              zero_shade = zero_shade, zero_gp = zero_gp, panel = panel),
              split_vertical = split_vertical,
              spacing = spacing,
              spacing_args = spacing_args,
              gp = gp,
              expected = expected,
              shade = shade,
              main = main,
              sub = sub,
              ...)
}

## old code: more elegant, but less performant
##
## struc_mosaic2 <- function(zero_size = 0.5, zero_split = FALSE,
##                          zero_shade = TRUE, zero_gp = gpar(col = 0))
##   function(residuals, observed, expected = NULL, spacing, gp, split_vertical, prefix = "") {
##     dn <- dimnames(observed)
##     dnn <- names(dn)
##     dx <- dim(observed)
##     dl <- length(dx)

##     ## split workhorse
##     zerostack <- character(0)
##     split <- function(x, i, name, row, col, zero) {
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
##       f <- if (i < dl)
##         function(m) {
##           co <- cotab[[m]]
##           z <- mean(co) <= .Machine$double.eps
##           if (z && !zero && !zero_split) zerostack <<- c(zerostack, name[m])
##           split(co, i + 1, name[m], row[m], col[m], z && !zero_split)
##         }
##       else
##         function(m) {
##           if (cotab[[m]] <= .Machine$double.eps && !zero)
##             zerostack <<- c(zerostack, name[m])
##           viewport(layout.pos.col = col[m], layout.pos.row = row[m],
##                    name = remove_trailing_comma(name[m]))
##         }
##       vpleaves <- structure(lapply(1:d, f), class = c("vpList", "viewport"))

##       vpTree(vproot, vpleaves)
##     }

##     ## start spltting on top, creates viewport-tree
##     pushViewport(split(observed + .Machine$double.eps,
##                        i = 1, name = paste(prefix, "cell:", sep = ""),
##                        row = 1, col = 1, zero = FALSE))

##     ## draw rectangles
##     mnames <-  apply(expand.grid(dn), 1,
##                      function(i) paste(dnn, i, collapse=",", sep = "=")
##                      )
##     zeros <- observed <= .Machine$double.eps

##     ## draw zero cell lines
##     for (i in remove_trailing_comma(zerostack)) {
##       seekViewport(i)
##       grid.lines(x = 0.5)
##       grid.lines(y = 0.5)
##       if (!zero_shade && zero_size > 0) {
##         grid.points(0.5, 0.5, pch = 19, size = unit(zero_size, "char"),
##                     gp = zero_gp,
##                     name = paste(prefix, "disc:", mnames[i], sep = ""))
##         grid.points(0.5, 0.5, pch = 1, size = unit(zero_size, "char"),
##                     name = paste(prefix, "circle:", mnames[i], sep = ""))
##       }
##     }

##     # draw boxes
##     for (i in seq_along(mnames)) {
##       seekViewport(paste(prefix, "cell:", mnames[i], sep = ""))
##       gpobj <- structure(lapply(gp, function(x) x[i]), class = "gpar")
##       if (!zeros[i]) {
##         grid.rect(gp = gpobj, name = paste(prefix, "rect:", mnames[i], sep = ""))
##       } else {
##         if (zero_shade && zero_size > 0) {
##           grid.points(0.5, 0.5, pch = 19, size = unit(zero_size, "char"),
##                       gp = gpar(col = gp$fill[i]),
##                       name = paste(prefix, "disc:", mnames[i], sep = ""))
##           grid.points(0.5, 0.5, pch = 1, size = unit(zero_size, "char"),
##                       name = paste(prefix, "circle:", mnames[i], sep = ""))
##         }
##       }
##     }
##   }
## class(struc_mosaic2) <- "grapcon_generator"

struc_mosaic <- function(zero_size = 0.5, zero_split = FALSE,
                         zero_shade = TRUE, zero_gp = gpar(col = 0),
                         panel = NULL)
    function(residuals, observed, expected = NULL,
             spacing, gp, split_vertical, prefix = "") {
        dn <- dimnames(observed)
        dnn <- names(dn)
        dx <- dim(observed)
        dl <- length(dx)

        zeros <- function(gp, name) {
            grid.lines(x = 0.5)
            grid.lines(y = 0.5)
            if (zero_size > 0) {
                grid.points(0.5, 0.5, pch = 19, size = unit(zero_size, "char"),
                            gp = gp, name = paste(prefix, "disc:", name, sep = ""))
                grid.points(0.5, 0.5, pch = 1, size = unit(zero_size, "char"),
                            name = paste(prefix, "circle:", name, sep = ""))
            }
        }

        ## split workhorse
        zerostack <- character(0)
        split <- function(x, i, name, row, col, zero, index) {
            cotab <- co_table(x, 1)
            margin <- sapply(cotab, sum)
            margin[margin == 0] <- .Machine$double.eps
                                        #      margin <- margin + .Machine$double.eps
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
                                                   remove_trailing_comma(name),
                                                   sep = "")))

            ## next level: either create further splits, or final viewports
            row <- col <- rep.int(1, d)
            if (v) col <- 2 * 1:d - 1 else row <- 2 * 1:d - 1
            for (m in 1:d) {
                nametmp <- paste(name, dnn[i], "=", dn[[i]][m], ",", sep = "")
                if (i < dl) {
                    co <- cotab[[m]]

                    ## zeros
                    z <- mean(co) <= .Machine$double.eps
                    split(co, i + 1, nametmp, row[m], col[m],
                          z && !zero_split, cbind(index, m))
                    if (z && !zero && !zero_split && !zero_shade && (zero_size > 0))
                        zeros(zero_gp, nametmp)
                } else {
                    pushViewport(viewport(layout.pos.col = col[m],
                                          layout.pos.row = row[m],
                                          name = paste(prefix, "cell:",
                                          remove_trailing_comma(nametmp), sep = "")))

                    ## zeros
                    if (cotab[[m]] <= .Machine$double.eps && !zero) {
                        zeros(if (!zero_shade) zero_gp else gpar(col = gp$fill[cbind(index,m)]), nametmp)
                    } else {
                        ## rectangles
                        gpobj <- structure(lapply(gp, function(x) x[cbind(index, m)]),
                                           class = "gpar")
                        nam <- paste(prefix, "rect:",
                                     remove_trailing_comma(nametmp), sep = "")
                        if (!is.null(panel))
                            panel(residuals, observed, expected, c(cbind(index, m)),
                                  gpobj, nam)
                        else
                            grid.rect(gp = gpobj, name = nam)
                    }
                }
                upViewport(1)
            }
        }

        ## start splitting on top, creates viewport-tree
        split(observed, i = 1, name = "", row = 1, col = 1,
              zero = FALSE, index = cbind())
    }
class(struc_mosaic) <- "grapcon_generator"
