#################################################################333
## assocplot

assoc <- function(x, ...)
    UseMethod("assoc")

assoc.formula <-
    function(formula, data = NULL, ..., subset = NULL, na.action = NULL,
             main = NULL, sub = NULL)
{
    if (is.logical(main) && main)
        main <- deparse(substitute(data))
    else if (is.logical(sub) && sub)
        sub <- deparse(substitute(data))

    assoc.default(structable(formula, data, subset = subset,
                             na.action = na.action),
                  main = main, sub = sub, ...)
}

assoc.default <- function(x,
                          row_vars = NULL, col_vars = NULL,
                          compress = TRUE, xlim = NULL, ylim = NULL,
                          spacing = spacing_conditional(sp = 0),
                          spacing_args = list(),
                          split_vertical = NULL,
                          keep_aspect_ratio = FALSE,
                          xscale = 0.9, yspace = unit(0.5, "lines"),
                          main = NULL,
                          sub = NULL,
                          ...,
			  residuals_type = "Pearson",
                          gp_axis = gpar(lty = 3)
                          ) {

    if (is.logical(main) && main)
        main <- deparse(substitute(x))
    else if (is.logical(sub) && sub)
        sub <- deparse(substitute(x))

    if (!inherits(x, "ftable"))
        x <- structable(x)

    tab <- as.table(x)
    dl <- length(dim(tab))

    ## spacing
    cond <- rep(TRUE, dl)
    cond[length(attr(x, "row.vars")) + c(0, length(attr(x, "col.vars")))] <- FALSE
    if (inherits(spacing, "grapcon_generator"))
        spacing <- do.call("spacing", spacing_args)
    spacing <- spacing(dim(tab), condvars = which(cond))

    ## splitting arguments
    if (is.null(split_vertical))
        split_vertical <- attr(x, "split_vertical")

    if(match.arg(tolower(residuals_type), "pearson") != "pearson")
        warning("Only Pearson residuals can be visualized with association plots.")

    strucplot(tab,
              spacing = spacing,
              split_vertical = split_vertical,
              core = struc_assoc(compress = compress, xlim = xlim, ylim = ylim,
              yspace = yspace, xscale = xscale, gp_axis = gp_axis),
              keep_aspect_ratio = keep_aspect_ratio,
              residuals_type = "Pearson",
              main = main,
              sub = sub,
              ...)
}

## old code: more elegant conceptually, but less performant
##
## struc_assoc2 <- function(compress = TRUE, xlim = NULL, ylim = NULL,
##                         yspace = unit(0.5, "lines"), xscale = 0.9,
##                         gp_axis = gpar(lty = 3))
##   function(residuals, observed = NULL, expected, spacing, gp, split_vertical, prefix = "") {
##     dn <- dimnames(expected)
##     dnn <- names(dn)
##     dx <- dim(expected)
##     dl <- length(dx)

##     ## axis limits
##     resid <- structable(residuals, split_vertical = split_vertical)
##     sexpected <- structable(sqrt(expected), split_vertical = split_vertical)
##     rfunc <- function(x) c(min(x, 0), max(x, 0))
##     if (is.null(ylim))
##       ylim <- if (compress)
##         matrix(apply(as.matrix(resid), 1, rfunc), nrow = 2)
##       else
##         rfunc(as.matrix(resid))
##     if (!is.matrix(ylim))
##       ylim <- matrix(as.matrix(ylim), nrow = 2, ncol = nrow(as.matrix(resid)))

##     attr(ylim, "split_vertical") <- rep(TRUE, sum(!split_vertical))
##     attr(ylim, "dnames") <- dn[!split_vertical]
##     class(ylim) <- "structable"

##     if(is.null(xlim))
##       xlim <- if (compress)
##         matrix(c(-0.5, 0.5) %o% apply(as.matrix(sexpected), 2, max), nrow = 2)
##       else
##         c(-0.5, 0.5) * max(sexpected)
##     if (!is.matrix(xlim))
##       xlim <- matrix(as.matrix(xlim), nrow = 2, ncol = ncol(as.matrix(resid)))
##     attr(xlim, "split_vertical") <- rep(TRUE, sum(split_vertical))
##     attr(xlim, "dnames") <- dn[split_vertical]
##     class(xlim) <- "structable"

##     ## split workhorse
##     split <- function(res, sexp, i, name, row, col) {
##       v <- split_vertical[i]
##       splitbase <- if (v) sexp else res
##       splittab <- lapply(seq(dx[i]), function(j) splitbase[[j]])
##       len <- sapply(splittab, function(x) sum(unclass(x)[1,] - unclass(x)[2,]))

##       d <- dx[i]

##       ## compute total cols/rows and build split layout
##       dist <- unit.c(unit(len, "null"), spacing[[i]]  + (1 * !v) * yspace)
##       idx <- matrix(1:(2 * d), nrow = 2, byrow = TRUE)[-2 * d]
##       layout <- if (v)
##         grid.layout(ncol = 2 * d - 1, widths = dist[idx])
##       else
##         grid.layout(nrow = 2 * d - 1, heights = dist[idx])
##       vproot <- viewport(layout.pos.col = col, layout.pos.row = row,
##                          layout = layout, name = remove_trailing_comma(name))

##       ## next level: either create further splits, or final viewports
##       name <- paste(name, dnn[i], "=", dn[[i]], ",", sep = "")
##       rows <- cols <- rep.int(1, d)
##       if (v) cols <- 2 * 1:d - 1 else rows <- 2 * 1:d - 1

##       f <- if (i < dl) {
##         if (v)
##           function(m) split(res, splittab[[m]], i + 1, name[m], rows[m], cols[m])
##         else
##           function(m) split(splittab[[m]], sexp, i + 1, name[m], rows[m], cols[m])
##       } else {
##         if (v)
##           function(m) viewport(layout.pos.col = cols[m], layout.pos.row = rows[m],
##                                name = remove_trailing_comma(name[m]),
##                                yscale = unclass(res)[,1],
##                                xscale = unclass(sexp)[,m], default.units = "null")
##         else
##           function(m) viewport(layout.pos.col = cols[m], layout.pos.row = rows[m],
##                                name = remove_trailing_comma(name[m]),
##                                yscale = unclass(res)[,m],
##                                xscale = unclass(sexp)[,1], default.units = "null")
##       }
##       vpleaves <- structure(lapply(1:d, f), class = c("vpList", "viewport"))

##       vpTree(vproot, vpleaves)
##     }

##     ## start spltting on top, creates viewport-tree
##     pushViewport(split(ylim, xlim, i = 1, name = paste(prefix, "cell:", sep = ""),
##                        row = 1, col = 1))

##     ## draw tiles
##     mnames <- paste(apply(expand.grid(dn), 1,
##                           function(i) paste(dnn, i, collapse = ",", sep = "=")
##                           )
##                     )
##     for (i in seq_along(mnames)) {
##       seekViewport(paste(prefix, "cell:", mnames[i], sep = ""))
##       grid.lines(y = unit(0, "native"), gp = gp_axis)
##       grid.rect(y = 0, x = 0,
##                 height = residuals[i],
##                 width = xscale * unit(sqrt(expected[i]), "native"),
##                 default.units = "native",
##                 gp = structure(lapply(gp, function(x) x[i]), class = "gpar"),
##                 just = c("center", "bottom"),
##                 name = paste(prefix, "rect:", mnames[i], sep = "")
##                 )
##     }

##   }
## class(struc_assoc2) <- "grapcon_generator"


struc_assoc <- function(compress = TRUE, xlim = NULL, ylim = NULL,
                        yspace = unit(0.5, "lines"), xscale = 0.9,
                        gp_axis = gpar(lty = 3))
    function(residuals, observed = NULL, expected, spacing,
             gp, split_vertical, prefix = "") {
        if(is.null(expected)) stop("Need expected values.")
        dn <- dimnames(expected)
        dnn <- names(dn)
        dx <- dim(expected)
        dl <- length(dx)

        ## axis limits
        resid <- structable(residuals, split_vertical = split_vertical)
        sexpected <- structable(sqrt(expected), split_vertical = split_vertical)
        rfunc <- function(x) c(min(x, 0), max(x, 0))
        if (is.null(ylim))
            ylim <- if (compress)
                matrix(apply(as.matrix(resid), 1, rfunc), nrow = 2)
            else
                rfunc(as.matrix(resid))
        if (!is.matrix(ylim))
            ylim <- matrix(as.matrix(ylim), nrow = 2, ncol = nrow(as.matrix(resid)))
        ylim[2,] <- ylim[2,] + .Machine$double.eps

        attr(ylim, "split_vertical") <- rep(TRUE, sum(!split_vertical))
        attr(ylim, "dnames") <- dn[!split_vertical]
        class(ylim) <- "structable"

        if(is.null(xlim))
            xlim <- if (compress)
                matrix(c(-0.5, 0.5) %o% apply(as.matrix(sexpected), 2, max), nrow = 2)
            else
                c(-0.5, 0.5) * max(sexpected)
        if (!is.matrix(xlim))
            xlim <- matrix(as.matrix(xlim), nrow = 2, ncol = ncol(as.matrix(resid)))
        attr(xlim, "split_vertical") <- rep(TRUE, sum(split_vertical))
        attr(xlim, "dnames") <- dn[split_vertical]
        class(xlim) <- "structable"

        ## split workhorse
        split <- function(res, sexp, i, name, row, col, index) {
            v <- split_vertical[i]
            d <- dx[i]

            splitbase <- if (v) sexp else res
            splittab <- lapply(seq(d), function(j) splitbase[[j]])
            len <- abs(sapply(splittab, function(x) sum(unclass(x)[1,] - unclass(x)[2,])))

            ## compute total cols/rows and build split layout
            dist <- if (d > 1)
                unit.c(unit(len, "null"), spacing[[i]]  + (1 * !v) * yspace)
            else
                unit(len, "null")
            idx <- matrix(1:(2 * d), nrow = 2, byrow = TRUE)[-2 * d]
            layout <- if (v)
                grid.layout(ncol = 2 * d - 1, widths = dist[idx])
            else
                grid.layout(nrow = 2 * d - 1, heights = dist[idx])
            pushViewport(viewport(layout.pos.col = col, layout.pos.row = row,
                                  layout = layout, name = paste(prefix, "cell:",
                                                   remove_trailing_comma(name), sep = "")))

            ## next level: either create further splits, or final viewports
            rows <- cols <- rep.int(1, d)
            if (v) cols <- 2 * 1:d - 1 else rows <- 2 * 1:d - 1

            for (m in 1:d) {
                nametmp <- paste(name, dnn[i], "=", dn[[i]][m], ",", sep = "")
                if (i < dl) {
                    if (v) sexp <- splittab[[m]] else res <- splittab[[m]]
                    split(res, sexp, i + 1, nametmp, rows[m], cols[m], cbind(index, m))
                } else {
                    pushViewport(viewport(layout.pos.col = cols[m], layout.pos.row = rows[m],
                                          name = paste(prefix, "cell:",
                                          remove_trailing_comma(nametmp), sep = ""),
                                          yscale = unclass(res)[,if (v) 1 else m],
                                          xscale = unclass(sexp)[,if (v) m else 1],
                                          default.units = "npc")
                                 )

                    ## draw tiles
                    grid.lines(y = unit(0, "native"), gp = gp_axis)
                    grid.rect(y = 0, x = 0,
                              height = residuals[cbind(index, m)],
                              width = xscale * unit(sqrt(expected[cbind(index, m)]), "native"),
                              default.units = "native",
                              gp = structure(lapply(gp, function(x) x[cbind(index,m)]),
                              class = "gpar"),
                              just = c("center", "bottom"),
                              name = paste(prefix, "rect:", remove_trailing_comma(nametmp), sep = "")
                              )
                }
                upViewport(1)
            }
        }
        split(ylim, xlim, i = 1, name = "", row = 1, col = 1, index = cbind())
    }
class(struc_assoc) <- "grapcon_generator"
