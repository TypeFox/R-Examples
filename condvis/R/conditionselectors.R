conditionselectors <-
function (Xc, type = "minimal", method = "default", Xc.cond = NULL,
    select.colour = "blue", select.cex = 1, select.lwd = 2, ...)
{
    C <- arrangeC(data = Xc, method = method)
    Xc.cond <- if (is.null(Xc.cond))
        Xc[1, , drop = FALSE]
    else Xc.cond
    if (identical(type, "minimal")){
        xcplots <- list()
        n.selector.cols <- ceiling(length(C) / 4L)
        close.screen(all.screens = T)
        selectors <- split.screen(figs = c(ceiling(length(C) / n.selector.cols),
            n.selector.cols))
        for (i in seq_along(C)){
            screen(selectors[i])
            xcplots[[i]] <- plotxc(xc = Xc[, C[[i]]], xc.cond = Xc.cond[1, C[[i
                ]]], name = C[[i]], select.colour = select.colour, select.cex =
                select.cex, ...)
        }
        output <- list(Xc = Xc, Xc.cond = Xc.cond, xcplots = xcplots,
            screens = selectors, type = type, method = method, select.colour =
            select.colour, select.cex = select.cex)
    } else {
        if (identical(type, "full")){
            factorindex <- vapply(Xc, is.factor, logical(1))
            Xc.num <- vapply(Xc, as.numeric, numeric(nrow(Xc)))
            #Xc.cond <- Xc[1, , drop = FALSE]
            Xc.cond.num <- vapply(Xc.cond, as.numeric, numeric(1L))
            close.screen(all.screens = TRUE)
            scr <- split.screen(c(ncol(Xc) + 2, ncol(Xc) + 2))
            scr2 <- as.vector(matrix(scr, ncol = ncol(Xc) + 2)[c(-1,
                -(ncol(Xc) + 2)), c(-1, -(ncol(Xc) + 2))])
            rows <- rep(1:ncol(Xc), each = ncol(Xc))
            cols <- rep(1:ncol(Xc), ncol(Xc))
            dev.hold()
            for (i in seq_along(scr2)){
                screen(scr2[i])
                par(mar = c(0.1, 0.1, 0.1, 0.1))
                par(mgp = c(3, 0.25, 0.15))
                plot(Xc.num[,cols[i]], Xc.num[,rows[i]], cex = 0.6, xlab = "",
                    ylab = "", xaxt = "n", yaxt = "n", col = if (identical(
                    rows[i], cols[i])) NULL else "black")
                if (!identical(rows[i], cols[i]))
                    abline(v = Xc.cond.num[cols[i]], h = Xc.cond.num[rows[i]],
                        col = select.colour, lwd = select.lwd)
                if (identical(rows[i], 1L) & (2 * (round(cols[i] / 2)) ==
                    cols[i]))
                    axis(3, cex.axis = 0.7, tcl = -0.2)
                if (identical(rows[i], ncol(Xc)) & !(2 * (round(cols[i] / 2)) ==
                    cols[i]))
                    axis(1, cex.axis = 0.7, tcl = -0.2)
                if (identical(cols[i], 1L) & (2 * (round(rows[i] / 2)) ==
                    rows[i]))
                    axis(2, cex.axis = 0.7, tcl = -0.2)
                if (identical(cols[i], ncol(Xc)) & !(2 * (round(rows[i] / 2)) ==
                    rows[i]))
                    axis(4, cex.axis = 0.7, tcl = -0.2)
                if (identical(rows[i], cols[i]))
                    text(x = mean(range(Xc.num[,rows[i]])), y = mean(range(
                        Xc.num[,cols[i]])), labels = colnames(Xc.num)[rows[i]])
            }
            coords <- data.frame(t(vapply(scr2,
                function(i) {
                    screen(i, new = F)
                    par("fig")
                }, numeric(4))))
            names(coords) <- c("xleft", "xright", "ybottom", "ytop")
            coords$xcplots.index <- scr2
            dev.flush()
            output <- list(Xc = Xc, Xc.cond = Xc.cond, rows = rows, cols = cols,
                scr2 = scr2, coords = coords, type = type, method = method,
                select.colour = select.colour, select.cex = select.cex)
        } else {
            if (identical(type, "pcp")){
                factorindex <- vapply(Xc, is.factor, logical(1))
                Xc.num <- vapply(Xc, as.numeric, numeric(nrow(Xc)))
                Xc.cond.num <- vapply(Xc.cond, as.numeric, numeric(1L))
                xcoord <- 1:ncol(Xc)
                ycoord <- (Xc.cond.num - apply(Xc.num, 2L, min))/(apply(Xc.num,
                    2L, max) - apply(Xc.num, 2L, min))
                parcoord(Xc.num, main = "Condition selector")
                points(xcoord, ycoord, col = select.colour, type = "l", lwd =
                    select.lwd)
                points(xcoord, ycoord, col = select.colour, pch = 16)
                output <- list(Xc = Xc, Xc.cond = Xc.cond, xcoord = xcoord,
                    ycoord = ycoord, type = type, method = method, select.colour
                    = select.colour, select.cex = select.cex)
            }
        }
    }

    class(output) <- "conditionselector"
    output
}
