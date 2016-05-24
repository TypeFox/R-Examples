update.xcplot <-
function (object, xclick, yclick, xc.cond = NULL, ...)
{
    if (dev.cur() != object$device)
        dev.set(object$device)
    screen(n = object$screen, new = FALSE)
    par(usr = object$usr)
    par(mar = object$mar)
    screen(n = object$screen, new = FALSE)
    if (is.null(xc.cond)){
        xclickconv <- grconvertX(xclick, "ndc", "user")
        yclickconv <- grconvertY(yclick, "ndc", "user")
    }
    if (identical(object$plot.type, "histogram")){
        if (is.null(xc.cond)){
            xc.cond.new <- max(min(xclickconv, max(object$xc, na.rm = TRUE)),
                min(object$xc, na.rm = TRUE), na.rm = TRUE)
        } else {
            xc.cond.new <- xc.cond
        }
        if (xc.cond.new != object$xc.cond.old){
            abline(v = object$xc.cond.old, lwd = 2 * object$select.lwd,
                col = "white")
            break4redraw <- which.min(abs(object$histmp$breaks -
                object$xc.cond.old))
            rect(xleft = object$histmp$breaks[break4redraw + c(-1, 0)], xright =
                object$histmp$breaks[break4redraw + c(0, 1)], ybottom = c(0, 0),
                ytop = object$histmp$counts[break4redraw + c(-1, 0)])
            abline(v = xc.cond.new, lwd = object$select.lwd, col =
                object$select.col)
            object$xc.cond.old <- xc.cond.new
        }
    } else if (identical(object$plot.type, "barplot")){
        if (is.null(xc.cond)){
            xc.cond.new <- as.factor(object$factorcoords$level)[which.min(abs(
                xclickconv - object$factorcoords$x))]
        } else {
            xc.cond.new <- xc.cond
        }
        if (xc.cond.new != object$xc.cond.old){
            barindex.old <- levels(object$xc) == object$xc.cond.old
            rect(xleft = object$bartmp$w.l[barindex.old], xright =
                object$bartmp$w.r[barindex.old], ybottom = 0, ytop =
                object$bartmp$height[barindex.old], col = "gray")
            barindex.new <- levels(object$xc) == xc.cond.new
            rect(xleft = object$bartmp$w.l[barindex.new], xright =
                object$bartmp$w.r[barindex.new], ybottom = 0, ytop =
                object$bartmp$height[barindex.new], col = object$select.colour,
                density = -1)
            object$xc.cond.old <- xc.cond.new
        }
    } else if (identical(object$plot.type, "scatterplot")){
        if (is.null(xc.cond)){
            xc.cond.new.x <- max(min(xclickconv, max(object$xc[, 1], na.rm =
                TRUE)), min(object$xc[, 1], na.rm = TRUE), na.rm = TRUE)
            xc.cond.new.y <- max(min(yclickconv, max(object$xc[, 2], na.rm =
                TRUE)), min(object$xc[, 2], na.rm = TRUE), na.rm = TRUE)
            xc.cond.new <- c(xc.cond.new.x, xc.cond.new.y)
        } else {
            xc.cond.new <- xc.cond
        }
        if (any(xc.cond.new != object$xc.cond.old)){
            if (nrow(object$xc) > 2000 && requireNamespace("gplots", quietly =
                TRUE)){
                par(bg = "white")
                dev.hold()
                screen(new = TRUE)
                object <- plotxc(xc = object$xc, xc.cond = xc.cond.new,
                    name = object$name, select.colour =
                    object$select.colour, select.lwd = object$select.lwd,
                    cex.axis = object$cex.axis, cex.lab = object$cex.lab,
                    tck = object$tck)
                dev.flush()
            } else {
            abline(v = object$xc.cond.old[1], h = object$xc.cond.old[2], lwd =
                2 * object$select.lwd, col = "white")
            xrange <- abs(diff(range(object$xc[, 1])))
            yrange <- abs(diff(range(object$xc[, 2])))
            redrawindex.x <- findInterval(object$xc[, 1], object$xc.cond.old[1]
                + xrange * c(-0.125, 0.125) ) == 1
            redrawindex.y <- findInterval(object$xc[, 2], object$xc.cond.old[2]
                + yrange * c(-0.125, 0.125) ) == 1
            points(object$xc[redrawindex.x | redrawindex.y, ], cex =
                object$select.cex)
            box()
            abline(v = xc.cond.new.x, h = xc.cond.new.y, lwd =
                object$select.lwd, col = object$select.colour)
            object$xc.cond.old <- xc.cond.new
            }
        }
    } else if (identical(object$plot.type, "boxplot")){
        if (is.null(xc.cond)){
            xc.cond.new.x <- as.factor(object$factorcoords$level)[
                which.min(abs(xclickconv - object$factorcoords$x))]
            xc.cond.new.y <- if (abs(yclickconv - object$xc.cond.old[, 2]) >
                0.025 * abs(diff(range(object$xc[, 2])))){
                max(min(yclickconv, max(object$xc[, 2], na.rm = TRUE)),
                min(object$xc[, 2], na.rm = TRUE), na.rm = TRUE)
            } else object$xc.cond.old[, 2]
            xc.cond.new <- c(xc.cond.new.x, xc.cond.new.y)
        } else {
            xc.cond.new <- xc.cond
        }
        if (any(xc.cond.new != object$xc.cond.old)){
            if (xc.cond.new.x != object$xc.cond.old[, 1]){
                abline(v = as.integer(object$xc.cond.old[, 1]), lwd = 2 *
                    object$select.lwd, col = "white")
            }

            if (xc.cond.new.y != object$xc.cond.old[, 2]) {
                abline(h = object$xc.cond.old[, 2], lwd = 2 * object$select.lwd,
                    col = "white")
            }
            par(new = TRUE)
            bxp(object$boxtmp, xaxt = "n", yaxt = "n")
            abline(v = as.integer(xc.cond.new.x), h = xc.cond.new.y, lwd =
                object$select.lwd, col = object$select.colour)
            xc.cond.new <- data.frame(xc.cond.new.x, xc.cond.new.y)
            names(xc.cond.new) <- names(object$xc.cond.old)
            object$xc.cond.old <- xc.cond.new
        }
    } else if (identical(object$plot.type, "spineplot")){
        if (is.null(xc.cond)){
            sptmp <- object$sptmp
            rectcoords <- data.frame(sptmp$xleft, sptmp$xright, sptmp$ybottom,
                sptmp$ytop)
            if (c(xclickconv, yclickconv) %inrectangle% c(min(sptmp$xleft),
                max(sptmp$xright) , min(sptmp$ybottom), max(sptmp$ytop)) ){
                comb.index <- apply(rectcoords, 1L, `%inrectangle%`, point =
                    c(xclickconv, yclickconv))
                if (any(comb.index)){
                    xc.cond.new <- data.frame(as.factor(sptmp$xnames)[
                        comb.index], as.factor(sptmp$ynames)[comb.index])
                    names(xc.cond.new) <- names(object$xc.cond.old)
                    if (any(xc.cond.new != object$xc.cond.old)){
                        object$xc.cond.old <- xc.cond.new
                        par(bg = "white")
                        screen(new = TRUE)
                        object <- plotxc(xc = object$xc, xc.cond = xc.cond.new,
                            name = object$name, select.colour =
                            object$select.colour, select.lwd = object$select.lwd
                            , cex.axis = object$cex.axis, cex.lab =
                            object$cex.lab, tck = object$tck)
                    }
                }
            }
        } else {
            xc.cond.new <- xc.cond
            names(xc.cond.new) <- names(object$xc.cond.old)
            if (any(xc.cond.new != object$xc.cond.old)){
                object$xc.cond.old <- xc.cond.new
                par(bg = "white")
                screen(new = TRUE)
                object <- plotxc(xc = object$xc, xc.cond = xc.cond.new,
                    name = object$name, select.colour =
                    object$select.colour, select.lwd = object$select.lwd,
                    cex.axis = object$cex.axis, cex.lab = object$cex.lab,
                    tck = object$tck)
            }
        }
    }
    object
}

update.xsplot <-
function (object, xc.cond = NULL, data.colour = NULL, data.order = NULL,
    view3d = NULL, theta3d = NULL, phi3d = NULL, xs.grid = NULL, prednew = NULL,
    ...)
{
    if (dev.cur() != object$device)
        dev.set(object$device)
    par(bg = "white")
    screen(n = object$screen, new = FALSE)
    view3d <- if (!is.null(view3d))
        view3d
    else object$view3d
    par(usr = object$usr)
    par(mar = object$mar)
    xc.cond <- if (!is.null(xc.cond))
        xc.cond
    else object$xc.cond
    data.colour <- if (!is.null(data.colour))
        data.colour
    else object$data.colour
    data.order <- if (!is.null(data.order))
        data.order
    else object$data.order
    theta3d <- if (!is.null(theta3d))
        theta3d
    else object$theta3d
    phi3d <- if (!is.null(phi3d))
        phi3d
    else object$phi3d
    conf <- object$conf
    if (any(xc.cond != object$xc.cond)){
        newdata <- makenewdata(xs = object$xs.grid, xc.cond = xc.cond)
        prednew <- lapply(object$model, predict1, newdata = newdata)
    } else {
        newdata <- object$newdata
        prednew <- object$prednew
    }
    color <- if (is.factor(object$y[, 1L])){
      if (identical(nlevels(object$y[, 1L]), 2L) && inherits(object$model[[1L]],
        "glm")){
        factor2color(as.factor(round(prednew[[1L]])))
      } else factor2color(as.factor(prednew[[1L]]))
    } else cont2color(prednew[[1L]], range(object$y[, 1L]))
    ybg <- if (length(data.order) > 0){
        if (is.factor(object$y[, 1L]))
	        factor2color(object$y[data.order, 1L])
	    else cont2color(object$y[data.order, 1L], range(object$y[, 1L]))
    } else NULL
    arefactorsxs <- vapply(object$xs, is.factor, logical(1L))

    if (object$plot.type %in% c("ff")){
        screen(n = object$screen, new = FALSE)
        dev.hold()
        rect(object$usr[1], object$usr[3], object$usr[2], object$usr[4], col =
            "white", border = NA)
        box()
        if (identical(nlevels(object$y[, 1L]), 2L)){
            if (length(data.order) > 0)
                points.default((as.numeric(object$xs[data.order, 1L])) +
                    rnorm(n = length(data.order), sd = 0.1), (as.integer(
                    object$y[data.order, 1L]) - 1) + rnorm(n = length(
                    data.order), sd = 0.01), col = data.colour[data.order])
            for (i in seq_along(object$model)){
                if ("glm" %in% class(object$model[[i]])){
                    points.default(object$xs.grid[, 1L], prednew[[i]],
                        type = 'l', col = object$model.colour[i],
                        lwd = object$model.lwd[i], lty = object$model.lty[i])
                } else {
                    points.default(object$xs.grid[, 1L], as.numeric(prednew[[i
                        ]]) - 1, type = 'l', col = object$model.colour[i], lwd =
                        object$model.lwd[i], lty = object$model.lty[i])
                        }
                    }
        } else {
            if (length(data.order) > 0)
                points(as.numeric(object$xs[data.order, 1L]), as.integer(
                    object$y[data.order, 1L]), col = data.colour[data.order])
            for (i in seq_along(object$model)){
                points.default(as.numeric(object$xs.grid[, 1L]),
                    as.integer(prednew[[i]]), type = 'l', col =
                    object$model.colour[i], lwd = object$model.lwd[i],
                    lty = object$model.lty[i])
            }
        }
        legend("topright", legend = object$model.name, col = object$model.colour
            , lwd = object$model.lwd, lty = object$model.lty)
        dev.flush()
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
    if (object$plot.type %in% c("cf")){
        screen(n = object$screen, new = FALSE)
        dev.hold()
        rect(object$usr[1], object$usr[3], object$usr[2], object$usr[4], col =
            "white", border = NA)
        box()

        if (length(data.order) > 0)
            points(object$xs[data.order, 1L], object$y[data.order, 1L], col =
                data.colour[data.order])
        if (conf){
            prednew2 <- lapply(object$model, confpred, newdata = newdata)
            for (i in seq_along(object$model)){
                points.default(object$xs.grid[, 1L], prednew[[i]], type = 'l',
                    col = object$model.colour[i], lwd = object$model.lwd[i], lty
                        = object$model.lty[i])
                if (all(c("lwr", "upr") %in% colnames(prednew2[[i]]))){
                    points.default(object$xs.grid[, 1L], prednew2[[i]][, "lwr"]
                        , type = 'l', lty = 2, col = object$model.colour[i],
                        lwd = max(0.8, 0.5 * object$model.lwd[i]))
                    points.default(object$xs.grid[, 1L], prednew2[[i]][, "upr"]
                        , type = 'l', lty = 2, col = object$model.colour[i],
                        lwd = max(0.8, 0.5 * object$model.lwd[i]))
                }
            }
        } else {
            for (i in seq_along(object$model)){
                points.default(object$xs.grid[, 1L], prednew[[i]], type = 'l',
                    col = object$model.colour[i], lwd = object$model.lwd[i], lty
                    = object$model.lty[i])
            }
        }
        legend("topright", legend = object$model.name, col = object$model.colour
            , lwd = object$model.lwd, lty = object$model.lty)
        dev.flush()
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
    if (object$plot.type %in% c("fc")){
        screen(n = object$screen, new = FALSE)
        dev.hold()
        rect(object$usr[1], object$usr[3], object$usr[2], object$usr[4], col =
            "white", border = NA)
        box()
        if (identical(nlevels(object$y[, 1L]), 2L)){
            if (length(data.order) > 0)
		        points.default(object$xs[data.order, 1L], as.integer(object$y[
                    data.order, 1L]) - 1, col = data.colour[data.order])
            for (i in seq_along(object$model)){
                if ("glm" %in% class(object$model[[i]])){
                    points.default(object$xs.grid[, 1L], prednew[[i]],
                        type = 'l', col = object$model.colour[i],
                        lwd = object$model.lwd[i], lty = object$model.lty[i])
                } else {
                    points.default(object$xs.grid[, 1L], as.numeric(prednew[[i]
                        ]) - 1, type = 'l', col = object$model.colour[i], lwd =
                        object$model.lwd[i], lty = object$model.lty[i])
                }
            }
        } else {
            if (length(data.order) > 0)
                points(object$xs[data.order, 1L], as.integer(object$y[data.order
                    , 1L]) , col = data.colour[data.order])
            for (i in seq_along(object$model)){
                points.default(as.numeric(object$xs.grid[, 1L]), as.integer(
                    prednew[[i]]), type = 'l', col = object$model.colour[i],
                    lwd = object$model.lwd[i], lty = object$model.lty[i])
            }
        }
        legend("topright", legend = object$model.name, col = object$model.colour
            , lwd = object$model.lwd, lty = object$model.lty)
        dev.flush()
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
    if (object$plot.type %in% c("cc")){
        screen(n = object$screen, new = FALSE)
        dev.hold()
        rect(object$usr[1], object$usr[3], object$usr[2], object$usr[4], col =
            "white", border = NA)
        box()
        if (length(data.order) > 0)
            points(object$xs[data.order, 1L], object$y[data.order, 1L], col =
                data.colour[data.order])
        for (i in seq_along(object$model)){
                points.default(object$xs.grid[, 1L], prednew[[i]], type = 'l',
                    col = object$model.colour[i], lwd = object$model.lwd[i], lty
                    = object$model.lty[i])
        }
        if (conf){
            prednew2 <- lapply(object$model, confpred, newdata = newdata)
            for (i in seq_along(object$model)){
                points.default(object$xs.grid[, 1L], prednew[[i]], type = 'l',
                    col = object$model.colour[i], lwd = object$model.lwd[i], lty
                    = object$model.lty[i])
                if (all(c("lwr", "upr") %in% colnames(prednew2[[i]]))){
                    points.default(object$xs.grid[, 1L], prednew2[[i]][, "lwr"],
                        type = 'l', lty = 2, col = object$model.colour[i], lwd =
                        max(0.8, 0.5 * object$model.lwd[i]))
                    points.default(object$xs.grid[, 1L], prednew2[[i]][, "upr"],
                        type = 'l', lty = 2, col = object$model.colour[i], lwd =
                        max(0.8, 0.5 * object$model.lwd[i]))
                }
            }
        }
        if (is.numeric(object$xs[, 1L])){
            pos <- if (cor(object$xs, object$y) < 0)
                "topright"
            else "bottomright"
            legend(pos, legend = object$model.name, col = object$model.colour,
                lwd = object$model.lwd, lty = object$model.lty)
            } else {
                legend("topright", legend = object$model.name, col =
                    object$model.colour, lwd = object$model.lwd, lty =
                    object$model.lty)
            }
        dev.flush()
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
    if (object$plot.type %in% c("fff", "cff")){
        screen(n = object$screen, new = FALSE)
	    xrect <- as.integer(object$xs.grid[, 1L])
		yrect <- as.integer(object$xs.grid[, 2L])
		xoffset <- abs(diff(unique(xrect)[1:2])) / 2.1
		yoffset <- abs(diff(unique(yrect)[1:2])) / 2.1
        dev.hold()
		rect(xleft = xrect - xoffset, xright = xrect + xoffset, ybottom = yrect
            - yoffset, ytop = yrect + yoffset, col = color)
        if (length(data.order) > 0)
       	    points(jitter(as.integer(object$xs[data.order, 1L]), amount = 0.6 *
                xoffset), jitter(as.integer(object$xs[data.order, 2L]), amount =
                0.6 * yoffset), bg = ybg, col = data.colour[data.order],
                pch = 21)
        dev.flush()
        object$data.colour <- data.colour
        object$data.order <- data.order
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
    if (object$plot.type %in% c("ffc", "cfc")){
        screen(n = object$screen, new = FALSE)
    	xrect <- object$xs.grid[, !arefactorsxs]
		yrect <- as.integer(object$xs.grid[, arefactorsxs])
		xoffset <- abs(diff(unique(xrect)[1:2])) / 2
		yoffset <- abs(diff(unique(yrect)[1:2])) / 2.1
        dev.hold()
		rect(xleft = xrect - xoffset, xright = xrect + xoffset, ybottom = yrect
            - yoffset, ytop = yrect + yoffset, col = color, border = NA)
        if (length(data.order) > 0)
            points(jitter(object$xs[data.order, !arefactorsxs]), jitter(
                as.integer(object$xs[data.order, arefactorsxs])), bg = ybg,
                col = data.colour[data.order], pch = 21)
        dev.flush()
        object$data.colour <- data.colour
        object$data.order <- data.order
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
    if (object$plot.type %in% c("fcc", "ccc")){
        screen(n = object$screen, new = FALSE)
        dev.hold()
        if (object$view3d & identical(object$plot.type, "ccc")){
            screen(n = object$screen, new = TRUE)
            z <- matrix(prednew[[1L]], ncol = 20L, byrow = FALSE)
            zfacet <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1]
                + z[-nrow(z), -ncol(z)]) / 4
            colorfacet <- cont2color(zfacet, range(object$y[, 1L]))
            par(mar = c(3, 3, 3, 3))
            persp.object <- suppressWarnings(persp(x =
                unique(object$xs.grid[, 1L]), y = unique(object$xs.grid[, 2L]),
                border = rgb(0.3, 0.3, 0.3), lwd = 0.1, z = z, col =
                colorfacet, zlim = range(object$y), xlab = colnames(object$xs)[
                1L], ylab = colnames(object$xs)[2L], zlab = colnames(object$y)[
                1L], d = 10, ticktype = "detailed", main =
                "Conditional expectation", theta = theta3d,
                phi = phi3d))
            if (length(data.order) > 0){
                points(trans3d(object$xs[data.order, 1L], object$xs[data.order,
                    2L], object$y[data.order, 1L], pmat = persp.object),
                    col = data.colour[data.order])
                linestarts <- trans3d(object$xs[data.order, 1L], object$xs[
                    data.order, 2L], object$y[data.order, 1L], pmat =
                    persp.object)
                lineends <- trans3d(object$xs[data.order, 1L], object$xs[
                    data.order, 2L], object$yhat[[1]][data.order], pmat =
                    persp.object)
                segments(x0 = linestarts$x, y0 = linestarts$y, x1 =
                    lineends$x, y1 = lineends$y, col = data.colour[
                    data.order])
            }
            object$data.colour <- data.colour
            object$data.order <- data.order
            object$theta3d <- theta3d
            object$phi3d <- phi3d
        } else {
            xoffset <- abs(diff(unique(object$xs.grid[, 1L])[1:2])) / 2
            yoffset <- abs(diff(unique(object$xs.grid[, 2L])[1:2])) / 2
            rect(xleft = object$xs.grid[, 1L] - xoffset, xright =
                object$xs.grid[, 1L] + xoffset, ybottom = object$xs.grid[, 2L]
                - yoffset, ytop = object$xs.grid[, 2L] + yoffset, col = color,
                border = NA)
            if (length(data.order) > 0)
                points(object$xs[data.order, , drop = FALSE], bg = ybg,
                    col = data.colour[data.order], pch = 21)
            }
        dev.flush()
        object$newdata <- newdata
        object$prednew <- prednew
        return(object)
    }
}
