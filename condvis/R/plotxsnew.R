plotxs1 <-
function (xs, y, xc.cond, model, model.colour = NULL, model.lwd = NULL,
    model.lty = NULL, model.name = NULL, yhat = NULL, mar = NULL, data.colour =
    NULL, data.order = NULL, view3d = FALSE, theta3d = 45, phi3d = 20, xs.grid =
    NULL, prednew = NULL, conf = FALSE)
{
    dev.hold()
    if (!(ncol(xs) %in% 1:2))
        stop("xs must be a dataframe with 1 or 2 columns")
    if (ncol(y) != 1)
        stop("y must be a dataframe with 1 column")
    model <- if (!is.list(model))
        list(model)
    else model
    model.colour <- if (is.null(model.colour)){
        if (requireNamespace("RColorBrewer", quietly = TRUE))
		    RColorBrewer::brewer.pal(n = max(length(model), 3L), name = "Dark2")
		else rainbow(max(length(model), 4L))
    } else rep(model.colour, length.out = length(model))
    model.lwd <- if (is.null(model.lwd))
        rep(2, length(model))
    else rep(model.lwd, length.out = length(model))
    model.lty <- if (is.null(model.lty))
        rep(1, length(model))
    else rep(model.lty, length.out = length(model))
    model.name <- if(!is.null(names(model)))
        names(model)
    else seq_along(model)
    data.order <- if (is.null(data.order))
        1:nrow(xs)
    else data.order
    data.colour <- if (is.null(data.colour))
        rep("gray", length(data.order))
    else data.colour
    par(mar = c(5, 4, 3, 2))
    if (identical(ncol(xs), 1L)){
        # xs has one column
        if (is.null(xs.grid)){
            xs.grid <- if (!is.factor(xs[, 1L]))
                data.frame(seq(min(xs[, 1L], na.rm = TRUE), max(xs[, 1L], na.rm
                    = TRUE), length.out = if (view3d) {20L} else 50L))
            else data.frame(as.factor(levels(xs[, 1L])))
            colnames(xs.grid) <- colnames(xs)
        }
        newdata <- makenewdata(xs = xs.grid, xc.cond = xc.cond)
        if (is.null(prednew))
	        prednew <- lapply(model, predict1, newdata = newdata)
        if (is.factor(xs[, 1L])){
            # xs is a factor
            if (is.factor(y[, 1L])){
                # y is factor
                plot.type <- "ff"
                if (identical(nlevels(y[, 1L]), 2L)){
                    plot(unique(xs[, 1L]), rep(-888, length(levels(xs[, 1L]))),
                        col = NULL, main = "Conditional expectation", ylab =
                        paste("Probability ", colnames(y)[1L], "=", levels(y[,
                        1L])[2L]), ylim = c(0, 1))
                    if (length(data.order) > 0)
				        points.default((as.numeric(xs[data.order, 1L])) +
                            rnorm(n = length(data.order), sd = 0.1),
                            (as.integer(y[data.order, 1L]) - 1) + rnorm(n =
                            length(data.order), sd = 0.01), col = data.colour[
                            data.order])
                    for (i in seq_along(model)){
                        if ("glm" %in% class(model[[i]])){
                            points.default(xs.grid[, 1L], prednew[[i]],
                                type = 'l', col = model.colour[i],
                                lwd = model.lwd[i], lty = model.lty[i])
                        } else {
                            points.default(xs.grid[, 1L],
                                as.numeric(prednew[[i]]) - 1, type = 'l',
                                col = model.colour[i], lwd = model.lwd[i],
                                lty = model.lty[i])
                        }
                    }
                } else {
                    plot(range(as.numeric(xs[, 1L])) + c(0, 0.1 * abs(diff(
                        range(as.numeric(xs[, 1L])))) ), range(as.integer(
                        y[, 1L])), col = NULL, xlab = colnames(xs)[1L], ylab =
                        colnames(y)[1L], yaxt = "n", main = "Conditional
                        expectation", xaxt = if (is.factor(xs[, 1L])) "n" else
                        NULL)
                    axis(2, at = 1:nlevels(y[, 1L]), labels = levels(y[, 1L]))
                    if (is.factor(xs[, 1L]))
                        axis(1, at = 1:nlevels(xs[, 1L]), labels = levels(xs[,
                            1L]))
                    if (length(data.order) > 0)
                        points(as.numeric(xs[data.order, 1L]), as.integer(y[
                            data.order, 1L]), col = data.colour[data.order])
                    for (i in seq_along(model)){
                        points.default(as.numeric(xs.grid[, 1L]), as.integer(
                            prednew[[i]]), type = 'l', col = model.colour[i],
                            lwd = model.lwd[i], lty = model.lty[i])
                    }
                }
                legend("topright", legend = model.name, col = model.colour,
                    lwd = model.lwd, lty = model.lty)
            } else {
                # y is continuous
                plot.type <- "cf"
                plot(unique(xs[, 1L]), rep(-888, length(levels(xs[, 1L]))), col
                    = NULL, main = "Conditional expectation", xlab = colnames(
                    xs)[1L], ylab = colnames(y)[1L], ylim = range(y[, 1L]))
                if (length(data.order) > 0)
                    points(xs[data.order, 1L], y[data.order, 1L], col =
                        data.colour[data.order])
                if (conf){
                    prednew2 <- lapply(model, confpred, newdata = newdata)
                    for (i in seq_along(model)){
                        points.default(xs.grid[, 1L], prednew[[i]], type = 'l',
                            col = model.colour[i], lwd = model.lwd[i], lty =
                            model.lty[i])
                        if (all(c("lwr", "upr") %in% colnames(prednew2[[i]]))){
                            points.default(xs.grid[, 1L], prednew2[[i]][, "lwr"]
                                , type = 'l', lty = 2, col = model.colour[i],
                                lwd = max(0.8, 0.5 * model.lwd[i]))
                            points.default(xs.grid[, 1L], prednew2[[i]][, "upr"]
                                , type = 'l', lty = 2, col = model.colour[i],
                                lwd = max(0.8, 0.5 * model.lwd[i]))
                        }
                    }
                } else {
                    for (i in seq_along(model)){
                        points.default(xs.grid[, 1L], prednew[[i]], type = 'l',
                            col = model.colour[i], lwd = model.lwd[i], lty =
                            model.lty[i])
                    }
                }
                legend("topright", legend = model.name, col = model.colour,
                    lwd = model.lwd, lty = model.lty)
            }
        } else {
            #xs is continuous
            if (is.factor(y[, 1L])){
                # y is factor
                plot.type <- "fc"
                if (identical(nlevels(y[, 1L]), 2L)){
                    plot(range(xs[, 1L]) + 0.1 * abs(diff(range(xs[, 1L]))),
                        c(0, 0), col = NULL, main = "Conditional expectation",
                        xlab = colnames(xs)[1L], ylab = paste("Probability ",
                        colnames(y)[1L], "=", levels(y[, 1L])[2L]), ylim =
                        c(0, 1))
                    if (length(data.order) > 0)
				        points.default(xs[data.order, 1L], as.integer(y[
                            data.order, 1L]) - 1, col = data.colour[data.order])
                    for (i in seq_along(model)){
                        if ("glm" %in% class(model[[i]])){
                            points.default(xs.grid[, 1L], prednew[[i]],
                                type = 'l', col = model.colour[i],
                                lwd = model.lwd[i], lty = model.lty[i])
                        } else {
                            points.default(xs.grid[, 1L],
                                as.numeric(prednew[[i]]) - 1, type = 'l',
                                col = model.colour[i], lwd = model.lwd[i],
                                lty = model.lty[i])
                        }
                    }
                } else {
                    plot(range(xs[, 1L]), range(as.integer(y[, 1L])), col = NULL
                        , xlab = colnames(xs)[1L], ylab = colnames(y)[1L], yaxt
                        = "n", main = "Conditional expectation", xaxt = if
                        (is.factor(xs[, 1L])) "n" else NULL)
                    axis(2, at = 1:nlevels(y[, 1L]), labels = levels(y[, 1L]))
                    if (is.factor(xs[, 1L]))
                        axis(1, at = 1:nlevels(xs[, 1L]), labels = levels(xs[,
                            1L]))
                    if (length(data.order) > 0)
                        points(xs[data.order, 1L], as.integer(y[data.order, 1L])
                            , col = data.colour[data.order])
                    for (i in seq_along(model)){
                        points.default(as.numeric(xs.grid[, 1L]), as.integer(
                            prednew[[i]]), type = 'l', col = model.colour[i],
                            lwd = model.lwd[i], lty = model.lty[i])
                    }
                }
                legend("topright", legend = model.name, col = model.colour,
                    lwd = model.lwd, lty = model.lty)
            } else {
                # y is continuous
                plot.type <- "cc"
                plot(range(xs[, 1L]), range(y[, 1L]), col = NULL,
                    main = "Conditional expectation", xlab = colnames(xs)[1L],
                    ylab = colnames(y)[1L], ylim = range(y[, 1L]))
                if (length(data.order) > 0)
                    points(xs[data.order, 1L], y[data.order, 1L], col =
                        data.colour[data.order])
                if (conf){
                    prednew2 <- lapply(model, confpred, newdata = newdata)
                    for (i in seq_along(model)){
                        points.default(xs.grid[, 1L], prednew[[i]], type = 'l',
                            col = model.colour[i], lwd = model.lwd[i], lty =
                            model.lty[i])
                        if (all(c("lwr", "upr") %in% colnames(prednew2[[i]]))){
                            points.default(xs.grid[, 1L], prednew2[[i]][, "lwr"]
                                , type = 'l', col = model.colour[i],
                                lwd = max(0.8, 0.5 * model.lwd[i]), lty = 2)
                            points.default(xs.grid[, 1L], prednew2[[i]][, "upr"]
                                , type = 'l', col = model.colour[i],
                                lwd = max(0.8, 0.5 * model.lwd[i]), lty = 2)
                        }
                    }
                } else {
                    for (i in seq_along(model)){
                        points.default(xs.grid[, 1L], prednew[[i]], type = 'l',
                            col = model.colour[i], lwd = model.lwd[i], lty =
                            model.lty[i])
                    }
                }
                pos <- if (cor(xs, y) < 0)
                    "topright"
                else "bottomright"
                legend(pos, legend = model.name, col = model.colour,
                    lwd = model.lwd, lty = model.lty)
            }
        }
    } else {
        # xs has two columns
        arefactorsxs <- vapply(xs, is.factor, logical(1L))
        if (is.null(xs.grid)){
            xs.grid1 <- if (!is.factor(xs[, 1L]))
                seq(min(xs[, 1L], na.rm = TRUE), max(xs[, 1L], na.rm = TRUE),
                    length.out = if (view3d) {20L} else 50L)
            else as.factor(levels(xs[, 1L]))
            xs.grid2 <- if (!is.factor(xs[, 2L]))
                seq(min(xs[, 2L], na.rm = TRUE), max(xs[, 2L], na.rm = TRUE),
                    length.out = if (view3d) {20L} else 50L)
            else as.factor(levels(xs[, 2L]))
            xs.grid <- data.frame(
                rep(xs.grid1, by = length(xs.grid2)),
		        rep(xs.grid2, each = length(xs.grid1)))
            colnames(xs.grid) <- colnames(xs)
        }
        newdata <- makenewdata(xs = xs.grid, xc.cond = xc.cond)
        if (is.null(prednew))
            prednew <- lapply(model, predict1, newdata = newdata)
		color <- if (is.factor(y[, 1L])){
        if (identical(nlevels(y[, 1L]), 2L) && inherits(model[[1L]], "glm")){
            factor2color(as.factor(round(prednew[[1L]])))
        } else factor2color(as.factor(prednew[[1L]]))
		} else cont2color(prednew[[1L]], range(y[, 1L]))
        ybg <- if (length(data.order) > 0){
            if (is.factor(y[, 1L]))
		        factor2color(y[data.order, 1L])
		    else cont2color(y[data.order, 1L], range(y[, 1L]))
        } else NULL
        if (all(arefactorsxs)){
            # xs are both factors
			xrect <- as.integer(xs.grid[, 1L])
			yrect <- as.integer(xs.grid[, 2L])
			xoffset <- abs(diff(unique(xrect)[1:2])) / 2.1
			yoffset <- abs(diff(unique(yrect)[1:2])) / 2.1
			plot(xrect, yrect, col = NULL, xlab = colnames(xs)[1L],
                ylab = colnames(xs)[2L], xlim = c(min(xrect) - xoffset,
                max(xrect) + xoffset), xaxt = "n", bty = "n", ylim =
                c(min(yrect) - yoffset, max(yrect) + yoffset), yaxt = "n",
                main = "Conditional expectation")
			rect(xleft = xrect - xoffset, xright = xrect + xoffset,
			     ybottom = yrect - yoffset, ytop = yrect + yoffset,
				 col = color)
            if (length(data.order) > 0)
          	    points(jitter(as.integer(xs[data.order, 1L]), amount = 0.6 *
                    xoffset), jitter(as.integer(xs[data.order, 2L]), amount =
                    0.6 * yoffset), bg = ybg, col = data.colour[data.order],
                    pch = 21)
		    axis(1L, at = unique(xrect), labels = levels(xs[, 1L]),
                tick = FALSE)
			axis(2L, at = unique(yrect), labels = levels(xs[, 2L]),
                tick = FALSE)
            if (is.factor(y[, 1L])){
                # y is factor
                plot.type <- "fff"
            } else {
                # y is continuous
                plot.type <- "cff"
            }
        } else {
            if (any(arefactorsxs)){
                # xs is one factor, one continuous
                plot.type <- if (is.factor(y[, 1L]))
                    "ffc" # y is factor
                else "cfc" # y is continuous
    	        xrect <- xs.grid[, !arefactorsxs]
			    yrect <- as.integer(xs.grid[, arefactorsxs])
			    xoffset <- abs(diff(unique(xrect)[1:2])) / 2
			    yoffset <- abs(diff(unique(yrect)[1:2])) / 2.1
		        plot(0, 0, col = NULL, xlab = colnames(xs)[!arefactorsxs],
                    ylab = colnames(xs)[arefactorsxs], xlim = c(min(xrect) -
                    xoffset, max(xrect) + xoffset), bty = "n",
                    main = "Conditional expectation", ylim = c(min(yrect) -
                    yoffset, max(yrect) + yoffset), yaxt = "n")
		        rect(xleft = xrect - xoffset, xright = xrect + xoffset,
		            ybottom = yrect - yoffset, ytop = yrect + yoffset,
			        col = color, border = NA)
                if (length(data.order) > 0)
                    points(jitter(xs[data.order, !arefactorsxs]), jitter(
                        as.integer(xs[data.order, arefactorsxs])), bg = ybg,
                        col = data.colour[data.order], pch = 21)
			    axis(2L, at = unique(yrect), labels = levels(xs[, arefactorsxs
                    ]), tick = FALSE)
            } else {
                # xs are both continuous
                if (is.factor(y[, 1L])){
                    # y is factor
                    plot.type <- "fcc"
                    xoffset <- abs(diff(unique(xs.grid[, 1L])[1:2])) / 2
                    yoffset <- abs(diff(unique(xs.grid[, 2L])[1:2])) / 2
                    plot(range(xs.grid[, 1L]), range(xs.grid[, 2L]), col = NULL,
                        xlab = colnames(xs)[1L], ylab = colnames(xs)[2L],
                        main = "Conditional expectation")
                    rect(xleft = xs.grid[, 1L] - xoffset, xright = xs.grid[, 1L]
                        + xoffset, ybottom = xs.grid[, 2L] - yoffset, ytop =
                        xs.grid[, 2L] + yoffset, col = color, border = NA)
                    if (length(data.order) > 0)
                        points(xs[data.order, , drop = FALSE], bg = ybg, col =
                            data.colour[data.order], pch = 21)
                } else {
                    # y is continuous
                    plot.type <- "ccc"
                    if (view3d){
                        yhat <- if (is.null(yhat))
                            lapply(model[1], predict1)
                        else yhat
                        z <- matrix(prednew[[1L]], ncol = 20L, byrow = FALSE)
                        zfacet <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1]
                            + z[-nrow(z), -ncol(z)]) / 4
                        colorfacet <- cont2color(zfacet, range(y[, 1L]))
                        par(mar = c(3, 3, 3, 3))
                        persp.object <- suppressWarnings(persp(x =
                            unique(xs.grid[, 1L]), y = unique(xs.grid[, 2L]),
                            border = rgb(0.3, 0.3, 0.3), lwd = 0.1, z = z, col =
                            colorfacet, zlim = range(y), xlab = colnames(xs)[
                            1L], ylab = colnames(xs)[2L], zlab = colnames(y)[
                            1L], d = 10, ticktype = "detailed", main =
                            "Conditional expectation", theta = theta3d,
                            phi = phi3d))
                        if (length(data.order) > 0){
                            points(trans3d(xs[data.order, 1L], xs[data.order,
                                2L], y[data.order, 1L], pmat = persp.object),
                                col = data.colour[data.order])
                            linestarts <- trans3d(xs[data.order, 1L], xs[
                                data.order, 2L], y[data.order, 1L], pmat =
                                persp.object)
                            lineends <- trans3d(xs[data.order, 1L], xs[
                                data.order, 2L], yhat[[1]][data.order], pmat =
                                persp.object)
                            segments(x0 = linestarts$x, y0 = linestarts$y, x1 =
                                lineends$x, y1 = lineends$y, col = data.colour[
                                data.order])
                        }
                    } else {
                        xoffset <- abs(diff(unique(xs.grid[, 1L])[1:2])) / 2
                        yoffset <- abs(diff(unique(xs.grid[, 2L])[1:2])) / 2
                        plot(range(xs.grid[, 1L]), range(xs.grid[, 2L]), col =
                            NULL, xlab = colnames(xs)[1L], ylab = colnames(xs)[
                            2L], main = "Conditional expectation")
                        rect(xleft = xs.grid[, 1L] - xoffset, xright = xs.grid[,
                            1L] + xoffset, ybottom = xs.grid[, 2L] - yoffset,
                            ytop = xs.grid[, 2L] + yoffset, col = color, border
                            = NA)
                        if (length(data.order) > 0)
                            points(xs[data.order, , drop = FALSE], bg = ybg,
                                col = data.colour[data.order], pch = 21)
                    }
                }
            }
        }
    }
    dev.flush()
    structure(list(xs = xs, y = y, xc.cond = xc.cond, model = model,
        model.colour = model.colour, model.lwd = model.lwd, model.lty =
        model.lty, model.name = model.name, yhat = yhat, mar = par("mar"),
        data.colour = data.colour, data.order = data.order, view3d = view3d,
        theta3d = theta3d, usr = par("usr"), phi3d = phi3d, plot.type = if
        (exists("plot.type")) plot.type else NULL, screen = screen(), device =
        dev.cur(), xs.grid = xs.grid, newdata = newdata, prednew = prednew,
        xs.grid = xs.grid, conf = conf), class = "xsplot")
}
