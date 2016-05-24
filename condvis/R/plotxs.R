plotxs <-
function (xs, y, xc.cond, model, model.colour = NULL, model.lwd = NULL, 
    model.lty = NULL, model.name = NULL, yhat = NULL, mar = NULL, 
    data.colour = NULL, data.order = NULL, view3d = FALSE, theta3d = 45, 
    phi3d = 20)
{
    model.colour <- if (is.null(model.colour)) 
        if (requireNamespace("RColorBrewer", quietly = TRUE))
		    RColorBrewer::brewer.pal(n = max(length(model), 3L), name = "Dark2")
		else rainbow(length(model))
    else rep(model.colour, length.out = length(model))
    model.lwd <- if (is.null(model.lwd)) 
        rep(2, length(model))
    else rep(model.lwd, length.out = length(model))
    model.lty <- if (is.null(model.lty)) 
        rep(1, length(model))
    else rep(model.lty, length.out = length(model))
    model.name <- if (is.null(model.name)) 
        vapply(model, function(x) tail(class(x), n = 1L), character(1))
    else model.name
    yhat <- if (is.null(yhat))
        lapply(model, predict1)
    else yhat
    data.colour <- if(is.null(data.colour))
        rep("gray", nrow(xs))
    else data.colour
    data.order <- if(is.null(data.order))
        1:nrow(xs)
    else data.order 
    xs.new <- xs[data.order,, drop = FALSE]
    y.new <- y[data.order,, drop = FALSE]
    yhat.new <- yhat[[1L]][data.order]
    data.colour <- data.colour[data.order]
    
    if (identical(ncol(xs), 2L)){
        xs.grid1 <- if (!is.factor(xs[, 1L]))
            seq(min(xs[, 1L], na.rm = TRUE), max(xs[, 1L], na.rm = TRUE), 
                length.out = if (view3d) {20L} else 50L)
        else as.factor(levels(xs[, 1L]))
        xs.grid2 <- if (!is.factor(xs[, 2L]))
            seq(min(xs[, 2L], na.rm = TRUE), max(xs[, 2L], na.rm = TRUE), 
                length.out = if (view3d) {20L} else 50L)
        else as.factor(levels(xs[, 2L]))
        xs.grid <- data.frame(rep(xs.grid1, by = length(xs.grid2)), 
		                      rep(xs.grid2, each = length(xs.grid1)))
    } else {
        xs.grid <- if (!is.factor(xs[, 1L]))
            data.frame(seq(min(xs[, 1L], na.rm = TRUE), max(xs[, 1L], na.rm = 
                TRUE), length.out = if (view3d) {20L} else 50L))
        else data.frame(as.factor(levels(xs[, 1L])))
    }
    colnames(xs.grid) <- colnames(xs)
    newdata <- makenewdata(xs = xs.grid, xc.cond = xc.cond)
	prednew <- lapply(model, predict1, newdata = newdata)
    if(identical(ncol(xs), 1L)){
        if (is.numeric(y[, 1L])){
            plot((xs[, 1L]), (y[, 1L]) + 10 * diff(range(y[, 1L])), col = NULL, 
            main = "Conditional expectation", xlab = colnames(xs)[1L], 
            ylab = colnames(y)[1L], ylim = range(y[, 1L]))
            if (nrow(xs.new) > 0)
                points(xs.new[, 1L], y.new[, 1L], col = data.colour)
            prednew2 <- lapply(model, confpred, newdata = newdata)    
            for (i in seq_along(model)){
                points.default(xs.grid[, 1L], prednew[[i]], type = 'l',
                    col = model.colour[i], lwd = model.lwd[i], lty = model.lty
                    [i])
                if (all(c("lwr", "upr") %in% colnames(prednew2[[i]]))){
                    points.default(xs.grid[, 1L], prednew2[[i]][, "lwr"], type = 
                        'l', lty = 3, col = model.colour[i], lwd = model.lwd[i])
                    points.default(xs.grid[, 1L], prednew2[[i]][, "upr"], type = 
                        'l', lty = 3, col = model.colour[i], lwd = model.lwd[i])    
                }
            }
            if (is.numeric(xs[, 1L])){
                pos <- if (cor(xs, y) < 0L)
                    "topright"
                else
                    "bottomright"
                legend(pos, legend = model.name, col = model.colour, 
                    lwd = model.lwd, lty = model.lty)
            } else {
                legend("topright", legend = model.name, col = model.colour, 
                    lwd = model.lwd, lty = model.lty)
            }
        } else {
            if (is.factor(y[, 1L])){
                if (identical(nlevels(y[, 1L]), 2L)){
                    plot(cbind(xs, as.numeric(y[, 1L]) + 10), col = data.colour, 
                        main = "Conditional expectation", 
                        ylab = paste("Probability ", colnames(y)[1L], "=", 
                        levels(y[, 1L])[2L]), ylim = c(0, 1))
				    points.default((as.numeric(xs.new[, 1L])), 
                        (as.integer(y.new[, 1L]) - 1), col = data.colour)
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
                    legend("topright", legend = model.name, col = model.colour, 
                        lwd = model.lwd, lty = model.lty)
                } else {
                    plot(range(as.numeric(xs[, 1L])) + c(0, 0.1 * abs(diff(range
                        (as.numeric(xs[, 1L])))) ), range(as.integer(y[, 1L])), 
                        col = NULL, xlab = colnames(xs)[1L], ylab = colnames(y)[
                        1L], yaxt = "n", main = "Conditional expectation", xaxt 
                        = if (is.factor(xs[, 1L])) "n" else NULL) 
                    axis(2, at = 1:nlevels(y[, 1L]), labels = levels(y[, 1L]))
                    if (is.factor(xs[, 1L])) 
                        axis(1, at = 1:nlevels(xs[, 1L]), labels = levels(xs[, 
                            1L]))              
                    if (nrow(xs.new) > 0) 
                        points(as.numeric(xs.new[, 1L]), as.integer(y.new[, 1L])
                            , col = data.colour)  
                    for (i in seq_along(model)){
                        points.default(as.numeric(xs.grid[, 1L]), as.integer(
                            prednew[[i]]), type = 'l', col = model.colour[i], 
                            lwd = model.lwd[i], lty = model.lty[i])
                    }
                }    
            }
        }
    } else {
        legendwidth <- 0.15
	    expectation <- split.screen(figs = matrix(c(0, 1 - legendwidth , 1 - 
            legendwidth , 1, 0, 0, 1, 1), ncol = 4))
		screen(expectation[1L])
        arefactorsxs <- vapply(xs, is.factor, logical(1L))
		fhat <- prednew[[1L]]
		color <- if (is.factor(y[, 1L]))
		    factor2color(fhat)
		else cont2color(fhat, range(y[, 1L]))
        ybg <- if (nrow(xs.new) > 0)
            if (is.factor(y.new[, 1L]))
		        factor2color(y.new[, 1L])
		    else cont2color(y.new[, 1L], range(y[, 1L]))  
        else NULL    
        if (all(arefactorsxs)){
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
            if (nrow(xs.new) > 0) 
          	    points(jitter(as.integer(xs.new[, 1L]), amount = 0.6 * xoffset), 
                    jitter(as.integer(xs.new[, 2L]), amount = 0.6 * yoffset), 
                    bg = ybg, col = data.colour, pch = 21)	 
		    axis(1L, at = unique(xrect), labels = levels(xs[, 1L]), 
                tick = FALSE)
			axis(2L, at = unique(yrect), labels = levels(xs[, 2L]),
                tick = FALSE)
		} else { 
            if (any(arefactorsxs)){ 
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
                if (nrow(xs.new) > 0)  
                    points(jitter(xs.new[, !arefactorsxs]), jitter(as.integer(
                        xs.new[, arefactorsxs])), bg = ybg, col = data.colour, 
                        pch = 21)	 
				axis(2L, at = unique(yrect), labels = levels(xs[, 
                    arefactorsxs]), tick = FALSE)
			} else {
                if (view3d){
                    z <- matrix(fhat, ncol = if (view3d) {20L} else 50L, byrow = 
                        FALSE)
                    zfacet <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + 
                        z[-nrow(z), -ncol(z)]) / 4
                    colorfacet <- cont2color(zfacet, range(y[, 1L]))
                    suppressWarnings(persp(x = unique(xs.grid[, 1L]), y = 
                        unique(xs.grid[, 2L]), border = rgb(0.3, 0.3, 0.3), lwd 
                        = 0.1, z = z, col = colorfacet, zlim = range(y), xlab = 
                        colnames(xs)[1L], ylab = colnames(xs)[2L], zlab = 
                        colnames(y)[1L], d = 10, ticktype = "detailed", main = 
                        "Conditional expectation", theta = theta3d, phi = phi3d
                        )) -> persp.object
                    if (nrow(xs.new) > 0){     
                        points(trans3d(xs.new[, 1], xs.new[, 2], y.new[, 1], 
                            pmat = persp.object), col = data.colour)  
                        linestarts <- trans3d(xs.new[, 1], xs.new[, 2], y.new[, 
                            1], pmat = persp.object)   
                        lineends <- trans3d(xs.new[, 1], xs.new[, 2], yhat.new, 
                            pmat = persp.object) 
                        segments(x0 = linestarts$x, y0 = linestarts$y, x1 = 
                            lineends$x, y1 = lineends$y, col = data.colour)                            
                    }                            
                } else {
                    xoffset <- abs(diff(unique(xs.grid[, 1L])[1:2])) / 2
                    yoffset <- abs(diff(unique(xs.grid[, 2L])[1:2])) / 2
                    plot(range(xs.grid[, 1L]), range(xs.grid[, 2L]), col = NULL, 
                        xlab = colnames(xs)[1L], ylab = colnames(xs)[2L], 
                        main = "Conditional expectation")
                    rect(xleft = xs.grid[, 1L] - xoffset, xright = xs.grid[, 1L] 
                        + xoffset, ybottom = xs.grid[, 2L] - yoffset, ytop = 
                        xs.grid[, 2L] + yoffset, col = color, border = NA)
                    if (nrow(xs.new) > 0)     
                        points(xs.new, bg = ybg, col = data.colour, pch = 21)
                }
            }
        }
		screen(expectation[2L])
        if (is.factor(y[, 1L])){
		    par(mar = c(0, 0, 0, 0))
            legend("left", legend = levels(y[, 1L]), 
                fill = factor2color(as.factor(levels(y[, 1L]))),
                title = colnames(y)[1L], bg = "white")				
		} else {
			par(mar = c(8,2.2,8,2.2))
            fullrange <- abs(diff(range(y[, 1L])))
			yrange <- seq(min(y[, 1L], na.rm = TRUE) - 0.15 * fullrange, 
                max(y[, 1L], na.rm = TRUE) + 0.15 * fullrange, length.out = 50L)
			spacing <- abs(diff(unique(yrange[1:2]))) / 2
            plot(0, 0, xaxt = "n", main = colnames(y)[1L], 
                ylab = "", col = NULL, pch = 16, xlab = "", bty = "n", 
                xlim = c(-0.5, 0.5), ylim = range(y[, 1L]))
            rect(xleft = -0.7, xright = 0.7, ybottom = yrange - spacing, 
			     ytop = yrange + spacing, col = cont2color(yrange, 
                 range(y[, 1L])), border = NA)
		    box()
		}
    }
    
    list(xs = xs, y = y, xc.cond = xc.cond, model = model, model.colour = 
        model.colour, model.lwd = model.lwd, model.lty = model.lty, 
        model.name = model.name, yhat = yhat, mar = mar, 
        data.colour = data.colour, data.order = data.order, view3d = view3d, 
        theta3d = theta3d, phi3d = phi3d)
}
