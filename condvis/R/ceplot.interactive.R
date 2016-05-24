ceplot.interactive <-
function (data, model, response = NULL, S = NULL, C = NULL, sigma = NULL,
    distance = "euclidean", cex.axis = NULL, cex.lab = NULL, tck = NULL,
    view3d = FALSE, Corder = "default", conf = FALSE, separate = TRUE,
    select.colour = "blue", select.cex = 1)
{
    uniqC <- unique(unlist(C))
    xc.cond <- data[1L, uniqC, drop = FALSE]
    xcplots <- list()
    coords <- matrix(ncol = 4L, nrow = length(C))
    plotlegend <- length(S) == 2
    n.selector.cols <- ceiling(length(C) / 4L)
    selector.colwidth <- 2
    height <- 8
    if (separate){
        width <- height + 0.5 * plotlegend
        opendev(width = width, height = height)
        devexp <- dev.cur()
        close.screen(all.screens = TRUE)

        legendwidth <- 1.15 / height
        xsscreens <- if (plotlegend){
            split.screen(figs = matrix(c(0, 1 - legendwidth, 1 - legendwidth, 1,
                0, 0, 1, 1), ncol = 4))
        } else split.screen()
        if (plotlegend){
            screen(xsscreens[2L])
            xslegend(data[, response], colnames(data)[response])
        }
        screen(xsscreens[1L])
        vw <- visualweight(xc = data[, uniqC, drop = FALSE], xc.cond = xc.cond,
            sigma = sigma, distance = distance)
        par(mar = c(3, 3, 3, 3))
        xsplot <- plotxs1(xs = data[, S, drop = FALSE], data[, response,
            drop = FALSE], xc.cond = xc.cond, model = model, data.colour = rgb(1
            - vw$k, 1 - vw$k, 1 - vw$k), data.order = vw$order, view3d = view3d,
            conf = conf)
        xscoords <- par("fig")
        xcwidth <- selector.colwidth * n.selector.cols
        n.selector.rows <- ceiling(length(C) / n.selector.cols)
        xcheight <- selector.colwidth * n.selector.rows
        opendev(height = xcheight, width = xcwidth)
        devcond <- dev.cur()
        close.screen(all.screens = TRUE)
        xcscreens <- split.screen(c(n.selector.rows, n.selector.cols))
        for (i in seq_along(C)){
            screen(xcscreens[i])
            xcplots[[i]] <- plotxc(xc = data[, C[[i]]], xc.cond = data[1L,
                C[[i]]], name = colnames(data[, C[[i]], drop = FALSE]),
                select.colour = select.colour, select.cex = select.cex)
            coords[i, ] <- par("fig")
        }
    } else {
        width <- height + 0.5 * plotlegend + selector.colwidth * n.selector.cols
        opendev(width = width, height = height)
        close.screen(all.screens = TRUE)
        xcwidth <- selector.colwidth * n.selector.cols / width
        mainscreens <- split.screen(figs = matrix(c(0, 1 - xcwidth, 1 - xcwidth,
            1, 0, 0, 1, 1), ncol = 4L))
        xcscreens <- split.screen(c(4L, n.selector.cols), screen = mainscreens[
            2L])
        for (i in seq_along(C)){
            screen(xcscreens[i])
            xcplots[[i]] <- plotxc(xc = data[, C[[i]]], xc.cond = data[1L,
                C[[i]]], name = colnames(data[, C[[i]], drop = FALSE]),
                select.colour = select.colour, select.cex = select.cex)
            coords[i, ] <- par("fig")
        }
        legendwidth <- 1.15 / height
        xsscreens <- if (plotlegend){
            split.screen(figs = matrix(c(0, 1 - legendwidth, 1 - legendwidth, 1,
                0, 0, 1, 1), ncol = 4), screen = mainscreens[1L])
        } else mainscreens[1L]
        if (plotlegend){
            screen(xsscreens[2L])
            xslegend(data[, response], colnames(data)[response])
        }
        screen(xsscreens[1L])
        vw <- visualweight(xc = data[, uniqC, drop = FALSE], xc.cond = xc.cond,
            sigma = sigma, distance = distance)
        par(mar = c(3, 3, 3, 3))
        xsplot <- plotxs1(xs = data[, S, drop = FALSE], data[, response,
            drop = FALSE], xc.cond = xc.cond, model = model, data.colour = rgb(
            1 - vw$k, 1 - vw$k, 1 - vw$k), data.order = vw$order, view3d =
            view3d, conf = conf)
        xscoords <- par("fig")
        xold <- NULL
        yold <- NULL
    }
    mouseclick <- function (separate = FALSE)
    {
        function (buttons, x, y)
        {
            plotindex <- which(apply(coords, 1, `%inrectangle%`, point =
                c(x, y)))
            if (length(plotindex) > 0 && if(exists("buttons")) 0 %in% buttons){
                xcplots[[plotindex]] <<- update(xcplots[[plotindex]], x, y)
                if (any(xc.cond[, xcplots[[plotindex]]$name] !=
                    xcplots[[plotindex]]$xc.cond.old)){
                    xc.cond[, xcplots[[plotindex]]$name] <<-
                        xcplots[[plotindex]]$xc.cond.old
                    vw <<- visualweight(xc = data[, uniqC, drop = FALSE],
                        xc.cond = xc.cond, sigma = vw$sigma, distance =
                        vw$distance)
                    xsplot <<- update(xsplot, xc.cond = xc.cond, data.colour =
                        rgb(1 - vw$k, 1 - vw$k, 1 - vw$k), data.order =
                        vw$order)
                }
            }
            if (all(!separate, findInterval(x, xscoords[1:2]) == 1, identical(
                xsplot$plot.type, "ccc"), xsplot$view3d, 0 %in% buttons)){
                if (!is.null(xold))
                    xsplot <<- update(xsplot, theta3d = xsplot$theta3d + 1 *
                        (xold > x) - 1 * (xold < x), phi3d = xsplot$phi3d + 1 *
                        (yold > y) - 1 * (yold < y), xs.grid = xsplot$xs.grid,
                        prednew = xsplot$prednew)
                xold <<- x
                yold <<- y
            }
        points(NULL)
        }
    }
    keystroke <- function ()
    {
        function (key)
        {
            if (identical(key, "q")){
                cat("\nInteractive session ended.\n")
                return(invisible(1))
            }
            if (identical(xsplot$plot.type, "ccc") & xsplot$view3d &
                key %in% c("Up", "Down", "Left", "Right")){
                xsplot <<- update(xsplot, theta3d = xsplot$theta3d - 2 *
                    (key == "Right") + 2 * (key == "Left"), phi3d = xsplot$phi3d
                    - 2 * (key == "Up") + 2 * (key == "Down"), xs.grid =
                    xsplot$xs.grid, prednew = xsplot$prednew)
            }
#            if (identical(xsplot$plot.type, "ccc") & identical(key, "3"))
#                xsplot <<- update(xsplot, view3d = !xsplot$view3d)
            if (key %in% c(",", ".")){
                sigma <- vw$sigma + 0.01 * vw$sigma * (key == ".") - 0.01 *
                    vw$sigma * (key == ",")
                vw <<- visualweight(xc = data[, uniqC, drop = FALSE],
                    xc.cond = xc.cond, sigma = sigma, distance = vw$distance)
                xsplot <<- update(xsplot, data.colour = rgb(1 - vw$k, 1 - vw$k,
                    1 - vw$k), data.order = vw$order, xs.grid = xsplot$xs.grid,
                    newdata = xsplot$newdata, prednew = xsplot$prednew)
            }
            if (identical(key, "s")){
                if (separate){
                    filename <- paste("snapshot_", gsub(":", ".", gsub(" ", "_",
                        Sys.time())), c("-expectation.pdf", "-condition.pdf"),
                        sep = "")
                    dev.set(devexp)
                    devexpsize <- dev.size()
                    pdf(file = filename[1L], width = devexpsize[1L], height =
                        devexpsize[2L])
                    close.screen(all.screens = TRUE)
                    xsscreens <- if (plotlegend){
                        split.screen(figs = matrix(c(0, 1 - legendwidth, 1 -
                            legendwidth, 1, 0, 0, 1, 1), ncol = 4L))
                    } else split.screen()
                    if (plotlegend){
                        screen(xsscreens[2L])
                        xslegend(data[, response], colnames(data)[response])
                    }
                    screen(xsscreens[1L])
                    plotxs1(xs = data[, S, drop = FALSE], data[, response,
                        drop = FALSE], xc.cond = xc.cond, model = model, data.colour
                        = rgb(1 - vw$k, 1 - vw$k, 1 - vw$k), data.order = vw$order,
                        view3d = xsplot$view3d, theta3d = xsplot$theta3d, phi3d =
                        xsplot$phi3d, conf = conf)
                    dev.off()
                    cat(paste("\nSnapshot saved: '", filename[1L],"'", sep = ""))
                    dev.set(devcond)
                    devcondsize <- dev.size()
                    pdf(file = filename[2L], width = devcondsize[1L], height =
                        devcondsize[2L])
                    close.screen(all.screens = TRUE)
                    xcscreens <- split.screen(c(n.selector.rows, n.selector.cols))
                    for (i in seq_along(C)){
                        screen(xcscreens[i])
                        plotxc(xc = xcplots[[i]]$xc,
                            xc.cond = xcplots[[i]]$xc.cond.old,
                            name = xcplots[[i]]$name,
                            select.colour = xcplots[[i]]$select.colour,
                            select.cex = xcplots[[i]]$select.cex)
                    }
                    dev.off()
                    cat(paste("\nSnapshot saved: '", filename[2L],"'", sep = ""))
                    cat("\n")
                } else {
                filename <- paste("snapshot_", gsub(":", ".", gsub(" ", "_",
                    Sys.time())), ".pdf", sep = "")
                pdf(file = filename, width = width, height = height)
                close.screen(all.screens = TRUE)
                xcwidth <- selector.colwidth * n.selector.cols / width
                mainscreens <- split.screen(figs = matrix(c(0, 1 - xcwidth, 1 -
                    xcwidth, 1, 0, 0, 1, 1), ncol = 4))
                xcscreens <- split.screen(c(4, n.selector.cols), screen =
                    mainscreens[2L])
                for (i in seq_along(C)){
                    screen(xcscreens[i])
                    plotxc(xc = xcplots[[i]]$xc,
                        xc.cond = xcplots[[i]]$xc.cond.old,
                        name = xcplots[[i]]$name,
                        select.colour = xcplots[[i]]$select.colour,
                        select.cex = xcplots[[i]]$select.cex)
                }
                xsscreens <- if (plotlegend){
                    split.screen(figs = matrix(c(0, 1 - legendwidth, 1 -
                        legendwidth, 1, 0, 0, 1, 1), ncol = 4L), screen =
                        mainscreens[1L])
                } else mainscreens[1L]
                if (plotlegend){
                    screen(xsscreens[2L])
                    xslegend(data[, response], colnames(data)[response])
                }
                screen(xsscreens[1L])
                plotxs1(xs = data[, S, drop = FALSE], data[, response,
                    drop = FALSE], xc.cond = xc.cond, model = model, data.colour
                    = rgb(1 - vw$k, 1 - vw$k, 1 - vw$k), data.order = vw$order,
                    view3d = xsplot$view3d, theta3d = xsplot$theta3d, phi3d =
                    xsplot$phi3d, conf = conf)
                dev.off()
                cat(paste("\nSnapshot saved: '", filename,"'", sep = ""))
                cat("\n")
                }
            }
            points(NULL)
        }
    }
    setGraphicsEventHandlers(
        onMouseDown = mouseclick(separate),
        onMouseMove = mouseclick(separate),
        onKeybd = keystroke())
    getGraphicsEventEnv()
    getGraphicsEvent()
}
