MDSgui <-
function () 
{
	tclRequire("Tktable")
tclRequire("BWidget")
    newUser <- function() {
        NewU = tktoplevel()
        tkwm.resizable(NewU, "0", "0")
        tkwm.deiconify(NewU)
        tkwm.title(NewU, "New User")
        tkwm.geometry(NewU, "280x100")
        NewUcanvas = tkcanvas(NewU, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(NewUcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = NewU)
        frameNU <- tkwidget(NewU, "TitleFrame", text = "User Name Entry", 
            background = "white")
        tkplace(frameNU, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheight = 0.96, `in` = NewU)
        tkplace(tklabel(frameNU, text = "Enter your Name", background = "white"), 
            relx = 0.05, rely = 0.3, `in` = frameNU)
        UNinput = tclVar("")
        UNtext = tkentry(NewU, width = 15, textvariable = UNinput)
        tkplace(UNtext, relx = 0.6, rely = 0.3, `in` = frameNU)
        On.Enter <- function() {
            UN = as.character(tclvalue(UNinput))
            if (tclvalue(MGvar$Active.UserName) == "None") {
                tclvalue(MGvar$Active.UserName) <<- UN
                tclvalue(MGvar$UserName) <- paste("User: ", tclvalue(MGvar$Active.UserName))
                tkdestroy(NewU)
            }
            else {
                surechange = tkmessageBox(message = "Are you sure you want to replace current user?", 
                  type = "yesno", default = "no")
                if (as.character(surechange) == "yes") {
                  tclvalue(MGvar$Active.UserName) <<- UN
                  tclvalue(MGvar$UserName) <- paste("User: ", 
                    tclvalue(MGvar$Active.UserName))
                  tkdestroy(NewU)
                }
                else {
                  tkdestroy(NewU)
                }
            }
        }
        tkplace(tkbutton(NewU, text = "Enter", width = 15, command = function() On.Enter()), 
            relx = 0.3, rely = 0.7, `in` = NewU)
        tkfocus(NewU)
        tkbind(NewU, "<Return>", On.Enter)
        tkwait.window(NewU)
    }
    plotting2D <- function(data, title = MGvar$activeplot.title, 
        Measure = MGvar$dMeas, showtitle = MGvar$activeplot.title.show, 
        showmeas = MGvar$activeplot.distmeas, xlabel = MGvar$activeplot.xlab, 
        ylabel = MGvar$activeplot.ylab, bgcol = MGvar$activeplot.bg, 
        pointcex = MGvar$activeplot.cex, showlabs = MGvar$activeplot.labs, 
        showpoints = MGvar$activeplot.showpoints, pointcol = MGvar$activeplot.pointcol, 
        pointshape = MGvar$activeplot.type, ymeas = MGvar$activeplot.yaxt, 
        xmeas = MGvar$activeplot.xaxt, axcol = MGvar$activeplot.axescol, 
        indexLabeled = MGvar$indexLabeled, zoomedcoords = MGvar$newCoords, 
        showreg = MGvar$activeplot.showreg, regcol = MGvar$activeplot.regcol, 
        showleg = MGvar$activeplot.showleg, Zrat = MGvar$zoominrat, 
        Mvup = MGvar$moveup, Mvdn = MGvar$movedown, Mvlt = MGvar$moveleft, 
        Mvrt = MGvar$moveright, PTcolsindex = MGvar$MDSmat.Cols, 
        showdist = MGvar$activeplot.showdist, distcol = MGvar$activeplot.distcol) {
        if (showtitle == "yes") {
            graphtitle = title
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showmeas == "yes") {
            distanceMeasure = MGvar$dMeas
        }
        if (showmeas == "no") {
            distanceMeasure = ""
        }
        if (showpoints == "yes") {
            pointtype = "p"
        }
        if (showpoints == "no") {
            pointtype = "n"
        }
        params <- par(bg = bgcol)
        par(mar = c(3, 3, 3, 3))
        par(cex.axis = 0.8)
        if (nrow(data) == 1 && ncol(data) == 1) {
            plot(data, type = "n", xaxt = "n", yaxt = "n", ylab = "n", 
                xlab = "n")
            par(params)
        }
        else {
            MGvar$X.Coords <<- data[, 1]
            MGvar$Y.Coords <<- data[, 2]
            tabzoomswitch <- "off"
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                if (MGvar$Tab1.zoomedswitch == "on") {
                  tabzoomswitch <- "on"
                  MGvar$X.Coords <<- zoomedcoords[, 1]
                  MGvar$Y.Coords <<- zoomedcoords[, 2]
                }
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                if (MGvar$Tab2.zoomedswitch == "on") {
                  tabzoomswitch <- "on"
                  MGvar$X.Coords <<- zoomedcoords[, 1]
                  MGvar$Y.Coords <<- zoomedcoords[, 2]
                }
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                if (MGvar$Tab3.zoomedswitch == "on") {
                  tabzoomswitch <- "on"
                  MGvar$X.Coords <<- zoomedcoords[, 1]
                  MGvar$Y.Coords <<- zoomedcoords[, 2]
                }
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                if (MGvar$Tab4.zoomedswitch == "on") {
                  tabzoomswitch <- "on"
                  MGvar$X.Coords <<- zoomedcoords[, 1]
                  MGvar$Y.Coords <<- zoomedcoords[, 2]
                }
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                if (MGvar$Tab5.zoomedswitch == "on") {
                  tabzoomswitch <- "on"
                  MGvar$X.Coords <<- zoomedcoords[, 1]
                  MGvar$Y.Coords <<- zoomedcoords[, 2]
                }
            }
            Xmin <- min(MGvar$X.Coords)
            Xmax <- max(MGvar$X.Coords)
            Ymin <- min(MGvar$Y.Coords)
            Ymax <- max(MGvar$Y.Coords)
            Xrange <- Xmax - Xmin
            Yrange <- Ymax - Ymin
            Xmin <- Xmin + (1 - Zrat) * Xrange/2 - (Mvlt - 1) * 
                Xrange + (Mvrt - 1) * Xrange
            Xmax <- Xmax - (1 - Zrat) * Xrange/2 - (Mvlt - 1) * 
                Xrange + (Mvrt - 1) * Xrange
            Ymin <- Ymin + (1 - Zrat) * Yrange/2 - (Mvdn - 1) * 
                Yrange + (Mvup - 1) * Yrange
            Ymax <- Ymax - (1 - Zrat) * Yrange/2 - (Mvdn - 1) * 
                Yrange + (Mvup - 1) * Yrange
            if (showleg == "yes") {
                Ymax = Ymax + 0.03 * Yrange
                Xmax = Xmax + 0.015 * Xrange
                Xmin = Xmin - 0.015 * Xrange
            }
            plot(MGvar$X.Coords, MGvar$Y.Coords, ylim = c(Ymin, 
                Ymax), xlim = c(Xmin, Xmax), type = "n", asp = 1, 
                xaxt = xmeas, yaxt = ymeas, main = graphtitle, 
                bg = bgcol, cex = pointcex, col = pointcol, fg = axcol, 
                pch = pointshape, ylab = "", xlab = "", plt = c(1, 
                  1, 1, 1))
            par(params)
            countcols = 1
            if (length(MGvar$remPcompindex) > 0) {
                for (i in 1:nrow(MGvar$activedata)) {
                  skip = FALSE
                  for (j in 1:length(MGvar$remPcompindex)) {
                    if (i == MGvar$remPcompindex[j]) {
                      skip = TRUE
                    }
                  }
                  if (!skip) {
                    MGvar$PTcolsindex[countcols] = MGvar$MDSmat.Cols[i]
                    countcols = countcols + 1
                  }
                }
            }
            if (showpoints == "yes") {
                for (i in 1:nrow(data)) {
                  points(MGvar$X.Coords[i], MGvar$Y.Coords[i], 
                    type = "p", cex = pointcex, col = PTcolsindex[i], 
                    pch = pointshape)
                }
            }
            mtext(text = xlabel, side = 1, line = 1.9, adj = 0.5, 
                cex = 1, col = "black")
            mtext(text = ylabel, side = 2, line = 1.9, adj = 0.5, 
                cex = 1, col = "black")
            if (showleg == "yes") {
                legend("topleft", pch = 15, bty = "o", col = MGvar$ClasTabCols, 
                  legend = names(MGvar$ClasTabCols))
            }
            if (showreg == "yes") {
                regdat = MGvar$activedata
                if (MGvar$removedpoints) {
                  regdat = MGvar$removedpointsactivedata
                }
                X = data
                labels = colnames(regdat)
                var.names <- labels
                Z.info <- biplot.check.Z(regdat, FALSE)
                Z <- Z.info$X
                unscaled.Z <- Z.info$unscaled.X
                means <- Z.info$means
                sd <- Z.info$sd
                Z <- as.matrix(Z)
                Bmat <- solve(t(X) %*% X) %*% t(X) %*% Z
                Br <- t(Bmat)
                usr <- par("usr")
                numax <- ncol(Z)
                ax.style <- biplot.ax.control(numax, list(var.names))
                ax <- ax.style
                if (ax$type == "prediction") 
                  if (nrow(Br) > 1) 
                    axes.direction <- solve(diag(diag(Br %*% 
                      t(Br)))) %*% Br
                  else axes.direction <- (1/(Br %*% t(Br))) %*% 
                    Br
                else axes.direction <- Br
                z.axes <- lapply(1:length(ax$which), calibrate.axis, 
                  unscaled.Z, means, sd, axes.direction, ax$which, 
                  ax$ticks, ax$orthogx, ax$orthogy, ax$oblique)
                colrs = c()
                for (i in 1:ncol(regdat)) {
                  if (i%%7 == 0) {
                    num = 7
                  }
                  if (i%%7 != 0) {
                    num = i%%7
                  }
                  colrs[i] = brewer.pal(7, "Dark2")[num]
                }
                regcurves <- vector("list", ncol(regdat))
                for (i in 1:ncol(regdat)) {
                  nextrun = FALSE
                  if (length(MGvar$remAxcompindex) > 0) {
                    for (j in 1:length(MGvar$remAxcompindex)) {
                      if (i == MGvar$remAxcompindex[j]) {
                        nextrun = TRUE
                      }
                    }
                  }
                  if (nextrun) {
                    next
                  }
                  y = regdat[, i]
                  B = solve(t(X) %*% X) %*% t(X) %*% y
                  m = B[2]/B[1]
                  grandmax = max(Ymax, Xmax)
                  grandmin = min(Ymin, Xmin)
                  if (m > 1 || m < -1) {
                    y1 = grandmin * 2
                    y2 = grandmax * 2
                    x1 = y1/m
                    x2 = y2/m
                  }
                  if (m <= 1 && m >= -1) {
                    x1 = grandmin * 2
                    x2 = grandmax * 2
                    y1 = m * x1
                    y2 = m * x2
                  }
                  Xvals = c(x1, x2)
                  Yvals = c(y1, y2)
                  lines(Xvals, Yvals, col = colrs[i])
                  theta = (atan(m) * 57.2957795)
                  if (m > 1 || m < -1) {
                    text(x1/2, Ymin, labels[i], srt = theta, 
                      pos = 3, cex = 0.5, col = colrs[i])
                  }
                  if (m <= 1 && m >= -1) {
                    text(Xmin, y1/2, labels[i], srt = theta, 
                      pos = 3, cex = 0.5, col = colrs[i])
                  }
                  marker.mat <- z.axes[[i]][z.axes[[i]][, 4] == 
                    1, 1:3]
                  x.vals <- marker.mat[, 1]
                  y.vals <- marker.mat[, 2]
                  lin.coef <- coefficients(lm(y.vals ~ x.vals))
                  invals <- x.vals < usr[2] & x.vals > usr[1] & 
                    y.vals < usr[4] & y.vals > usr[3]
                  std.markers <- zapsmall(marker.mat[invals, 
                    3])
                  x.vals <- x.vals[invals]
                  y.vals <- y.vals[invals]
                  if (ax.style$tick.label[i]) 
                    label.on.off <- rep(1, sum(invals))
                  else rep(0, sum(invals))
                  if (!ax.style$tick.label[i]) 
                    label.on.off[c(1, length(label.on.off))] <- 1
                  apply(cbind(x.vals, y.vals, std.markers, label.on.off), 
                    1, .marker.func, coef = lin.coef, col = colrs[i], 
                    tick.size = ax.style$tick.size[i], side = ax.style$tick.label.side[i], 
                    pos = ax.style$tick.label.pos[i], offset = ax.style$tick.label.offset[i], 
                    label.col = ax.style$tick.label.col[i], cex = ax.style$tick.label.cex[i])
                }
            }
            if (showlabs == "yes") {
                if (showpoints == "no") {
                  for (i in 1:nrow(data)) {
                    text(MGvar$X.Coords[i], MGvar$Y.Coords[i], 
                      rownames(data)[i], cex = pointcex, col = PTcolsindex[i])
                  }
                }
                else {
                  for (i in 1:nrow(data)) {
                    text(MGvar$X.Coords[i], MGvar$Y.Coords[i], 
                      rownames(data)[i], cex = pointcex, pos = 3)
                  }
                }
            }
            if (showmeas == "yes") {
                mtext(text = Measure, side = 3, line = 2.4, adj = 0.01, 
                  cex = 0.7, col = "black")
            }
            MGvar$parPlotSize <<- par("plt")
            MGvar$usrCoords <<- par("usr")
            MGvar$labelsVec <<- rownames(data)
            MGvar$ZlabelsVec <<- rownames(zoomedcoords)
            if (length(MGvar$indexLabeled) > 0) {
                for (i in (1:length(MGvar$indexLabeled))) {
                  indexClosest <- MGvar$indexLabeled[i]
                  if (tabzoomswitch == "off") {
                    text(MGvar$X.Coords[indexClosest], MGvar$Y.Coords[indexClosest], 
                      labels = MGvar$labelsVec[indexClosest], 
                      pos = MGvar$GenSet.ILPos)
                  }
                  if (tabzoomswitch == "on") {
                    text(MGvar$X.Coords[indexClosest], MGvar$Y.Coords[indexClosest], 
                      labels = MGvar$ZlabelsVec[indexClosest], 
                      pos = MGvar$GenSet.ILPos)
                  }
                }
            }
            Shepindexmat <- matrix(0, nrow = nrow(MGvar$MDSmat), 
                ncol = nrow(MGvar$MDSmat))
            rownames(Shepindexmat) <- rownames(MGvar$MDSmat)
            colnames(Shepindexmat) <- rownames(MGvar$MDSmat)
            count = 1
            for (i in 1:(nrow(MGvar$MDSmat) - 1)) {
                for (j in (i + 1):nrow(MGvar$MDSmat)) {
                  Shepindexmat[i, j] = count
                  count = count + 1
                }
            }
            truepos = 0
            if (showdist == "yes") {
                if (length(MGvar$Shep.indexLabeled) > 0) {
                  cls = rainbow(length(MGvar$Shep.indexLabeled))
                  for (i in (1:length(MGvar$Shep.indexLabeled))) {
                    for (j in 1:length(MGvar$Shepx)) if (MGvar$Shep.indexLabeled[i] == 
                      MGvar$ShepPointindex[j, 1]) {
                      truepos = MGvar$ShepPointindex[j, 2]
                    }
                    for (k in 1:nrow(MGvar$MDSmat)) {
                      for (m in 1:nrow(MGvar$MDSmat)) {
                        if (Shepindexmat[k, m] == truepos) {
                          Shepobj1 = rownames(Shepindexmat)[k]
                          for (p in 1:nrow(MGvar$MDSmat)) {
                            if (Shepobj1 == rownames(MGvar$MDSmat)[p]) {
                              obj1x = data[p, 1]
                              obj1y = data[p, 2]
                            }
                          }
                          Shepobj2 = colnames(Shepindexmat)[m]
                          for (r in 1:nrow(MGvar$MDSmat)) {
                            if (Shepobj2 == rownames(MGvar$MDSmat)[r]) {
                              obj2x = data[r, 1]
                              obj2y = data[r, 2]
                            }
                          }
                          segments(obj1x, obj1y, obj2x, obj2y, 
                            col = cls[i])
                        }
                      }
                    }
                  }
                }
            }
            segments(MGvar$first.xCoord, MGvar$first.yCoord, 
                MGvar$first.xCoord, MGvar$latest.yCoord, col = "blue")
            segments(MGvar$first.xCoord, MGvar$first.yCoord, 
                MGvar$latest.xCoord, MGvar$first.yCoord, col = "blue")
            segments(MGvar$latest.xCoord, MGvar$first.yCoord, 
                MGvar$latest.xCoord, MGvar$latest.yCoord, col = "blue")
            segments(MGvar$first.xCoord, MGvar$latest.yCoord, 
                MGvar$latest.xCoord, MGvar$latest.yCoord, col = "blue")
        }
    }
    .marker.func <- function(vec, coef, col, tick.size, side, 
        pos, offset, label.col, cex) {
        x <- vec[1]
        y <- vec[2]
        marker.val <- vec[3]
        label.on.off <- vec[4]
        if (is.na(coef[2])) 
            .marker.label.cm(x, y, grad = "h", marker.val, expand = tick.size, 
                col = col, label.on.off = label.on.off, side = side, 
                pos = pos, offset = offset, label.col = label.col, 
                cex = cex)
        else if (coef[2] == 0) 
            .marker.label.cm(x, y, grad = "v", marker.val, expand = tick.size, 
                col = col, label.on.off = label.on.off, side = side, 
                pos = pos, offset = offset, label.col = label.col, 
                cex = cex)
        else .marker.label.cm(x, y, grad = -1/coef[2], marker.val, 
            expand = tick.size, col = col, label.on.off = label.on.off, 
            side = side, pos = pos, offset = offset, label.col = label.col, 
            cex = cex)
    }
    calibrate.axis <- function(j, unscaled.X, means, sd, axes.rows, 
        ax.which, ax.tickvec, ax.orthogxvec, ax.orthogyvec, ax.oblique) {
        ax.num <- ax.which[j]
        tick <- ax.tickvec[j]
        ax.direction <- axes.rows[j, ]
        r <- ncol(axes.rows)
        ax.orthog <- rbind(ax.orthogxvec, ax.orthogyvec)
        if (nrow(ax.orthog) < r) 
            ax.orthog <- rbind(ax.orthog, 0)
        if (nrow(axes.rows) > 1) 
            phi.vec <- diag(1/diag(axes.rows %*% t(axes.rows))) %*% 
                axes.rows %*% ax.orthog[, j]
        else phi.vec <- (1/(axes.rows %*% t(axes.rows))) %*% 
            axes.rows %*% ax.orthog[, j]
        number.points <- 100
        std.ax.tick.label <- pretty(unscaled.X[, ax.num], n = tick)
        std.range <- c(min(std.ax.tick.label), max(std.ax.tick.label))
        std.ax.tick.label.min <- std.ax.tick.label - (std.range[2] - 
            std.range[1])
        std.ax.tick.label.max <- std.ax.tick.label + (std.range[2] - 
            std.range[1])
        std.ax.tick.label <- c(std.ax.tick.label, std.ax.tick.label.min, 
            std.ax.tick.label.max)
        interval <- (std.ax.tick.label - means[ax.num])/sd[ax.num]
        axis.vals <- seq(from = min(interval), to = max(interval), 
            length = number.points)
        axis.vals <- sort(unique(c(axis.vals, interval)))
        number.points <- length(axis.vals)
        axis.points <- matrix(0, nrow = number.points, ncol = r)
        for (i in 1:r) axis.points[, i] <- ax.orthog[i, ax.num] + 
            (axis.vals - phi.vec[ax.num]) * ax.direction[i]
        axis.points <- cbind(axis.points, axis.vals * sd[ax.num] + 
            means[ax.num], 0)
        for (i in 1:number.points) if (any(zapsmall(axis.points[i, 
            r + 1] - std.ax.tick.label) == 0)) 
            axis.points[i, r + 2] <- 1
        axis.points
    }
    biplot.ax.control <- function(p, X.names, which = 1:p, type = "prediction", 
        col = "grey", lwd = 1, lty = 1, label = "Orthog", label.col = col, 
        label.cex = 0.75, label.dist = 0, ticks = 5, tick.col = col, 
        tick.size = 1, tick.label = T, tick.label.col = tick.col, 
        tick.label.cex = 0.6, tick.label.side = "left", tick.label.offset = 0.5, 
        tick.label.pos = 1, predict.col = col, predict.lwd = lwd, 
        predict.lty = lty, ax.names = X.names, rotate = NULL, 
        orthogx = 0, orthogy = 0, oblique = NULL) {
        if (!all(is.numeric(which))) 
            which <- match(which, X.names, nomatch = 0)
        which <- which[which <= p]
        which <- which[which > 0]
        ax.num <- length(which)
        if (type != "prediction" & type != "interpolation") 
            stop("Incorrect type of biplot axes specified")
        while (length(col) < ax.num) col <- c(col, col)
        col <- as.vector(col[1:ax.num])
        while (length(lwd) < ax.num) lwd <- c(lwd, lwd)
        lwd <- as.vector(lwd[1:ax.num])
        while (length(lty) < ax.num) lty <- c(lty, lty)
        lty <- as.vector(lty[1:ax.num])
        if (label != "Orthog" & label != "Hor" & label != "Paral") 
            stop("Incorrect specification of axis label direction")
        while (length(label.col) < ax.num) label.col <- c(label.col, 
            label.col)
        label.col <- as.vector(label.col[1:ax.num])
        while (length(label.cex) < ax.num) label.cex <- c(label.cex, 
            label.cex)
        label.cex <- as.vector(label.cex[1:ax.num])
        while (length(label.dist) < ax.num) label.dist <- c(label.dist, 
            label.dist)
        label.dist <- as.vector(label.dist[1:ax.num])
        while (length(ticks) < ax.num) ticks <- c(ticks, ticks)
        ticks <- as.vector(ticks[1:ax.num])
        while (length(tick.col) < ax.num) tick.col <- c(tick.col, 
            tick.col)
        tick.col <- as.vector(tick.col[1:ax.num])
        while (length(tick.size) < ax.num) tick.size <- c(tick.size, 
            tick.size)
        tick.size <- as.vector(tick.size[1:ax.num])
        while (length(tick.label) < ax.num) tick.label <- c(tick.label, 
            tick.label)
        tick.label <- as.vector(tick.label[1:ax.num])
        while (length(tick.label.col) < ax.num) tick.label.col <- c(tick.label.col, 
            tick.label.col)
        tick.label.col <- as.vector(tick.label.col[1:ax.num])
        while (length(tick.label.cex) < ax.num) tick.label.cex <- c(tick.label.cex, 
            tick.label.cex)
        tick.label.cex <- as.vector(tick.label.cex[1:ax.num])
        while (length(tick.label.side) < ax.num) tick.label.side <- c(tick.label.side, 
            tick.label.side)
        tick.label.side <- as.vector(tick.label.side[1:ax.num])
        while (length(tick.label.offset) < ax.num) tick.label.offset <- c(tick.label.offset, 
            tick.label.offset)
        tick.label.offset <- as.vector(tick.label.offset[1:ax.num])
        while (length(tick.label.pos) < ax.num) tick.label.pos <- c(tick.label.pos, 
            tick.label.pos)
        tick.label.pos <- as.vector(tick.label.pos[1:ax.num])
        while (length(predict.col) < ax.num) predict.col <- c(predict.col, 
            predict.col)
        predict.col <- as.vector(predict.col[1:ax.num])
        while (length(predict.lwd) < ax.num) predict.lwd <- c(predict.lwd, 
            predict.lwd)
        predict.lwd <- as.vector(predict.lwd[1:ax.num])
        while (length(predict.lty) < ax.num) predict.lty <- c(predict.lty, 
            predict.lty)
        predict.lty <- as.vector(predict.lty[1:ax.num])
        ax.names <- ax.names[which]
        while (length(ax.names) < p) ax.names <- c(ax.names, 
            "")
        ax.names <- as.vector(ax.names[1:ax.num])
        if (!is.null(oblique)) 
            if (length(oblique) != p) 
                stop("For oblique translations values must be specified for each variable")
        while (length(orthogx) < p) orthogx <- c(orthogx, orthogx)
        orthogx <- as.vector(orthogx[1:p])
        while (length(orthogy) < p) orthogy <- c(orthogy, orthogy)
        orthogy <- as.vector(orthogy[1:p])
        list(which = which, type = type, col = col, lwd = lwd, 
            lty = lty, label = label, label.col = label.col, 
            label.cex = label.cex, label.dist = label.dist, ticks = ticks, 
            tick.col = tick.col, tick.size = tick.size, tick.label = tick.label, 
            tick.label.col = tick.label.col, tick.label.cex = tick.label.cex, 
            tick.label.side = tick.label.side, tick.label.offset = tick.label.offset, 
            tick.label.pos = tick.label.pos, predict.col = predict.col, 
            predict.lty = predict.lty, predict.lwd = predict.lwd, 
            names = ax.names, rotate = rotate, orthogx = orthogx, 
            orthogy = orthogy, oblique = oblique)
    }
    biplot.check.Z <- function(X, scaled.mat) {
        X <- as.matrix(X)
        unscaled.X <- X
        means <- apply(X, 2, mean)
        sd <- sqrt(apply(X, 2, var))
        if (scaled.mat) 
            X <- scale(X)
        else {
            X <- scale(X, scale = FALSE)
            sd <- rep(1, ncol(X))
        }
        if (is.null(dimnames(X))) 
            dimnames(X) <- list(paste(1:nrow(X)), paste("V", 
                1:ncol(X), sep = ""))
        if (length(dimnames(X)[[1]]) == 0) 
            dimnames(X)[[1]] <- paste(1:nrow(X))
        if (length(dimnames(X)[[2]]) == 0) 
            dimnames(X)[[2]] <- paste("V", 1:ncol(X), sep = "")
        list(X = X, unscaled.X = unscaled.X, means = means, sd = sd)
    }
    .marker.label.cm <- function(x, y, grad, marker.val, expand = 1, 
        col, label.on.off, side, pos, offset, label.col, cex) {
        usr <- par("usr")
        uin <- par("pin")/c(usr[2] - usr[1], usr[4] - usr[3])
        mm <- 1/(uin[1] * 25.4)
        d <- expand * mm
        if (grad == "v") {
            lines(rep(x, 2), c(y - d, y + d), col = col)
            if (label.on.off == 1) 
                text(x, y - d, label = marker.val, pos = pos, 
                  offset = offset, col = label.col, cex = cex)
        }
        if (grad == "h") {
            lines(c(x - d, x + d), rep(y, 2), col = col)
            if (label.on.off == 1) 
                if (side == "right") 
                  text(x + d, y, label = marker.val, pos = pos, 
                    offset = offset, col = label.col, cex = cex)
                else text(x - d, y, label = marker.val, pos = pos, 
                  offset = offset, col = label.col, cex = cex)
        }
        if (is.numeric(grad)) {
            b <- d * sqrt(1/(1 + grad * grad))
            a <- b * grad
            lines(c(x - b, x + b), c(y - a, y + a), col = col)
            if (label.on.off == 1) 
                if (side == "right") 
                  text(x + b, y + a, label = marker.val, pos = pos, 
                    offset = offset, col = label.col, cex = cex)
                else text(x - b, y - a, label = marker.val, pos = pos, 
                  offset = offset, col = label.col, cex = cex)
        }
    }
    labelClosestPoint <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance.P <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        indexClosest <- which.min(squared.Distance.P)
        MGvar$indexLabeled <<- c(MGvar$indexLabeled, indexClosest)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$indexLabeled.T1 <<- c(MGvar$indexLabeled.T1, 
                indexClosest)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$indexLabeled.T2 <<- c(MGvar$indexLabeled.T2, 
                indexClosest)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$indexLabeled.T3 <<- c(MGvar$indexLabeled.T3, 
                indexClosest)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$indexLabeled.T4 <<- c(MGvar$indexLabeled.T4, 
                indexClosest)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$indexLabeled.T5 <<- c(MGvar$indexLabeled.T5, 
                indexClosest)
        }
        tabplot()
    }
    OnPlotLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$X.Coords - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$Y.Coords - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        labelClosestPoint(xClick, yClick, imgXcoords, imgYcoords)
    }
    LabelSpecificPoint <- function() {
        Labtt = tktoplevel()
        tkwm.resizable(Labtt, "0", "0")
        tkwm.deiconify(Labtt)
        tkwm.title(Labtt, "Point Label")
        tkwm.geometry(Labtt, "280x100")
        LabCanvas = tkcanvas(Labtt, width = 280, height = 100, 
            bg = col.sec)
        tkplace(LabCanvas, relx = 0, `in` = Labtt)
        frameLab <- tkwidget(Labtt, "TitleFrame", text = "Label a Specific Point", 
            background = "white")
        tkplace(frameLab, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheight = 0.96, `in` = Labtt)
        tkplace(tklabel(frameLab, text = "Choose the Point Name", 
            background = "white"), relx = 0.05, rely = 0.3, `in` = frameLab)
        objnames = rownames(MGvar$activedata)
        LabInput.val = tclVar(objnames[1])
        LabInput.ComboBox <- tkwidget(Labtt, "ComboBox", editable = FALSE, 
            values = objnames, width = 10, textvariable = LabInput.val)
        tkplace(LabInput.ComboBox, relx = 0.6, rely = 0.3, `in` = frameLab)
        On.Enter <- function() {
            Lab = as.character(tclvalue(LabInput.val))
            datlength = nrow(MGvar$MDSmat)
            indexlength = length(MGvar$indexLabeled)
            pointpos = 0
            for (i in 1:datlength) {
                if (rownames(MGvar$MDSmat)[i] == Lab) {
                  pointpos = i
                }
            }
            if (pointpos == 0) {
                tkmessageBox(title = "Point Label Error", message = "The name does not match any of the points. \nPlease check spelling", 
                  icon = "error", type = "ok")
            }
            if (pointpos > 0) {
                MGvar$indexLabeled[indexlength + 1] <<- pointpos
                names(MGvar$indexLabeled)[indexlength + 1] <<- rownames(MGvar$MDSmat)[pointpos]
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$indexLabeled.T1 <<- c(MGvar$indexLabeled.T1, 
                    MGvar$indexLabeled[indexlength + 1])
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$indexLabeled.T2 <<- c(MGvar$indexLabeled.T2, 
                    MGvar$indexLabeled[indexlength + 1])
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$indexLabeled.T3 <<- c(MGvar$indexLabeled.T3, 
                    MGvar$indexLabeled[indexlength + 1])
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$indexLabeled.T4 <<- c(MGvar$indexLabeled.T4, 
                    MGvar$indexLabeled[indexlength + 1])
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$indexLabeled.T5 <<- c(MGvar$indexLabeled.T5, 
                    MGvar$indexLabeled[indexlength + 1])
                }
            }
            tabplot()
            tkdestroy(Labtt)
        }
        tkplace(tkbutton(Labtt, text = "Enter", width = 15, command = function() On.Enter()), 
            relx = 0.3, rely = 0.7, `in` = frameLab)
        tkfocus(Labtt)
        tkbind(Labtt, "<Return>", On.Enter)
    }
    ChangePointCol <- function() {
        MGcomp$CPCol <<- tktoplevel()
        tkwm.resizable(MGcomp$CPCol, "0", "0")
        tkwm.title(MGcomp$CPCol, "Choose Point")
        tkwm.geometry(MGcomp$CPCol, "280x90")
        CPColcanvas = tkcanvas(MGcomp$CPCol, width = 300, height = 120, 
            bg = col.sec)
        tkplace(CPColcanvas, relx = 0, rely = 0, `in` = MGcomp$CPCol)
        frameCPC <- tkwidget(MGcomp$CPCol, "TitleFrame", text = "Active Cursor", 
            background = "white")
        tkplace(frameCPC, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = MGcomp$CPCol)
        tkplace(tklabel(frameCPC, text = "The cursor is now active. Please\nuse the cursor to select the point\nwhose colour you wish to change.", 
            background = "white"), relx = 0.1, rely = 0.2, `in` = frameCPC)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", GetCoordsLeftClick.PointCol)
            tkconfigure(img, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", GetCoordsLeftClick.PointCol)
            tkconfigure(img2, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", GetCoordsLeftClick.PointCol)
            tkconfigure(img3, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", GetCoordsLeftClick.PointCol)
            tkconfigure(img4, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", GetCoordsLeftClick.PointCol)
            tkconfigure(img5, cursor = "crosshair")
        }
        tkwait.window(MGcomp$CPCol)
        ActiveArrowCursor()
    }
    GetClosestPointChangeCol <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        Closest <- which.min(squared.Distance)
        PCol <- tktoplevel()
        tkwm.resizable(PCol, "0", "0")
        tkwm.title(PCol, "Point Colour")
        tkwm.geometry(PCol, "250x90")
        PColcanvas = tkcanvas(PCol, width = 300, height = 120, 
            bg = col.sec)
        tkplace(PColcanvas, relx = 0, rely = 0, `in` = PCol)
        framePC <- tkwidget(PCol, "TitleFrame", text = "Choose Point Colour", 
            background = "white")
        tkplace(framePC, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = PCol)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            ptcol <- MGvar$activeplot.pointcol.T1
            ptcol.temp <- MGvar$activeplot.pointcol.T1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            ptcol <- MGvar$activeplot.pointcol.T2
            ptcol.temp <- MGvar$activeplot.pointcol.T2
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            ptcol <- MGvar$activeplot.pointcol.T3
            ptcol.temp <- MGvar$activeplot.pointcol.T3
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            ptcol <- MGvar$activeplot.pointcol.T4
            ptcol.temp <- MGvar$activeplot.pointcol.T4
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            ptcol <- MGvar$activeplot.pointcol.T5
            ptcol.temp <- MGvar$activeplot.pointcol.T5
        }
        tkplace(tklabel(framePC, text = paste(rownames(MGvar$MDSmat)[Closest]), 
            background = "white"), relx = 0.1, rely = 0.35, `in` = framePC)
        ChangePtCol <- function() {
            ptcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = ptcol.temp, title = "Choose a Colour"))))
            if (nchar(ptcol.temp) > 0) 
                tkconfigure(Ptcol.but, bg = ptcol.temp)
            ptcol <<- ptcol.temp
            tkdestroy(PCol)
        }
        Ptcol.but <- tkbutton(PCol, text = "", width = 2, height = 1, 
            bg = ptcol, command = function() ChangePtCol())
        tkplace(Ptcol.but, relx = 0.8, rely = 0.35, `in` = framePC)
        tkwait.window(PCol)
        MGvar$MDSmat.Cols[Closest] <<- ptcol
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDSmat.Cols.T1[Closest] <<- ptcol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDSmat.Cols.T2[Closest] <<- ptcol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDSmat.Cols.T3[Closest] <<- ptcol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDSmat.Cols.T4[Closest] <<- ptcol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDSmat.Cols.T5[Closest] <<- ptcol
        }
        tabplot()
    }
    GetCoordsLeftClick.PointCol <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$MDSmat[, 1] - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$MDSmat[, 2] - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        GetClosestPointChangeCol(xClick, yClick, imgXcoords, 
            imgYcoords)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", OnPlotLeftClick)
        }
        tkdestroy(MGcomp$CPCol)
    }
    EnlargedActivePlot <- function() {
        ReusePOW()
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$EnActivePlot.switch.T1 <<- "on"
            MGvar$EnActivePlot.switch.T2 <<- "off"
            MGvar$EnActivePlot.switch.T3 <<- "off"
            MGvar$EnActivePlot.switch.T4 <<- "off"
            MGvar$EnActivePlot.switch.T5 <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$EnActivePlot.switch.T1 <<- "off"
            MGvar$EnActivePlot.switch.T2 <<- "on"
            MGvar$EnActivePlot.switch.T3 <<- "off"
            MGvar$EnActivePlot.switch.T4 <<- "off"
            MGvar$EnActivePlot.switch.T5 <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$EnActivePlot.switch.T1 <<- "off"
            MGvar$EnActivePlot.switch.T2 <<- "off"
            MGvar$EnActivePlot.switch.T3 <<- "on"
            MGvar$EnActivePlot.switch.T4 <<- "off"
            MGvar$EnActivePlot.switch.T5 <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$EnActivePlot.switch.T1 <<- "off"
            MGvar$EnActivePlot.switch.T2 <<- "off"
            MGvar$EnActivePlot.switch.T3 <<- "off"
            MGvar$EnActivePlot.switch.T4 <<- "on"
            MGvar$EnActivePlot.switch.T5 <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$EnActivePlot.switch.T1 <<- "off"
            MGvar$EnActivePlot.switch.T2 <<- "off"
            MGvar$EnActivePlot.switch.T3 <<- "off"
            MGvar$EnActivePlot.switch.T4 <<- "off"
            MGvar$EnActivePlot.switch.T5 <<- "on"
        }
        MGcomp$EActive <<- tktoplevel()
        tkwm.title(MGcomp$EActive, "Enlarged Active Plot")
        tkwm.geometry(MGcomp$EActive, "650x650")
        EActiveCanvas <- tkcanvas(MGcomp$EActive, width = 650, 
            height = 650, bg = col.sec)
        tkplace(EActiveCanvas, `in` = MGcomp$EActive)
        EC.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EActive)))
        EC.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EActive)))
        WidthScale = 650/1.5
        POhscale <- EC.width/WidthScale
        dimrat = 1.5/1.4
        POvscale <- POhscale/dimrat
        MGcomp$POimg <<- tkrplot(MGcomp$EActive, function() plotting2D(MGvar$MDSmat, 
            indexLabeled = MGvar$indexLabeled), hscale = POhscale, 
            vscale = POvscale)
        tkplace(MGcomp$POimg, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            `in` = MGcomp$EActive)
        tabplot()
        tkbind(MGcomp$POimg, "<Button-1>", OnPOPlotLeftClick)
        tkconfigure(MGcomp$POimg, cursor = "hand2")
        CopyAPToClip <- function() {
            tkrreplot(MGcomp$POimg)
        }
        tkplace(tkbutton(MGcomp$EActive, text = "Copy to Clipboard", 
            width = 20, command = function() CopyAPToClip()), 
            relx = 0.19, rely = 0.93, `in` = MGcomp$EActive)
        tkplace(tkbutton(MGcomp$EActive, text = "Plot Options", 
            width = 20, command = function() ConfPlotOptions()), 
            relx = 0.6, rely = 0.93, `in` = MGcomp$EActive)
        tkbind(MGcomp$EActive, "<Destroy>", function() {
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$EnActivePlot.switch.T1 <<- "off"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$EnActivePlot.switch.T2 <<- "off"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$EnActivePlot.switch.T3 <<- "off"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$EnActivePlot.switch.T4 <<- "off"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$EnActivePlot.switch.T5 <<- "off"
            }
        })
        MainPlotPOPMenu <- tkmenu(MGcomp$POimg, tearoff = FALSE)
        tkadd(MainPlotPOPMenu, "command", label = "Clear Added Point Labels", 
            command = ClearMainPOPPoints)
        tkadd(MainPlotPOPMenu, "command", label = "Label Specific Point", 
            command = LabelSpecificPoint)
        RightClickMainPOP <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", MGcomp$POimg))
            rooty <- as.integer(tkwinfo("rooty", MGcomp$POimg))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", MainPlotPOPMenu, xTxt, yTxt)
        }
        tkbind(MGcomp$POimg, "<Button-3>", RightClickMainPOP)
        tkbind(MGcomp$EActive, "<Configure>", resize.EConf)
    }
    resize.EConf <- function() {
        EC.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EActive)))
        EC.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EActive)))
        WidthScale = 650/1.5
        EShscale <- EC.width/WidthScale
        dimrat = 1.5/1.4
        ESvscale <- EShscale/dimrat
        tkrreplot(MGcomp$POimg, hscale = EShscale, vscale = ESvscale)
    }
    ClearMainPOPPoints <- function() {
        MGvar$indexLabeled <<- c()
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$indexLabeled.T1 <<- c()
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$indexLabeled.T2 <<- c()
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$indexLabeled.T3 <<- c()
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$indexLabeled.T4 <<- c()
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$indexLabeled.T5 <<- c()
        }
        tkrreplot(MGcomp$POimg)
        tabplot()
    }
    OnPOPlotLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", MGcomp$POimg)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", MGcomp$POimg)))
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$X.Coords - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$Y.Coords - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        msg <- paste("Label the point closest to these approximate plot coordinates: \n", 
            "x =", format(xPlotCoord, digits = 2), ",y =", format(yPlotCoord, 
                digits = 2), "?")
        mbval <- tkmessageBox(title = "Label Point Closest to These Approximate Plot Coordinates", 
            message = msg, type = "yesno", icon = "question")
        if (tclvalue(mbval) == "yes") 
            labelClosestPoint(xClick, yClick, imgXcoords, imgYcoords)
    }
    plotting3D <- function(data, title = MGvar$rglplot.title, 
        showtitle = MGvar$rglplot.title.show, xlabel = MGvar$rglplot.xlab, 
        ylabel = MGvar$rglplot.ylab, zlabel = MGvar$rglplot.zlab, 
        showpoints = MGvar$rglplot.showpoints, showlabels = MGvar$rglplot.showlabels, 
        colour = MGvar$rglplot.ptcol, ptsize = MGvar$rglplot.ptsize, 
        axmeas = MGvar$rglplot.axmeas) {
        if (showtitle == "yes") {
            graphtitle = title
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showpoints == "yes") {
            pointtype = "p"
        }
        if (showpoints == "no") {
            pointtype = "n"
        }
        if (axmeas == "yes") {
            AM = TRUE
        }
        if (axmeas == "no") {
            AM = FALSE
        }
        plot3d(data, aspect = TRUE, type = pointtype, col = colour, 
            size = ptsize, xlab = "", ylab = "", zlab = "", axes = AM)
        if (showlabels == "yes") {
            for (i in 1:nrow(data)) {
                text3d(data[i, ], texts = rownames(data)[i], 
                  col = MGvar$MDSmat.Cols[i], adj = c(0, 1), 
                  cex = 1)
            }
        }
        decorate3d(main = MGvar$sthreeDplot.title, xlab = xlabel, 
            ylab = ylabel, zlab = zlabel, axes = AM)
    }
    plotting3Dstatic <- function(data, title = MGvar$sthreeDplot.title, 
        Measure = MGvar$dMeas, showtitle = MGvar$sthreeDplot.title.show, 
        showmeas = MGvar$sthreeDplot.distmeas, xlabel = MGvar$sthreeDplot.xlab, 
        ylabel = MGvar$sthreeDplot.ylab, zlabel = MGvar$sthreeDplot.zlab, 
        showleg = MGvar$sthreeDplot.leg.show, bgcol = MGvar$sthreeDplot.bg, 
        showpoints = MGvar$sthreeDplot.showpoints, showlabs = MGvar$sthreeDplot.showlabels, 
        pointcex = MGvar$sthreeDplot.cex, pointshape = MGvar$sthreeDplot.type, 
        pointcol = MGvar$sthreeDplot.pointcol, showregX = MGvar$sthreeDplot.showregX, 
        regXcol = MGvar$sthreeDplot.regXcol, showregZ = MGvar$sthreeDplot.showregZ, 
        regZcol = MGvar$sthreeDplot.regZcol, showregY = MGvar$sthreeDplot.showregY, 
        regYcol = MGvar$sthreeDplot.regYcol, showaxes = MGvar$sthreeDplot.showaxes, 
        axcol = MGvar$sthreeDplot.axescol, showgrid = MGvar$sthreeDplot.showgrid, 
        Ang = MGvar$sthreeDplot.angle, HL = MGvar$sthreeDplot.HL, 
        yscale = MGvar$sthreeDplot.yscale) {
        if (showtitle == "yes") {
            graphtitle = title
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showmeas == "yes") {
            distanceMeasure = MGvar$dMeas
        }
        if (showmeas == "no") {
            distanceMeasure = ""
        }
        if (showpoints == "yes") {
            pointtype = "p"
        }
        if (showpoints == "no") {
            pointtype = "n"
        }
        if (showaxes == "yes") {
            tdbox = TRUE
        }
        if (showaxes == "no") {
            tdbox = FALSE
        }
        if (showgrid == "yes") {
            sgrid = TRUE
        }
        if (showgrid == "no") {
            sgrid = FALSE
        }
        if (HL == "yes") {
            HLight = TRUE
        }
        if (HL == "no") {
            HLight = FALSE
        }
        par(mar = c(3, 3, 3, 3))
        par(cex.axis = 0.8)
        params <- par(bg = bgcol)
        if (nrow(data) == 1 && ncol(data) == 1) {
            data = matrix(0, ncol = 2)
            scatterplot3d(data, type = "n")
            par(params)
        }
        else {
            if (!HLight) {
                statplot = scatterplot3d(data[, 1], data[, 2], 
                  data[, 3], type = pointtype, xlab = xlabel, 
                  ylab = ylabel, zlab = zlabel, main = graphtitle, 
                  pch = pointshape, if (!HLight) {
                    color = pointcol
                  }, cex.symbols = pointcex, box = tdbox, col.axis = axcol, 
                  grid = sgrid, angle = Ang, highlight.3d = HLight, 
                  scale.y = yscale)
            }
            if (HLight) {
                statplot = scatterplot3d(data[, 1], data[, 2], 
                  data[, 3], type = pointtype, xlab = xlabel, 
                  ylab = ylabel, zlab = zlabel, main = graphtitle, 
                  pch = pointshape, cex.symbols = pointcex, box = tdbox, 
                  col.axis = axcol, grid = sgrid, angle = Ang, 
                  highlight.3d = HLight, scale.y = yscale)
            }
            statplot.coords <- statplot$xyz.convert(data[, 1], 
                data[, 2], data[, 3])
            if (showlabs == "yes" && showpoints == "no") {
                for (i in 1:nrow(data)) {
                  text(statplot.coords$x[i], statplot.coords$y[i], 
                    labels = rownames(data)[i], cex = pointcex, 
                    col = MGvar$MDSmat.Cols[i])
                }
            }
            if (showlabs == "yes" && showpoints == "yes") {
                for (i in 1:nrow(data)) {
                  text(statplot.coords$x[i], statplot.coords$y[i], 
                    labels = rownames(data)[i], cex = 0.6 * pointcex, 
                    col = MGvar$MDSmat.Cols[i])
                }
            }
            par(params)
            if (showmeas == "yes") {
                mtext(text = Measure, side = 3, line = 2.8, adj = 0, 
                  cex = 0.7, col = "black")
            }
            if (showregX == "yes") {
            }
            if (showregY == "yes") {
            }
            if (showregZ == "yes") {
            }
        }
    }
    EnlargedStat3DPlot <- function() {
        ReusePOW()
        MGvar$EStat.switch <<- "on"
        MGcomp$EStat <<- tktoplevel()
        tkwm.title(MGcomp$EStat, "Enlarged Static 3D Plot")
        tkwm.geometry(MGcomp$EStat, "650x650")
        EStatCanvas <- tkcanvas(MGcomp$EStat, width = 650, height = 650, 
            bg = col.sec)
        tkplace(EStatCanvas, `in` = MGcomp$EStat)
        POhscale <- 1.5
        POvscale <- 1.4
        MGcomp$POstatimg <<- tkrplot(MGcomp$EStat, function() plotting3Dstatic(MGvar$MDSmat), 
            hscale = POhscale, vscale = POvscale)
        tkplace(MGcomp$POstatimg, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            `in` = MGcomp$EStat)
        CopyToClip <- function() {
            tkrreplot(MGcomp$POstatimg)
        }
        tkplace(tkbutton(MGcomp$EStat, text = "Copy to Clipboard", 
            width = 20, command = function() tkrreplot(MGcomp$POstatimg)), 
            relx = 0.19, rely = 0.93, `in` = MGcomp$EStat)
        tkplace(tkbutton(MGcomp$EStat, text = "Plot Options", 
            width = 20, command = function() threeDPlotOptions()), 
            relx = 0.6, rely = 0.93, `in` = MGcomp$EStat)
        tkbind(MGcomp$EStat, "<Configure>", resize.EStat)
    }
    resize.EStat <- function() {
        ESta.height <- as.numeric(tclvalue(tkwinfo("height", 
            MGcomp$EStat)))
        ESta.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EStat)))
        WidthScale = 650/1.5
        EShscale <- ESta.width/WidthScale
        dimrat = 1.5/1.4
        ESvscale <- EShscale/dimrat
        tkrreplot(MGcomp$POstatimg, hscale = EShscale, vscale = ESvscale)
    }
    plotShepard <- function(shepdists, conf, showtitle = MGvar$shepplot.title.show, 
        showlabs = MGvar$shepplot.labs.show, showleg = MGvar$shepplot.leg.show, 
        bgcol = MGvar$shepplot.bg, showpoints = MGvar$shepplot.showpoints, 
        pointcex = MGvar$shepplot.cex, pointshape = MGvar$shepplot.type, 
        pointcol = MGvar$shepplot.pointcol, showcurve = MGvar$shepplot.curve.show, 
        curvetype = MGvar$shepplot.curve.type, curvecol = MGvar$shepplot.curvecol, 
        xmeas = MGvar$shepplot.Axes.xaxt, bold = MGvar$shepplot.bold, 
        ymeas = MGvar$shepplot.Axes.yaxt, axcol = MGvar$shepplot.Axescol, 
        Tab = tclvalue(MGvar$ActivePlottingTab), showplabs = MGvar$shepplot.showlabels) {
        if (showtitle == "yes") {
            MGvar$activeshepplot.title <<- paste("Shepard Plot for", 
                Tab)
            graphtitle <- MGvar$activeshepplot.title
        }
        if (showtitle == "no") {
            graphtitle <- ""
        }
        if (showpoints == "yes") {
            pointtype = "p"
        }
        if (showpoints == "no") {
            pointtype = "n"
        }
        par(mar = c(3, 3, 2, 5))
        params <- par(bg = bgcol)
        if (nrow(conf) == 1 && ncol(conf) == 1) {
            dumx = c(1, 1)
            plot(dumx, type = "n", asp = 1, xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n")
            par(params)
        }
        else {
            MGvar$Shepx <<- as.vector(0)
            MGvar$Shepy <<- as.vector(0)
            ShepComp <- Shepard(as.dist(shepdists), conf)
            dmat <- matrix(0, nrow(shepdists), ncol(shepdists))
            for (i in 1:nrow(shepdists)) {
                for (j in 1:nrow(shepdists)) {
                  comp = 0
                  for (k in 1:ncol(conf)) {
                    comp = comp + (conf[i, k] - conf[j, k])^2
                  }
                  dmat[i, j] = sqrt(comp)
                }
            }
            tShepy <- as.vector(0)
            len = nrow(shepdists)
            count = 1
            for (i in 1:(len - 1)) {
                for (j in (i + 1):len) {
                  tShepy[count] <- dmat[j, i]
                  names(tShepy)[count] <- paste(rownames(shepdists)[i], 
                    ":", colnames(shepdists)[j])
                  count = count + 1
                }
            }
            count = 1
            for (i in 1:(len - 1)) {
                for (j in (i + 1):len) {
                  MGvar$tShepx[count] <<- shepdists[j, i]
                  names(MGvar$tShepx)[count] <<- paste(rownames(shepdists)[i], 
                    ":", colnames(shepdists)[j])
                  count = count + 1
                }
            }
            Shepvals <- data.frame(tShepx = MGvar$tShepx, tShepy, 
                seq(1:length(MGvar$tShepx)))
            rownames(Shepvals) = names(MGvar$tShepx)
            OrderedVals = Shepvals[order(Shepvals$tShepx), ]
            MGvar$OrderedVals <<- as.matrix(OrderedVals)
            MGvar$Shepx <<- MGvar$OrderedVals[, 1]
            MGvar$Shepy <<- MGvar$OrderedVals[, 2]
            names(MGvar$Shepx) <<- rownames(MGvar$OrderedVals)
            if (MGvar$shep.firstrun == "yes") {
                ShepCompfunc = ShepFirstRun(ShepComp, tShepx = MGvar$tShepx, 
                  Shepx = MGvar$Shepx)
                MGvar$ShepPointindex <<- ShepCompfunc$ShepPointindex
                MGvar$DistFunc <<- ShepCompfunc$f.DistFunc
                ActiveDistFunc()
                MGvar$shep.firstrun <<- "no"
            }
            Ymin <- min(MGvar$Shepy)
            Ymax <- max(MGvar$Shepy)
            Xmin <- min(MGvar$Shepx)
            Xmax <- max(MGvar$Shepx)
            Gmax <- max(Ymax, Xmax)
            Xrange <- Xmax - Xmin
            Yrange <- Ymax - Ymin
            if (showleg == "yes") {
                Ymax = Ymax + 0.15 * Yrange
            }
            plot(MGvar$Shepx, MGvar$Shepy, type = pointtype, 
                ylim = c(0, Gmax), xlim = c(0, Gmax), cex = pointcex, 
                xlab = "", ylab = "", main = graphtitle, pch = pointshape, 
                col = pointcol, xaxt = xmeas, yaxt = ymeas, fg = axcol, 
                cex.axis = 0.7, mex = 1, mgp = c(2.5, 0.5, 0))
            if (showplabs == "yes") {
                text(MGvar$Shepx, ShepComp$y, names(MGvar$Shepx), 
                  cex = pointcex, col = pointcol)
            }
            if (showcurve == "yes") {
                if (bold) {
                  lnwid = 2
                }
                if (!bold) {
                  lnwid = 1
                }
                if (!MGvar$is.Metric) {
                  lines(ShepComp$x, ShepComp$yf, type = "S", 
                    col = curvecol, lty = curvetype, lwd = lnwid)
                }
                if (MGvar$is.Metric) {
                  lines(ShepComp$x, ShepComp$x, col = curvecol, 
                    lty = curvetype, lwd = lnwid)
                }
            }
            par(params)
            if (showlabs == "yes") {
                mtext(text = "Observed Dissimilarity", side = 1, 
                  line = 1.2, adj = 0.5, cex = 1, col = "black")
                mtext(text = "Ordination Distance", side = 2, 
                  line = 2.1, adj = 0.5, cex = 1, col = "black")
            }
            if (showleg == "yes") {
                legend(Xmin, 1.06 * Ymax, " Line Label", bty = "n", 
                  cex = 0.8, pt.bg = "white", lty = curvetype, 
                  col = curvecol)
                legend(Xmin + 0.02 * Xrange, Ymax, "    Point Label", 
                  bty = "n", cex = 0.8, pch = pointshape, pt.bg = "white", 
                  col = pointcol)
            }
            MGvar$Shep.parPlotSize <<- par("plt")
            MGvar$Shep.usrCoords <<- par("usr")
            MGvar$Shep.labelsVec <<- rownames(MGvar$OrderedVals)
            if (length(MGvar$Shep.indexLabeled) > 0) {
                clrs = rainbow(length(MGvar$Shep.indexLabeled))
                points(MGvar$Shepx, MGvar$Shepy, type = pointtype, 
                  ylim = c(0, Ymax), xlim = c(0, Xmax), cex = pointcex, 
                  xlab = "", ylab = "", main = graphtitle, pch = pointshape, 
                  col = "grey", xaxt = xmeas, yaxt = ymeas, fg = axcol, 
                  cex.axis = 0.7, mex = 1, mgp = c(2.5, 0.5, 
                    0))
                if (showcurve == "yes") {
                  if (bold) {
                    lnwid = 2
                  }
                  if (!bold) {
                    lnwid = 1
                  }
                  if (!MGvar$is.Metric) {
                    lines(ShepComp$x, ShepComp$yf, type = "S", 
                      col = curvecol, lty = curvetype, lwd = lnwid)
                  }
                  if (MGvar$is.Metric) {
                    lines(ShepComp$x, ShepComp$x, col = curvecol, 
                      lty = curvetype, lwd = lnwid)
                  }
                }
                for (i in (1:length(MGvar$Shep.indexLabeled))) {
                  indexClosest <- MGvar$Shep.indexLabeled[i]
                  if (length(MGvar$Shep.indexLabeled) < 10) {
                    text(MGvar$Shepx[indexClosest], MGvar$Shepy[indexClosest], 
                      labels = MGvar$Shep.labelsVec[indexClosest], 
                      cex = 0.6, pos = MGvar$GenSet.ILPos)
                  }
                  points(MGvar$Shepx[indexClosest], MGvar$Shepy[indexClosest], 
                    cex = 2 * pointcex, col = clrs[i], pch = 16)
                }
            }
            segments(MGvar$first.xPlotCoord, MGvar$first.yPlotCoord, 
                MGvar$first.xPlotCoord, MGvar$latest.yPlotCoord, 
                col = "blue")
            segments(MGvar$first.xPlotCoord, MGvar$first.yPlotCoord, 
                MGvar$latest.xPlotCoord, MGvar$first.yPlotCoord, 
                col = "blue")
            segments(MGvar$latest.xPlotCoord, MGvar$first.yPlotCoord, 
                MGvar$latest.xPlotCoord, MGvar$latest.yPlotCoord, 
                col = "blue")
            segments(MGvar$first.xPlotCoord, MGvar$latest.yPlotCoord, 
                MGvar$latest.xPlotCoord, MGvar$latest.yPlotCoord, 
                col = "blue")
        }
    }
    Internal.ShepFirstRun <- function(ShepComp) {
        if (MGvar$shep.firstrun == "yes") {
            time1 = proc.time()
            MGvar$ShepPointindex <<- matrix(nrow = length(MGvar$tShepx), 
                ncol = 2)
            rownames(MGvar$ShepPointindex) <<- names(MGvar$tShepx)
            MGvar$ShepPointindex[, 2] <<- seq(1:length(MGvar$tShepx))
            for (i in 1:length(MGvar$tShepx)) {
                for (j in 1:length(MGvar$Shepx)) {
                  if (names(MGvar$tShepx)[i] == names(MGvar$Shepx)[j]) {
                    MGvar$ShepPointindex[i, 1] <<- j
                    MGvar$DistFunc[j] <<- ShepComp$yf[i]
                    names(MGvar$DistFunc)[j] <<- names(MGvar$Shepx)[j]
                  }
                }
            }
            MGvar$shep.firstrun <<- "no"
            ActiveDistFunc()
            time2 = proc.time()
            print(time2 - time1)
        }
    }
    Shep.labelClosestPoint <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance.Sh <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        indexClosest <- which.min(squared.Distance.Sh)
        MGvar$Shep.indexLabeled <<- c(MGvar$Shep.indexLabeled, 
            indexClosest)
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
        tabplot()
        if (MGvar$EnShep.switch == "on") {
            tkrreplot(MGcomp$imgEShep)
        }
        tclvalue(Shep.WhichShow.var) <<- "Active Plot"
        tkconfigure(Shep.WhichShow.ComboBox, textvariable = Shep.WhichShow.var)
    }
    Shep.OnLeftClick <- function(x, y) {
        if (MGvar$GenSet.CalcShep == "yes") {
            xClick <- x
            yClick <- y
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                imgshep)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                imgshep)))
            xMin <- MGvar$Shep.parPlotSize[1] * width
            xMax <- MGvar$Shep.parPlotSize[2] * width
            yMin <- MGvar$Shep.parPlotSize[3] * height
            yMax <- MGvar$Shep.parPlotSize[4] * height
            rangeX <- MGvar$Shep.usrCoords[2] - MGvar$Shep.usrCoords[1]
            rangeY <- MGvar$Shep.usrCoords[4] - MGvar$Shep.usrCoords[3]
            imgXcoords <- (MGvar$Shepx - MGvar$Shep.usrCoords[1]) * 
                (xMax - xMin)/rangeX + xMin
            imgYcoords <- (MGvar$Shepy - MGvar$Shep.usrCoords[3]) * 
                (yMax - yMin)/rangeY + yMin
            xClick <- as.numeric(xClick) + 0.5
            yClick <- as.numeric(yClick) + 0.5
            yClick <- height - yClick
            xPlotCoord <- MGvar$Shep.usrCoords[1] + (xClick - 
                xMin) * rangeX/(xMax - xMin)
            yPlotCoord <- MGvar$Shep.usrCoords[3] + (yClick - 
                yMin) * rangeY/(yMax - yMin)
            Shep.labelClosestPoint(xClick, yClick, imgXcoords, 
                imgYcoords)
        }
    }
    Shep.LabelSpecificPoint <- function() {
        LabStt = tktoplevel()
        tkwm.resizable(LabStt, "0", "0")
        tkwm.deiconify(LabStt)
        tkwm.title(LabStt, "Shepard Point Label")
        tkwm.geometry(LabStt, "300x210")
        LabScanvas = tkcanvas(LabStt, width = 300, height = 250, 
            bg = col.sec)
        tkplace(LabScanvas, `in` = LabStt)
        frameLabS <- tkwidget(LabStt, "TitleFrame", text = "Label a Specific Point", 
            background = "white")
        tkplace(frameLabS, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = LabStt)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frameLabS, text = "Each point on the Shepard Plot represents the MDS \nfound distance between each object pairing in the \nactive data. In order to choose a point please select \nthe names of two objects in the data.", 
            font = fontsmall), relx = 0.05, rely = 0.1, `in` = frameLabS)
        tkplace(tklabel(frameLabS, text = "Enter the Object1", 
            background = "white"), relx = 0.08, rely = 0.5, `in` = frameLabS)
        objnames = rownames(MGvar$MDSmat)
        LabInput1.val = tclVar(objnames[1])
        LabInput1.ComboBox <- tkwidget(LabStt, "ComboBox", editable = FALSE, 
            values = objnames, width = 12, textvariable = LabInput1.val)
        tkplace(LabInput1.ComboBox, relx = 0.6, rely = 0.5, `in` = frameLabS)
        tkplace(tklabel(frameLabS, text = "Enter the Object2", 
            background = "white"), relx = 0.08, rely = 0.65, 
            `in` = frameLabS)
        LabInput2.val = tclVar(objnames[2])
        LabInput2.ComboBox <- tkwidget(LabStt, "ComboBox", editable = FALSE, 
            values = objnames, width = 12, textvariable = LabInput2.val)
        tkplace(LabInput2.ComboBox, relx = 0.6, rely = 0.65, 
            `in` = frameLabS)
        On.SEnter <- function() {
            Lab1 = as.character(tclvalue(LabInput1.val))
            Lab2 = as.character(tclvalue(LabInput2.val))
            datlength = nrow(MGvar$MDSmat)
            indexlength = length(MGvar$Shep.indexLabeled)
            pointpos1 = 0
            for (i in 1:datlength) {
                if (rownames(MGvar$MDSmat)[i] == Lab1) {
                  pointpos1 = i
                }
            }
            pointpos2 = 0
            for (j in 1:datlength) {
                if (rownames(MGvar$MDSmat)[j] == Lab2) {
                  pointpos2 = j
                }
            }
            if (pointpos2 < pointpos1) {
                dummy = pointpos1
                pointpos1 = pointpos2
                pointpos2 = dummy
            }
            if (pointpos1 == 0 || pointpos2 == 0) {
                tkmessageBox(title = "Point Label Error", message = "At least one of the names does not match any of the points. \nPlease check spelling", 
                  icon = "error", type = "ok")
            }
            if (pointpos1 > 0 && pointpos2 > 0) {
                indexmat = matrix(nrow = datlength, ncol = datlength)
                count = 1
                for (i in 1:(datlength - 1)) {
                  for (j in (i + 1):datlength) {
                    indexmat[i, j] = count
                    count = count + 1
                  }
                }
                pointname <- names(MGvar$tShepx)[indexmat[pointpos1, 
                  pointpos2]]
                for (i in 1:length(MGvar$Shepx)) {
                  if (pointname == names(MGvar$Shepx)[i]) {
                    MGvar$Shep.indexLabeled[indexlength + 1] <<- i
                    names(MGvar$Shep.indexLabeled)[indexlength + 
                      1] <<- pointname
                  }
                }
            }
            tkdestroy(LabStt)
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            tabplot()
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEShep)
            }
            tclvalue(Shep.WhichShow.var) <<- "Active Plot"
            tkconfigure(Shep.WhichShow.ComboBox, textvariable = Shep.WhichShow.var)
        }
        tkplace(tkbutton(LabStt, text = "Enter", width = 15, 
            command = function() On.SEnter()), relx = 0.3, rely = 0.82, 
            `in` = frameLabS)
        tkfocus(LabStt)
        tkbind(LabStt, "<Return>", On.SEnter)
    }
    plotScree <- function(stress, current, best, Cdim, Odim, 
        showtitle = MGvar$screeplot.title.show, showlabs = MGvar$screeplot.labs.show, 
        showleg = MGvar$screeplot.leg.show, bgcol = MGvar$screeplot.bg, 
        pointtype = MGvar$screeplot.points.show, showCdim = MGvar$screeplot.Cdim.show, 
        showOdim = MGvar$screeplot.Odim.show, CPcol = MGvar$screeplot.Ccol, 
        OPcol = MGvar$screeplot.Ocol, showcurve = MGvar$screeplot.curve.show, 
        curvetype = MGvar$screeplot.curve.type, curvecol = MGvar$screeplot.curvecol, 
        Cline = MGvar$screeplot.Cline.show, Oline = MGvar$screeplot.Oline.show, 
        xmeas = MGvar$screeplot.Axes.xaxt, ymeas = MGvar$screeplot.Axes.yaxt, 
        axcol = MGvar$screeplot.Axescol, Tab = tclvalue(MGvar$ActivePlottingTab)) {
        if (showtitle == "yes") {
            MGvar$activescreeplot.title <<- paste("Scree Plot for", 
                Tab)
            if (Tab == "Tab1" && MGvar$MDStype.T1 == "ClasScal") {
                MGvar$activescreeplot.title <<- paste("Eigen-Plot for", 
                  Tab)
            }
            if (Tab == "Tab2" && MGvar$MDStype.T2 == "ClasScal") {
                MGvar$activescreeplot.title <<- paste("Eigen-Plot for", 
                  Tab)
            }
            if (Tab == "Tab3" && MGvar$MDStype.T3 == "ClasScal") {
                MGvar$activescreeplot.title <<- paste("Eigen-Plot for", 
                  Tab)
            }
            if (Tab == "Tab4" && MGvar$MDStype.T4 == "ClasScal") {
                MGvar$activescreeplot.title <<- paste("Eigen-Plot for", 
                  Tab)
            }
            if (Tab == "Tab5" && MGvar$MDStype.T5 == "ClasScal") {
                MGvar$activescreeplot.title <<- paste("Eigen-Plot for", 
                  Tab)
            }
            graphtitle = MGvar$activescreeplot.title
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showcurve == "yes") {
            curve = "l"
        }
        if (showcurve == "no") {
            curve = "n"
        }
        par(mar = c(4, 3, 2, 2))
        params <- par(bg = bgcol)
        if (length(stress) == 1) {
            plot(stress, type = "n", asp = 1, xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n")
            par(params)
        }
        else {
            cap = min((MGvar$maxdims - 2), 10)
            plot(stress[1:cap], type = "n", main = graphtitle, 
                col = curvecol, xaxt = xmeas, yaxt = ymeas, fg = axcol, 
                sub = "", ylab = "", xlab = "")
            if (showcurve == "yes") {
                lines(stress[1:2], lty = 2)
                linseq = seq(2, cap)
                lines(linseq, stress[2:cap], lty = curvetype)
            }
            par(params)
            if (pointtype == "yes") {
                points(stress[1:cap])
            }
            if (showCdim == "yes") {
                points(current, cex = 2, pch = 16, col = CPcol)
            }
            if (showOdim == "yes") {
                points(best, cex = 1.6, pch = 16, col = OPcol)
            }
            if (Cline == "yes") {
                abline(v = Cdim, lty = 3, col = CPcol)
            }
            if (Oline == "yes") {
                abline(v = Odim + 0.05, lty = 3, col = OPcol)
            }
            if (showleg == "yes") {
                if (showCdim == "yes") {
                  legend(7, (0.99 * max(stress)), "Current", 
                    pch = 16, bty = "n", pt.bg = "white", lty = 2, 
                    col = CPcol)
                }
                else {
                  legend(7, (0.99 * max(stress)), "Current", 
                    bty = "n", pt.bg = "white", lty = 2, col = CPcol)
                }
                if (showOdim == "yes") {
                  legend(7, (0.91 * max(stress)), "Optimal", 
                    pch = 16, bty = "n", pt.bg = "white", lty = 2, 
                    col = OPcol)
                }
                else {
                  legend(7, (0.91 * max(stress)), "Optimal", 
                    bty = "n", pt.bg = "white", lty = 2, col = OPcol)
                }
            }
            if (showlabs == "yes") {
                CS = "no"
                mtext(text = "Dimensions", side = 1, line = 2, 
                  adj = 0.5, cex = 1, col = "black")
                if (Tab == "Tab1" && MGvar$MDStype.T1 == "ClasScal") {
                  mtext(text = "Eigen Values", side = 2, line = 2.1, 
                    adj = 0.5, cex = 1, col = "black")
                  CS = "yes"
                }
                if (Tab == "Tab2" && MGvar$MDStype.T2 == "ClasScal") {
                  mtext(text = "Eigen Values", side = 2, line = 2.1, 
                    adj = 0.5, cex = 1, col = "black")
                  CS = "yes"
                }
                if (Tab == "Tab3" && MGvar$MDStype.T3 == "ClasScal") {
                  mtext(text = "Eigen Values", side = 2, line = 2.1, 
                    adj = 0.5, cex = 1, col = "black")
                  CS = "yes"
                }
                if (Tab == "Tab4" && MGvar$MDStype.T4 == "ClasScal") {
                  mtext(text = "Eigen Values", side = 2, line = 2.1, 
                    adj = 0.5, cex = 1, col = "black")
                  CS = "yes"
                }
                if (Tab == "Tab5" && MGvar$MDStype.T5 == "ClasScal") {
                  mtext(text = "Eigen Values", side = 2, line = 2.1, 
                    adj = 0.5, cex = 1, col = "black")
                  CS = "yes"
                }
                if (CS == "no") {
                  mtext(text = "Stress", side = 2, line = 2.1, 
                    adj = 0.5, cex = 1, col = "black")
                }
            }
        }
    }
    StressPlot <- function(showtitle = MGvar$stressplot.title.show, 
        showtime = MGvar$stressplot.time.show, showlabs = MGvar$stressplot.labs.show, 
        bgcol = MGvar$stressplot.bg, curvetype = MGvar$stressplot.curve.type, 
        curvecol = MGvar$stressplot.curvecol, xmeas = MGvar$stressplot.Axes.xaxt, 
        ymeas = MGvar$stressplot.Axes.yaxt, axcol = MGvar$stressplot.Axescol) {
        if (showtitle == "yes") {
            graphtitle = "Stress Plot"
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        par(mar = c(4, 3, 2, 2))
        params <- par(bg = bgcol)
        if (length(MGvar$stressitervec) == 1) {
            plot(MGvar$stressitervec, type = "n", col = curvecol, 
                xaxt = "n", yaxt = "n", ylab = "", xlab = "")
            par(params)
        }
        else {
            len = length(MGvar$stressitervec)
            plotseq = cbind(seq(1:(len)), MGvar$stressitervec)
            if (len > 50) {
                plotseq = cbind(plotseq[(len - 49):len, 1], MGvar$stressitervec[(len - 
                  49):len])
            }
            tline = max(plotseq[, 2])
            tline = tline + 0.4 * (max(plotseq[, 2]) - min(plotseq[, 
                2]))
            bline = max(0, min(plotseq[, 2]) - 0.6 * (max(plotseq[, 
                2]) - min(plotseq[, 2])))
            plot(plotseq, type = "l", main = graphtitle, ylim = c(bline, 
                tline), lty = curvetype, col = curvecol, xaxt = xmeas, 
                yaxt = ymeas, fg = axcol, sub = "", ylab = "", 
                xlab = "")
            if (showtime == "yes") {
                timeval <- paste("Time: ", round(MGvar$proctime, 
                  1), " seconds")
                stressval <- paste("Stress: ", round(MGvar$stressitervec[len], 
                  4))
                iterval <- paste("Iterations: ", len)
                legend("topright", c(timeval, stressval, iterval), 
                  cex = 0.7)
            }
            if (showlabs == "yes") {
                mtext(text = "Iterations", side = 1, line = 2, 
                  adj = 0.5, cex = 1, col = "black")
                mtext(text = "Stress", side = 2, line = 2.1, 
                  adj = 0.5, cex = 1, col = "black")
            }
            par(params)
        }
    }
    StressPlot2 <- function(showtitle = MGvar$stressplot.title.show, 
        showtime = MGvar$stressplot.time.show, showlabs = MGvar$stressplot.labs.show, 
        bgcol = MGvar$stressplot.bg, curvetype = MGvar$stressplot.curve.type, 
        curvecol = MGvar$stressplot.curvecol, xmeas = MGvar$stressplot.Axes.xaxt, 
        ymeas = MGvar$stressplot.Axes.yaxt, axcol = MGvar$stressplot.Axescol) {
        if (showtitle == "yes") {
            graphtitle = "Stress Plot (Logged Differences)"
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        par(mar = c(4, 3, 2, 2))
        params <- par(bg = bgcol)
        if (length(MGvar$stressitervec) == 1) {
            plot(MGvar$stressitervec, type = "n", col = curvecol, 
                xaxt = "n", yaxt = "n", ylab = "", xlab = "")
            par(params)
        }
        else {
            len = length(MGvar$stressitervec)
            plotseq = cbind(seq(1:(len)), MGvar$stressitervec)
            if (len > 50) {
                plotseq = cbind(plotseq[(len - 49):len, 1], MGvar$stressitervec[(len - 
                  49):len])
            }
            MGvar$diffvec <<- plotseq[, 2]
            MGvar$diffvec <<- log(-(MGvar$diffvec[-1] - MGvar$diffvec[-length(MGvar$diffvec)]))
            tline = max(MGvar$diffvec)
            tline = tline + 0.4 * (max(MGvar$diffvec) - min(MGvar$diffvec))
            bline = min(0, min(MGvar$diffvec) - 0.6 * (max(MGvar$diffvec) - 
                min(MGvar$diffvec)))
            plot(plotseq[-1, 1], MGvar$diffvec, type = "l", main = graphtitle, 
                ylim = c(bline, tline), lty = curvetype, col = curvecol, 
                xaxt = xmeas, yaxt = ymeas, fg = axcol, sub = "", 
                ylab = "", xlab = "")
            if (showtime == "yes") {
                timeval <- paste("Time: ", round(MGvar$proctime, 
                  1), " seconds")
                stressval <- paste("Stress: ", round(MGvar$stressitervec[len], 
                  4))
                iterval <- paste("Iterations: ", len)
                legend("topright", c(timeval, stressval, iterval), 
                  cex = 0.7)
            }
            if (showlabs == "yes") {
                mtext(text = "Iterations", side = 1, line = 2, 
                  adj = 0.5, cex = 1, col = "black")
                mtext(text = "Logged Stress Differences", side = 2, 
                  line = 2.1, adj = 0.5, cex = 1, col = "black")
            }
            par(params)
        }
    }
    EnlargedStress <- function() {
        ReusePOW()
        MGvar$EnStress.switch <<- "on"
        MGcomp$EStress <<- tktoplevel()
        tkwm.geometry(MGcomp$EStress, "773x575")
        tkwm.title(MGcomp$EStress, "Enlarged Stress Plot")
        EStresscanvas = tkcanvas(MGcomp$EStress, width = 818, 
            height = 580, bg = col.sec)
        tkplace(EStresscanvas, `in` = MGcomp$EStress)
        ES.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EStress)))
        ES.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EStress)))
        WidthScale = 773/1.8
        EShscale <- ES.width/WidthScale
        dimrat = 1.8/1.3
        ESvscale <- EShscale/dimrat
        MGcomp$imgEStress <<- tkrplot(MGcomp$EStress, function() StressPlot(), 
            hscale = EShscale, vscale = ESvscale)
        tkplace(MGcomp$imgEStress, relx = 0.05, rely = 0.02, 
            relwidth = 0.9, `in` = MGcomp$EStress)
        CopyStToClip <- function() {
            tkrreplot(MGcomp$imgEStress)
        }
        tkplace(tkbutton(MGcomp$EStress, text = "Copy to Clipboard", 
            width = 20, command = function() CopyStToClip()), 
            relx = 0.22, rely = 0.93, `in` = MGcomp$EStress)
        tkplace(tkbutton(MGcomp$EStress, text = "Plot Options", 
            width = 20, command = function() StressPlotOps()), 
            relx = 0.6, rely = 0.93, `in` = MGcomp$EStress)
        tkbind(MGcomp$imgEStress, "<Destroy>", function() {
            MGvar$EnStress.switch <<- "off"
        })
        tkbind(MGcomp$EStress, "<Configure>", resize.EStress)
    }
    resize.EStress <- function() {
        ES.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EStress)))
        ES.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EStress)))
        WidthScale = 773/1.8
        MGvar$ESt1hscale <<- ES.width/WidthScale
        dimrat = 1.8/1.3
        MGvar$ESt1vscale <<- MGvar$ESt1hscale/dimrat
        tkrreplot(MGcomp$imgEStress, hscale = MGvar$ESt1hscale, 
            vscale = MGvar$ESt1vscale)
    }
    Replot.imgEStress <- function() {
        resize.EStress()
        resize.EStress2()
    }
    EnlargedStress2 <- function() {
        ReusePOW()
        MGvar$EnStress2.switch <<- "on"
        MGcomp$EStress2 <<- tktoplevel()
        tkwm.geometry(MGcomp$EStress2, "773x575")
        tkwm.title(MGcomp$EStress2, "Enlarged Stress Plot 2")
        EStress2canvas = tkcanvas(MGcomp$EStress2, width = 818, 
            height = 580, bg = col.sec)
        tkplace(EStress2canvas, `in` = MGcomp$EStress2)
        EShscale = 1.8
        ESvscale = 1.3
        MGcomp$imgEStress2 <<- tkrplot(MGcomp$EStress2, function() StressPlot2(), 
            hscale = EShscale, vscale = ESvscale)
        tkplace(MGcomp$imgEStress2, relx = 0.05, rely = 0.02, 
            relwidth = 0.9, `in` = MGcomp$EStress2)
        CopyStToClip <- function() {
            tkrreplot(MGcomp$imgEStress2)
        }
        tkplace(tkbutton(MGcomp$EStress2, text = "Copy to Clipboard", 
            width = 20, command = function() CopyStToClip()), 
            relx = 0.22, rely = 0.93, `in` = MGcomp$EStress2)
        tkplace(tkbutton(MGcomp$EStress2, text = "Plot Options", 
            width = 20, command = function() StressPlotOps()), 
            relx = 0.6, rely = 0.93, `in` = MGcomp$EStress2)
        tkbind(MGcomp$imgEStress2, "<Destroy>", function() {
            MGvar$EnStress2.switch <<- "off"
        })
        tkbind(MGcomp$EStress2, "<Configure>", resize.EStress2)
    }
    resize.EStress2 <- function() {
        ES2.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EStress2)))
        ES2.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EStress2)))
        WidthScale = 773/1.8
        MGvar$ESt2hscale <<- ES2.width/WidthScale
        dimrat = 1.8/1.3
        MGvar$ESt2vscale <<- MGvar$ESt2hscale/dimrat
        tkrreplot(MGcomp$imgEStress2, hscale = MGvar$ESt2hscale, 
            vscale = MGvar$ESt2vscale)
    }
    Procrust <- function(X, Y, labels = NULL, col1 = "sea green", 
        col2 = "blue", ...) {
        X <- scale(X, scale = F)
        Y <- scale(Y, scale = F)
        q <- ncol(X)
        p <- ncol(Y)
        while (q < p) {
            X <- cbind(X, 0)
            q <- q + 1
        }
        Cmat <- t(Y) %*% X
        swd <- svd(Cmat)
        Amat <- swd$v %*% t(swd$u)
        rho <- sum(eigen(Amat %*% Cmat)$values)/sum(eigen(t(X) %*% 
            X)$values)
        Xtilde <- rho * X %*% Amat
    }
    ProcrustesMDSGUIplot <- function(conf1, conf2, showtitle = MGvar$procplot.title.show, 
        showleg = MGvar$procplot.leg.show, ylabel = MGvar$procplot.ylab, 
        xlabel = MGvar$procplot.xlab, bgcol = MGvar$procplot.bg, 
        showpoints1 = MGvar$procplot.showpoints1, showlabs1 = MGvar$procplot.labs1, 
        pointcex1 = MGvar$procplot.cex1, pointshape1 = MGvar$procplot.type1, 
        pointcol1 = MGvar$procplot.point1col, showpoints2 = MGvar$procplot.showpoints2, 
        showlabs2 = MGvar$procplot.labs2, pointcex2 = MGvar$procplot.cex2, 
        pointshape2 = MGvar$procplot.type2, pointcol2 = MGvar$procplot.point2col, 
        ymeas = MGvar$procplot.yaxt, xmeas = MGvar$procplot.xaxt, 
        axcol = MGvar$procplot.axescol, tab1 = "Plot1", tab2 = "Plot2", 
        showreg1 = MGvar$procplot.showreg1, showreg2 = MGvar$procplot.showreg2, 
        regcol1 = MGvar$procplot.regcol1, regcol2 = MGvar$procplot.regcol2, 
        Mvup = MGvar$proc.moveup, Mvdn = MGvar$proc.movedown, 
        Mvlt = MGvar$proc.moveleft, Mvrt = MGvar$proc.moveright, 
        Zrat = MGvar$proc.zoominrat) {
        if (showtitle == "yes") {
            graphtitle = paste("Procrustes Analysis for", tab1, 
                "on", tab2)
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showpoints1 == "yes") {
            pointtype1 = "p"
        }
        if (showpoints1 == "no") {
            pointtype1 = "n"
        }
        if (showpoints2 == "yes") {
            pointtype2 = "p"
        }
        if (showpoints2 == "no") {
            pointtype2 = "n"
        }
        par(mar = c(3, 3, 3, 3))
        par(cex.axis = 0.8)
        params <- par(bg = bgcol)
        if (nrow(conf1) == 1 && ncol(conf1) == 1) {
            plot(conf1, type = "n", xaxt = "n", yaxt = "n", ylab = "n", 
                xlab = "n", xaxt = "n", yaxt = "n")
            par(params)
        }
        else {
            proc = Procrust(conf1, conf2, labels = dimnames(MGvar$activedata)[[1]])
            Xminproc <- min(proc[, 1])
            Yminproc <- min(proc[, 2])
            Xminconf2 <- min(conf2[, 1])
            Yminconf2 <- min(conf2[, 2])
            Xmin <- min(Xminproc, Xminconf2)
            Ymin <- min(Yminproc, Yminconf2)
            Xmaxproc <- max(proc[, 1])
            Ymaxproc <- max(proc[, 2])
            Xmaxconf2 <- max(conf2[, 1])
            Ymaxconf2 <- max(conf2[, 2])
            Xmax <- max(Xmaxproc, Xmaxconf2)
            Ymax <- max(Ymaxproc, Ymaxconf2)
            Xrange <- Xmax - Xmin
            Yrange <- Ymax - Ymin
            Xmin <- Xmin + (1 - Zrat) * Xrange/2 - (Mvlt - 1) * 
                Xrange + (Mvrt - 1) * Xrange
            Xmax <- Xmax - (1 - Zrat) * Xrange/2 - (Mvlt - 1) * 
                Xrange + (Mvrt - 1) * Xrange
            Ymin <- Ymin + (1 - Zrat) * Yrange/2 - (Mvdn - 1) * 
                Yrange + (Mvup - 1) * Yrange
            Ymax <- Ymax - (1 - Zrat) * Yrange/2 - (Mvdn - 1) * 
                Yrange + (Mvup - 1) * Yrange
            if (showleg == "yes") {
                Ymax = Ymax + 0.1 * Yrange
            }
            if (showleg == "yes") {
                Ymax = Ymax + 0.1 * Yrange
                Xmax = Xmax + 0.05 * Xrange
                Xmin = Xmin - 0.05 * Xrange
            }
            plot(proc, type = "n", ylim = c(Ymin, Ymax), xlim = c(Xmin, 
                Xmax), main = graphtitle, ylab = "", xlab = "", 
                xaxt = xmeas, yaxt = ymeas, fg = axcol)
            if (showleg == "yes") {
                legend(Xmin, 1.09 * Ymax, tab1, pch = pointshape1, 
                  bty = "n", cex = 0.8, pt.bg = "white", col = pointcol1)
                legend(Xmin, 1.03 * Ymax, tab2, pch = pointshape2, 
                  bty = "n", cex = 0.8, pt.bg = "white", col = pointcol2)
                if (showreg1 == "yes") {
                  leglab1 = paste("Regression Line: ", tab1)
                  legend(Xmin + 0.25 * Xrange, 1.09 * Ymax, leglab1, 
                    bty = "n", cex = 0.8, pt.bg = "white", lty = 1, 
                    col = regcol1)
                }
                if (showreg2 == "yes") {
                  leglab2 = paste("Regression Line: ", tab2)
                  legend(Xmin + 0.25 * Xrange, 1.03 * Ymax, leglab2, 
                    bty = "n", cex = 0.8, pt.bg = "white", lty = 1, 
                    col = regcol2)
                }
            }
            points(proc, type = pointtype1, col = pointcol1, 
                cex = pointcex1, pch = pointshape1)
            points(conf2, type = pointtype2, col = pointcol2, 
                cex = pointcex2, pch = pointshape2)
            if (showlabs1 == "yes") {
                text(proc, rownames(proc), cex = pointcex1, col = pointcol1)
            }
            if (showlabs2 == "yes") {
                text(conf2, rownames(conf2), cex = pointcex2, 
                  col = pointcol2)
            }
            if (showreg1 == "yes") {
                regline1 = lm(proc[, 2] ~ proc[, 1])
                abline(regline1, col = regcol1)
            }
            if (showreg2 == "yes") {
                regline2 = lm(conf2[, 2] ~ conf2[, 1])
                abline(regline2, col = regcol2)
            }
            mtext(text = xlabel, side = 1, line = 1.9, adj = 0.5, 
                cex = 1, col = "black")
            mtext(text = ylabel, side = 2, line = 1.9, adj = 0.5, 
                cex = 1, col = "black")
            params <- par(bg = bgcol)
            MGvar$Proc.parPlotSize <<- par("plt")
            MGvar$Proc.usrCoords <<- par("usr")
            MGvar$Proc.labelsVec <<- rownames(conf1)
            if (length(MGvar$Proc.indexLabeled) > 0) {
                for (i in (1:length(MGvar$Proc.indexLabeled))) {
                  indexClosest <- MGvar$Proc.indexLabeled[i]
                  text(proc[indexClosest, 1], proc[indexClosest, 
                    2], labels = MGvar$Proc.labelsVec[indexClosest], 
                    pos = MGvar$GenSet.ILPos, col = pointcol1)
                  text(conf2[indexClosest, 1], conf2[indexClosest, 
                    2], labels = MGvar$Proc.labelsVec[indexClosest], 
                    pos = MGvar$GenSet.ILPos, col = pointcol2)
                }
            }
        }
    }
    Proc.labelClosestPoint <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance.P <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        indexClosest <- which.min(squared.Distance.P)
        MGvar$Proc.indexLabeled <<- c(MGvar$Proc.indexLabeled, 
            indexClosest)
        tkrreplot(procimg)
        if (MGvar$EnProcPlot.switch == "on") {
            tkrreplot(MGcomp$POprocimg)
        }
    }
    Proc.OnLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", procimg)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", procimg)))
        xMin <- MGvar$Proc.parPlotSize[1] * width
        xMax <- MGvar$Proc.parPlotSize[2] * width
        yMin <- MGvar$Proc.parPlotSize[3] * height
        yMax <- MGvar$Proc.parPlotSize[4] * height
        rangeX <- MGvar$Proc.usrCoords[2] - MGvar$Proc.usrCoords[1]
        rangeY <- MGvar$Proc.usrCoords[4] - MGvar$Proc.usrCoords[3]
        imgXcoords <- (ProcConf2[, 1] - MGvar$Proc.usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (ProcConf2[, 2] - MGvar$Proc.usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$Proc.usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$Proc.usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        Proc.labelClosestPoint(xClick, yClick, imgXcoords, imgYcoords)
    }
    Proc.LabelSpecificPoint <- function() {
        Labtt = tktoplevel()
        tkwm.resizable(Labtt, "0", "0")
        tkwm.deiconify(Labtt)
        tkwm.title(Labtt, "Procrutes Point Label")
        tkwm.geometry(Labtt, "280x100")
        LabCanvas = tkcanvas(Labtt, width = 280, height = 100, 
            bg = col.sec)
        tkplace(LabCanvas, relx = 0, `in` = Labtt)
        frameLab <- tkwidget(Labtt, "TitleFrame", text = "Label a Specific Point", 
            background = "white")
        tkplace(frameLab, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheight = 0.96, `in` = Labtt)
        tkplace(tklabel(frameLab, text = "Enter the Point Name", 
            background = "white"), relx = 0.05, rely = 0.3, `in` = frameLab)
        LabInput.val = tclVar("")
        LabInput = tkentry(Labtt, width = 15, textvariable = LabInput.val)
        tkplace(LabInput, relx = 0.6, rely = 0.3, `in` = frameLab)
        On.Enter <- function() {
            Lab = as.character(tclvalue(LabInput.val))
            datlength = nrow(ProcConf2)
            indexlength = length(MGvar$Proc.indexLabeled)
            pointpos = 0
            for (i in 1:datlength) {
                if (rownames(ProcConf2)[i] == Lab) {
                  pointpos = i
                }
            }
            if (pointpos == 0) {
                tkmessageBox(title = "Point Label Error", message = "The name does not match any of the points. \nPlease check spelling", 
                  icon = "error", type = "ok")
            }
            if (pointpos > 0) {
                MGvar$Proc.indexLabeled[indexlength + 1] <<- pointpos
                names(MGvar$Proc.indexLabeled)[indexlength + 
                  1] <<- rownames(ProcConf2)[pointpos]
            }
            tkrreplot(procimg)
            tkdestroy(Labtt)
        }
        tkplace(tkbutton(Labtt, text = "Enter", width = 15, command = function() On.Enter()), 
            relx = 0.3, rely = 0.7, `in` = frameLab)
        tkfocus(Labtt)
        tkbind(Labtt, "<Return>", On.Enter)
    }
    EnlargedProcPlot <- function() {
        ReusePOW()
        MGvar$EnProcPlot.switch <<- "on"
        MGcomp$EProc <<- tktoplevel()
        tkwm.title(MGcomp$EProc, "Enlarged Procrustes Plot")
        tkwm.geometry(MGcomp$EProc, "650x650")
        EProcCanvas <- tkcanvas(MGcomp$EProc, width = 650, height = 650, 
            bg = col.sec)
        tkplace(EProcCanvas, `in` = MGcomp$EProc)
        POhscale <- 1.5
        POvscale <- 1.4
        MGcomp$POprocimg <<- tkrplot(MGcomp$EProc, function() ProcrustesMDSGUIplot(conf1 = ProcConf1, 
            conf2 = ProcConf2), hscale = POhscale, vscale = POvscale)
        tkplace(MGcomp$POprocimg, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            `in` = MGcomp$EProc)
        tkbind(MGcomp$POprocimg, "<Button-1>", OnPOProcLeftClick)
        tkconfigure(MGcomp$POprocimg, cursor = "hand2")
        CopyAPToClip <- function() {
            tkrreplot(MGcomp$POprocimg)
        }
        tkplace(tkbutton(MGcomp$EProc, text = "Copy to Clipboard", 
            width = 20, command = function() CopyAPToClip()), 
            relx = 0.19, rely = 0.93, `in` = MGcomp$EProc)
        tkplace(tkbutton(MGcomp$EProc, text = "Plot Options", 
            width = 20, command = function() ProcPlotOps()), 
            relx = 0.6, rely = 0.93, `in` = MGcomp$EProc)
        tkbind(MGcomp$POprocimg, "<Destroy>", function() {
            MGvar$EnActivePlot.switch.T1 <<- "off"
        })
        ProcPlotPOPMenu <- tkmenu(MGcomp$POprocimg, tearoff = FALSE)
        tkadd(ProcPlotPOPMenu, "command", label = "Clear Added Point Labels", 
            command = ClearProcPoints)
        RightClickProcPOP <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", MGcomp$POprocimg))
            rooty <- as.integer(tkwinfo("rooty", MGcomp$POprocimg))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", ProcPlotPOPMenu, xTxt, yTxt)
        }
        tkbind(MGcomp$POprocimg, "<Button-3>", RightClickProcPOP)
        tkbind(MGcomp$EProc, "<Configure>", resize.EProc)
    }
    resize.EProc <- function() {
        EP.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EProc)))
        EP.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EProc)))
        WidthScale = 650/1.5
        EShscale <- EP.width/WidthScale
        dimrat = 1.5/1.4
        ESvscale <- EShscale/dimrat
        tkrreplot(MGcomp$POprocimg, hscale = EShscale, vscale = ESvscale)
    }
    OnPOProcLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", MGcomp$POprocimg)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", MGcomp$POprocimg)))
        xMin <- MGvar$Proc.parPlotSize[1] * width
        xMax <- MGvar$Proc.parPlotSize[2] * width
        yMin <- MGvar$Proc.parPlotSize[3] * height
        yMax <- MGvar$Proc.parPlotSize[4] * height
        rangeX <- MGvar$Proc.usrCoords[2] - MGvar$Proc.usrCoords[1]
        rangeY <- MGvar$Proc.usrCoords[4] - MGvar$Proc.usrCoords[3]
        imgXcoords <- (ProcConf2[, 1] - MGvar$Proc.usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (ProcConf2[, 2] - MGvar$Proc.usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$Proc.usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$Proc.usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        msg <- paste("Label the point closest to these approximate plot coordinates: \n", 
            "x =", format(xPlotCoord, digits = 2), ",y =", format(yPlotCoord, 
                digits = 2), "?")
        mbval <- tkmessageBox(title = "Label Point Closest to These Approximate Plot Coordinates", 
            message = msg, type = "yesno", icon = "question")
        if (tclvalue(mbval) == "yes") 
            Proc.labelClosestPoint(xClick, yClick, imgXcoords, 
                imgYcoords)
    }
    savetext <- function() {
        file <- tclvalue(tkgetSaveFile(initialfile = tclvalue(tclfile.tail(wfile)), 
            initialdir = tclvalue(tclfile.dir(wfile))))
        if (!length(file)) 
            return()
        chn <- tclopen(file, "w")
        tclputs(chn, tclvalue(tkget(Notetxt, "0.0", "end")))
        tclclose(chn)
        wfile <<- file
    }
    loadtext <- function() {
        file <- tclvalue(tkgetOpenFile())
        if (!length(file)) 
            return()
        chn <- tclopen(file, "r")
        tkinsert(Notetxt, "0.0", tclvalue(tclread(chn)))
        tclclose(chn)
        wfile <<- file
    }
    runtextascode <- function() {
        code <- tclvalue(tkget(Notetxt, "0.0", "end"))
        e <- try(parse(text = code))
        if (inherits(e, "try-error")) {
            tkmessageBox(message = "Syntax error", icon = "error")
            return()
        }
        cat("Executing from MDSGUI note tab window:", "-----", 
            code, "result:", sep = "\n")
        print(eval(e))
    }
    tableplace <- function() {
        myConfarray <- c("Tab", "Plot1", "Plot2", "Plot3", "Plot4", 
            "Plot5", "Stat 3D", "RGL 3D", "Measure", MGvar$dMeas.T1, 
            MGvar$dMeas.T2, MGvar$dMeas.T3, MGvar$dMeas.T4, MGvar$dMeas.T5, 
            MGvar$dMeas.3S, MGvar$dMeas.3R, "MDS", MGvar$MDStype.T1, 
            MGvar$MDStype.T2, MGvar$MDStype.T3, MGvar$MDStype.T4, 
            MGvar$MDStype.T5, MGvar$MDStype.3S, MGvar$MDStype.3R, 
            "Dims", MGvar$MDS.dimensions.T1, MGvar$MDS.dimensions.T2, 
            MGvar$MDS.dimensions.T3, MGvar$MDS.dimensions.T4, 
            MGvar$MDS.dimensions.T5, MGvar$MDS.dimensions.3S, 
            MGvar$MDS.dimensions.3R, MGvar$StCalc, MGvar$MDSStress.T1, 
            MGvar$MDSStress.T2, MGvar$MDSStress.T3, MGvar$MDSStress.T4, 
            MGvar$MDSStress.T5, MGvar$MDSStress.3S, MGvar$MDSStress.3R, 
            "Plot.Dims", MGvar$TabDims.T1, MGvar$TabDims.T2, 
            MGvar$TabDims.T3, MGvar$TabDims.T4, MGvar$TabDims.T5, 
            MGvar$TabDims.3S, MGvar$TabDims.3R, "Tolerence", 
            MGvar$MDS.tol.T1, MGvar$MDS.tol.T2, MGvar$MDS.tol.T3, 
            MGvar$MDS.tol.T4, MGvar$MDS.tol.T5, MGvar$MDS.tol.3S, 
            MGvar$MDS.tol.3R, "Iterations", MGvar$MDS.iter.T1, 
            MGvar$MDS.iter.T2, MGvar$MDS.iter.T3, MGvar$MDS.iter.T4, 
            MGvar$MDS.iter.T5, MGvar$MDS.iter.3S, MGvar$MDS.iter.3R)
        MGcomp$Confarray <<- tclArray()
        dim(myConfarray) <- c(8, 8)
        for (i in 0:7) for (j in 0:7) {
            MGcomp$Confarray[[i, j]] <- myConfarray[i + 1, j + 
                1]
        }
        table.Conf <- tkwidget(frameConfTab, "table", rows = 8, 
            cols = 8, variable = MGcomp$Confarray, xscrollcommand = function(...) tkset(tablescr, 
                ...), titlerows = "1", titlecols = "1", selectmode = "extended", 
            background = "white", width = 6, height = 8)
        tkgrid(table.Conf)
        tcl(table.Conf, "width", 0, 6)
        tcl(table.Conf, "width", 1, 19)
        tcl(table.Conf, "width", 2, 15)
        tcl(table.Conf, "width", 3, 7)
        tcl(table.Conf, "width", 4, 11)
        tcl(table.Conf, "width", 5, 11)
        tablescr <- tkscrollbar(frameConfTab, orient = "horizontal", 
            command = function(...) tkxview(table.Conf, ...))
        tkgrid(tablescr, sticky = "new", columnspan = 6)
        tkgrid.configure(tablescr, sticky = "new", columnspan = 5)
        tkconfigure(table.Conf, variable = MGcomp$Confarray, 
            background = "white", selectmode = "extended")
        MGcomp$table.Conf <<- table.Conf
    }
    tableupdate <- function() {
        myConfarray <- c("Tab", "Plot1", "Plot2", "Plot3", "Plot4", 
            "Plot5", "Stat 3D", "RGL 3D", "Measure", MGvar$dMeas.T1, 
            MGvar$dMeas.T2, MGvar$dMeas.T3, MGvar$dMeas.T4, MGvar$dMeas.T5, 
            MGvar$dMeas.3S, MGvar$dMeas.3R, "MDS", MGvar$MDStype.T1, 
            MGvar$MDStype.T2, MGvar$MDStype.T3, MGvar$MDStype.T4, 
            MGvar$MDStype.T5, MGvar$MDStype.3S, MGvar$MDStype.3R, 
            "Dims", MGvar$MDS.dimensions.T1, MGvar$MDS.dimensions.T2, 
            MGvar$MDS.dimensions.T3, MGvar$MDS.dimensions.T4, 
            MGvar$MDS.dimensions.T5, MGvar$MDS.dimensions.3S, 
            MGvar$MDS.dimensions.3R, MGvar$StCalc, MGvar$MDSStress.T1, 
            MGvar$MDSStress.T2, MGvar$MDSStress.T3, MGvar$MDSStress.T4, 
            MGvar$MDSStress.T5, MGvar$MDSStress.3S, MGvar$MDSStress.3R, 
            "Plot.Dims", MGvar$TabDims.T1, MGvar$TabDims.T2, 
            MGvar$TabDims.T3, MGvar$TabDims.T4, MGvar$TabDims.T5, 
            MGvar$TabDims.3S, MGvar$TabDims.3R, "Tolerence", 
            MGvar$MDS.tol.T1, MGvar$MDS.tol.T2, MGvar$MDS.tol.T3, 
            MGvar$MDS.tol.T4, MGvar$MDS.tol.T5, MGvar$MDS.tol.3S, 
            MGvar$MDS.tol.3R, "Iterations", MGvar$MDS.iter.T1, 
            MGvar$MDS.iter.T2, MGvar$MDS.iter.T3, MGvar$MDS.iter.T4, 
            MGvar$MDS.iter.T5, MGvar$MDS.iter.3S, MGvar$MDS.iter.3R)
        MGcomp$Confarray <<- tclArray()
        dim(myConfarray) <- c(8, 8)
        for (i in 0:7) for (j in 0:7) {
            MGcomp$Confarray[[i, j]] <- myConfarray[i + 1, j + 
                1]
        }
        tkconfigure(MGcomp$table.Conf, variable = MGcomp$Confarray)
    }
    GOtoActiveTab <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tk2notetab.select(myPlottingNB, "Plot1")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tk2notetab.select(myPlottingNB, "Plot2")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tk2notetab.select(myPlottingNB, "Plot3")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tk2notetab.select(myPlottingNB, "Plot4")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tk2notetab.select(myPlottingNB, "Plot5")
        }
    }
    Incomplete = function() {
        tkmessageBox(message = "this code is incomplete!", icon = "error", 
            type = "ok")
    }
    LargeDataOps <- function() {
        LDO = tktoplevel()
        tkwm.resizable(LDO, "0", "0")
        tkwm.deiconify(LDO)
        tkwm.title(LDO, "Large Data")
        tkwm.geometry(LDO, "350x250")
        LDOcanvas = tkcanvas(LDO, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(LDOcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = LDO)
        frameLDO <- tkwidget(LDO, "TitleFrame", text = "Large Data Options", 
            background = "white")
        tkplace(frameLDO, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheight = 0.96, `in` = LDO)
        tkplace(tklabel(frameLDO, text = "The Data that you have uploaded has a number of objects\n greater than 50. Some processes of the MDS-GUI may\n run very slowly with data of this size. You are advised\n to deactivate the following.", 
            background = "white"), relx = 0.05, rely = 0.15, 
            `in` = frameLDO)
        tkplace(tklabel(frameLDO, text = "Scree Plot", background = "white"), 
            relx = 0.15, rely = 0.5, `in` = frameLDO)
        DScree = tclVar(1)
        DScree.CB <- tk2checkbutton(LDO)
        tkconfigure(DScree.CB, variable = DScree)
        tkplace(DScree.CB, relx = 0.75, rely = 0.5, `in` = frameLDO)
        tkplace(tklabel(frameLDO, text = "Shepard Plot", background = "white"), 
            relx = 0.15, rely = 0.6, `in` = frameLDO)
        DShep = tclVar(1)
        DShep.CB <- tk2checkbutton(LDO)
        tkconfigure(DShep.CB, variable = DShep)
        tkplace(DShep.CB, relx = 0.75, rely = 0.6, `in` = frameLDO)
        On.Deactivate <- function() {
            if (as.numeric(tclvalue(DScree)) == 1) {
                MGvar$GenSet.CalcScree <<- "no"
                tclvalue(MGvar$GS.Scree.val) <<- 0
            }
            if (as.numeric(tclvalue(DShep)) == 1) {
                MGvar$GenSet.CalcShep <<- "no"
                tclvalue(MGvar$GS.Shep.val) <<- 0
            }
            tkdestroy(LDO)
        }
        tkplace(tkbutton(LDO, text = "Deactivate", width = 15, 
            command = function() On.Deactivate()), relx = 0.35, 
            rely = 0.8, `in` = LDO)
        tkfocus(LDO)
        tkbind(LDO, "<Return>", On.Deactivate)
        tkwait.window(LDO)
    }
    LoadDataSettxt = function() {
        fileName <- tclvalue(tkgetOpenFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            loadeddata = read.table(fileName)
            MGvar$activedata <<- loadeddata
            namingtt <- tktoplevel()
            tkwm.resizable(namingtt, "0", "0")
            tkwm.deiconify(namingtt)
            tkwm.title(namingtt, "New Active Dataset Options")
            tkwm.geometry(namingtt, "350x400")
            Loadcanvas = tkcanvas(namingtt, width = "1128", height = "756", 
                bg = col.sec)
            tkplace(Loadcanvas, relx = 0, rely = 0, relwidth = 1, 
                relheight = 1, `in` = namingtt)
            frameNaming <- tkwidget(namingtt, "TitleFrame", text = "Dataset Name", 
                background = "white")
            tkplace(frameNaming, relx = 0.02, rely = 0.02, relwidth = 0.96, 
                relheight = 0.16, `in` = namingtt)
            tkplace(tklabel(frameNaming, text = "Enter the name of your DataSet", 
                background = "white"), relx = 0.08, rely = 0.4, 
                `in` = frameNaming)
            ChangingName = tclVar("")
            NewNameBox = tkentry(namingtt, width = 15, textvariable = ChangingName)
            tkplace(NewNameBox, relx = 0.65, rely = 0.4, `in` = frameNaming)
            frameTrans <- tkwidget(namingtt, "TitleFrame", text = "Dataset Transpose", 
                background = "white")
            tkplace(frameTrans, relx = 0.02, rely = 0.19, relwidth = 0.96, 
                relheight = 0.3, `in` = namingtt)
            fontsmall <- tkfont.create(family = "times", size = 9)
            tkplace(tklabel(frameTrans, text = "All procedures in this package require that the active data have\n objects as rows and variables as columns. If your data is not\n in this format then please transpose.", 
                font = fontsmall), relx = 0.05, rely = 0.15, 
                `in` = frameTrans)
            tkplace(tklabel(frameTrans, text = "Transpose Active Data?", 
                background = "white"), relx = 0.08, rely = 0.65, 
                `in` = frameTrans)
            cbtrans <- tk2checkbutton(namingtt)
            cbTValue <- tclVar("0")
            tkconfigure(cbtrans, variable = cbTValue)
            tkplace(cbtrans, relx = 0.75, rely = 0.65, `in` = frameTrans)
            frameScale <- tkwidget(namingtt, "TitleFrame", text = "Dataset Scale", 
                background = "white")
            tkplace(frameScale, relx = 0.02, rely = 0.5, relwidth = 0.96, 
                relheight = 0.2, `in` = namingtt)
            tkplace(tklabel(frameScale, text = "Scaling Data will scale the columns of your data between 0 and 1.", 
                font = fontsmall), relx = 0.04, rely = 0.25, 
                `in` = frameScale)
            tkplace(tklabel(frameScale, text = "Scale your active data?", 
                background = "white"), relx = 0.08, rely = 0.6, 
                `in` = frameScale)
            ScDat.val <- tclVar(0)
            ScDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ScDat.CB, variable = ScDat.val)
            tkplace(ScDat.CB, relx = 0.75, rely = 0.6, `in` = frameScale)
            frameCol <- tkwidget(namingtt, "TitleFrame", text = "Dataset Colours", 
                background = "white")
            tkplace(frameCol, relx = 0.02, rely = 0.71, relwidth = 0.96, 
                relheight = 0.2, `in` = namingtt)
            tkplace(tklabel(frameCol, text = "Does the data contain a column of object category information?", 
                font = fontsmall), relx = 0.02, rely = 0.25, 
                `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Yes", background = "white"), 
                relx = 0.1, rely = 0.6, `in` = frameCol)
            ColDat <- tclVar(0)
            ColDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ColDat.CB, variable = ColDat)
            tkplace(ColDat.CB, relx = 0.25, rely = 0.6, `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Which Column?", 
                background = "white"), relx = 0.4, rely = 0.6, 
                `in` = frameCol)
            ColColumn <- tclVar("First")
            ColColumn.ComboBox <- tkwidget(namingtt, "ComboBox", 
                editable = FALSE, values = c("First", "Last"), 
                width = 6, textvariable = ColColumn)
            tkplace(ColColumn.ComboBox, relx = 0.7, rely = 0.6, 
                `in` = frameCol)
            OnOk.Name <- function() {
                MGvar$ClasTabCols <<- c()
                NName = as.character(tclvalue(ChangingName))
                MGvar$datatitle <<- NName
                tclvalue(labelText) <- paste("Active Dataset is", 
                  MGvar$datatitle)
                tkentryconfigure(dataMenu, 7, state = "disabled")
                YesCol <- as.character(tclvalue(ColDat))
                if (YesCol == "1") {
                  if (as.character(tclvalue(ColColumn)) == "First") {
                    MGvar$ClasVec <<- MGvar$activedata[, 1]
                    MGvar$activedata <<- MGvar$activedata[, -1]
                  }
                  if (as.character(tclvalue(ColColumn)) == "Last") {
                    MGvar$ClasVec <<- MGvar$activedata[, ncol(MGvar$activedata)]
                    MGvar$activedata <<- MGvar$activedata[, -ncol(MGvar$activedata)]
                  }
                  MGvar$ClasTab <<- table(MGvar$ClasVec)
                  potcols <- c(brewer.pal(9, "Set1"), brewer.pal(12, 
                    "Paired"))
                  for (i in 1:nrow(MGvar$activedata)) {
                    for (j in 1:length(MGvar$ClasTab)) {
                      if (MGvar$ClasVec[i] == names(MGvar$ClasTab)[j]) {
                        if (j == length(potcols)) {
                          num = length(potcols)
                        }
                        if (j != length(potcols)) {
                          num = j%%length(potcols)
                        }
                        MGvar$MDSmat.Cols[i] <<- potcols[num]
                      }
                    }
                  }
                  for (i in 1:length(MGvar$ClasTab)) {
                    MGvar$ClasTabCols <<- c(MGvar$ClasTabCols, 
                      potcols[i])
                    names(MGvar$ClasTabCols)[i] <<- names(MGvar$ClasTab)[i]
                  }
                  MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
                  tkentryconfigure(dataMenu, 7, state = "active")
                }
                else {
                  FirstPointColInitialise()
                  PointColInitialise()
                }
                TOp = as.character(tclvalue(cbTValue))
                if (TOp == "1") {
                  MGvar$activedata <<- t(MGvar$activedata)
                  tkmessageBox(message = paste("Your data has been transposed"))
                }
                Scl <- as.character(tclvalue(ScDat.val))
                if (Scl == "1") {
                  MGvar$activedata <<- stand.1.range(MGvar$activedata)
                }
                if (is.numeric(as.matrix(MGvar$activedata))) {
                  if (0 %in% as.matrix(MGvar$activedata)) {
                    tkentryconfigure(MDSMenuDistCal, 10, state = "disabled")
                    if (MGvar$dMeas == "Wave-Hedges") {
                      MGvar$dMeas <<- "Euclidean.Distance"
                    }
                    tkentryconfigure(MDSMenuDistCal, 6, state = "disabled")
                    if (MGvar$dMeas == "Divergence") {
                      MGvar$dMeas <<- "Euclidean.Distance"
                    }
                  }
                  else {
                    tkentryconfigure(MDSMenuDistCal, 10, state = "active")
                    tkentryconfigure(MDSMenuDistCal, 6, state = "active")
                  }
                  new.dissim.meas(MGvar$activedata)
                  MGvar$originaldistmat <<- MGvar$distmat
                  MGvar$distmat.T1 <<- MGvar$distmat
                  MGvar$distmat.T2 <<- MGvar$distmat
                  MGvar$distmat.T3 <<- MGvar$distmat
                  MGvar$distmat.T4 <<- MGvar$distmat
                  MGvar$distmat.T5 <<- MGvar$distmat
                  MGvar$maxdims <<- nrow(MGvar$distmat) - 1
                  MGvar$removedpointsactivedata <<- MGvar$activedata
                  tkentryconfigure(functionMenu, 5, state = "active")
                  tkentryconfigure(MainPlotMenu1, 15, state = "active")
                  tkentryconfigure(MainPlotMenu2, 15, state = "active")
                  tkentryconfigure(MainPlotMenu3, 15, state = "active")
                  tkentryconfigure(MainPlotMenu4, 15, state = "active")
                  tkentryconfigure(MainPlotMenu5, 15, state = "active")
                }
                else {
                  tkmessageBox(message = "Non-Numeric values in the data. Please retry data upload.", 
                    type = "ok")
                  MGvar$activedata <<- as.matrix(0)
                  tclvalue(labelText) <- paste("No Active Dataset")
                  MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                  tkentryconfigure(dataMenu, 7, state = "disabled")
                }
                tkdestroy(namingtt)
                if (nrow(MGvar$activedata) > 50) {
                  LargeDataOps()
                }
            }
            tkplace(tkbutton(namingtt, text = "OK", width = 10, 
                command = function() OnOk.Name()), relx = 0.38, 
                rely = 0.925, `in` = namingtt)
            tkbind(namingtt, "<Return>", OnOk.Name)
            tkentryconfigure(EditDataMenu, 0, state = "active")
            MGvar$tShepx <<- as.vector(0)
            ClearRemindex()
        }
    }
	LoadDataSetcsv = function() {
        fileName <- tclvalue(tkgetOpenFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            loadeddata = read.csv(fileName)
            MGvar$activedata <<- loadeddata
		rownames(MGvar$activedata) <<- MGvar$activedata[,1]
            MGvar$activedata <<- MGvar$activedata[, -1]
            namingtt <- tktoplevel()
            tkwm.resizable(namingtt, "0", "0")
            tkwm.deiconify(namingtt)
            tkwm.title(namingtt, "New Active Dataset Options")
            tkwm.geometry(namingtt, "350x400")
            Loadcanvas = tkcanvas(namingtt, width = "1128", height = "756", 
                bg = col.sec)
            tkplace(Loadcanvas, relx = 0, rely = 0, relwidth = 1, 
                relheight = 1, `in` = namingtt)
            frameNaming <- tkwidget(namingtt, "TitleFrame", text = "Dataset Name", 
                background = "white")
            tkplace(frameNaming, relx = 0.02, rely = 0.02, relwidth = 0.96, 
                relheight = 0.16, `in` = namingtt)
            tkplace(tklabel(frameNaming, text = "Enter the name of your DataSet", 
                background = "white"), relx = 0.08, rely = 0.4, 
                `in` = frameNaming)
            ChangingName = tclVar("")
            NewNameBox = tkentry(namingtt, width = 15, textvariable = ChangingName)
            tkplace(NewNameBox, relx = 0.65, rely = 0.4, `in` = frameNaming)
            frameTrans <- tkwidget(namingtt, "TitleFrame", text = "Dataset Transpose", 
                background = "white")
            tkplace(frameTrans, relx = 0.02, rely = 0.19, relwidth = 0.96, 
                relheight = 0.3, `in` = namingtt)
            fontsmall <- tkfont.create(family = "times", size = 9)
            tkplace(tklabel(frameTrans, text = "All procedures in this package require that the active data have\n objects as rows and variables as columns. If your data is not\n in this format then please transpose.", 
                font = fontsmall), relx = 0.05, rely = 0.15, 
                `in` = frameTrans)
            tkplace(tklabel(frameTrans, text = "Transpose Active Data?", 
                background = "white"), relx = 0.08, rely = 0.65, 
                `in` = frameTrans)
            cbtrans <- tk2checkbutton(namingtt)
            cbTValue <- tclVar("0")
            tkconfigure(cbtrans, variable = cbTValue)
            tkplace(cbtrans, relx = 0.75, rely = 0.65, `in` = frameTrans)
            frameScale <- tkwidget(namingtt, "TitleFrame", text = "Dataset Scale", 
                background = "white")
            tkplace(frameScale, relx = 0.02, rely = 0.5, relwidth = 0.96, 
                relheight = 0.2, `in` = namingtt)
            tkplace(tklabel(frameScale, text = "Scaling Data will scale the columns of your data between 0 and 1.", 
                font = fontsmall), relx = 0.04, rely = 0.25, 
                `in` = frameScale)
            tkplace(tklabel(frameScale, text = "Scale your active data?", 
                background = "white"), relx = 0.08, rely = 0.6, 
                `in` = frameScale)
            ScDat.val <- tclVar(0)
            ScDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ScDat.CB, variable = ScDat.val)
            tkplace(ScDat.CB, relx = 0.75, rely = 0.6, `in` = frameScale)
            frameCol <- tkwidget(namingtt, "TitleFrame", text = "Dataset Colours", 
                background = "white")
            tkplace(frameCol, relx = 0.02, rely = 0.71, relwidth = 0.96, 
                relheight = 0.2, `in` = namingtt)
            tkplace(tklabel(frameCol, text = "Does the data contain a column of object category information?", 
                font = fontsmall), relx = 0.02, rely = 0.25, 
                `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Yes", background = "white"), 
                relx = 0.1, rely = 0.6, `in` = frameCol)
            ColDat <- tclVar(0)
            ColDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ColDat.CB, variable = ColDat)
            tkplace(ColDat.CB, relx = 0.25, rely = 0.6, `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Which Column?", 
                background = "white"), relx = 0.4, rely = 0.6, 
                `in` = frameCol)
            ColColumn <- tclVar("First")
            ColColumn.ComboBox <- tkwidget(namingtt, "ComboBox", 
                editable = FALSE, values = c("First", "Last"), 
                width = 6, textvariable = ColColumn)
            tkplace(ColColumn.ComboBox, relx = 0.7, rely = 0.6, 
                `in` = frameCol)
            OnOk.Name <- function() {
                MGvar$ClasTabCols <<- c()
                NName = as.character(tclvalue(ChangingName))
                MGvar$datatitle <<- NName
                tclvalue(labelText) <- paste("Active Dataset is", 
                  MGvar$datatitle)
                tkentryconfigure(dataMenu, 7, state = "disabled")
                YesCol <- as.character(tclvalue(ColDat))
                if (YesCol == "1") {
                  if (as.character(tclvalue(ColColumn)) == "First") {
                    MGvar$ClasVec <<- MGvar$activedata[, 1]
                    MGvar$activedata <<- MGvar$activedata[, -1]
                  }
                  if (as.character(tclvalue(ColColumn)) == "Last") {
                    MGvar$ClasVec <<- MGvar$activedata[, ncol(MGvar$activedata)]
                    MGvar$activedata <<- MGvar$activedata[, -ncol(MGvar$activedata)]
                  }
                  MGvar$ClasTab <<- table(MGvar$ClasVec)
                  potcols <- c(brewer.pal(9, "Set1"), brewer.pal(12, 
                    "Paired"))
                  for (i in 1:nrow(MGvar$activedata)) {
                    for (j in 1:length(MGvar$ClasTab)) {
                      if (MGvar$ClasVec[i] == names(MGvar$ClasTab)[j]) {
                        if (j == length(potcols)) {
                          num = length(potcols)
                        }
                        if (j != length(potcols)) {
                          num = j%%length(potcols)
                        }
                        MGvar$MDSmat.Cols[i] <<- potcols[num]
                      }
                    }
                  }
                  for (i in 1:length(MGvar$ClasTab)) {
                    MGvar$ClasTabCols <<- c(MGvar$ClasTabCols, 
                      potcols[i])
                    names(MGvar$ClasTabCols)[i] <<- names(MGvar$ClasTab)[i]
                  }
                  MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
                  tkentryconfigure(dataMenu, 7, state = "active")
                }
                else {
                  FirstPointColInitialise()
                  PointColInitialise()
                }
                TOp = as.character(tclvalue(cbTValue))
                if (TOp == "1") {
                  MGvar$activedata <<- t(MGvar$activedata)
                  tkmessageBox(message = paste("Your data has been transposed"))
                }
                Scl <- as.character(tclvalue(ScDat.val))
                if (Scl == "1") {
                  MGvar$activedata <<- stand.1.range(MGvar$activedata)
                }
                if (is.numeric(as.matrix(MGvar$activedata))) {
                  if (0 %in% as.matrix(MGvar$activedata)) {
                    tkentryconfigure(MDSMenuDistCal, 10, state = "disabled")
                    if (MGvar$dMeas == "Wave-Hedges") {
                      MGvar$dMeas <<- "Euclidean.Distance"
                    }
                    tkentryconfigure(MDSMenuDistCal, 6, state = "disabled")
                    if (MGvar$dMeas == "Divergence") {
                      MGvar$dMeas <<- "Euclidean.Distance"
                    }
                  }
                  else {
                    tkentryconfigure(MDSMenuDistCal, 10, state = "active")
                    tkentryconfigure(MDSMenuDistCal, 6, state = "active")
                  }
                  new.dissim.meas(MGvar$activedata)
                  MGvar$originaldistmat <<- MGvar$distmat
                  MGvar$distmat.T1 <<- MGvar$distmat
                  MGvar$distmat.T2 <<- MGvar$distmat
                  MGvar$distmat.T3 <<- MGvar$distmat
                  MGvar$distmat.T4 <<- MGvar$distmat
                  MGvar$distmat.T5 <<- MGvar$distmat
                  MGvar$maxdims <<- nrow(MGvar$distmat) - 1
                  MGvar$removedpointsactivedata <<- MGvar$activedata
                  tkentryconfigure(functionMenu, 5, state = "active")
                  tkentryconfigure(MainPlotMenu1, 15, state = "active")
                  tkentryconfigure(MainPlotMenu2, 15, state = "active")
                  tkentryconfigure(MainPlotMenu3, 15, state = "active")
                  tkentryconfigure(MainPlotMenu4, 15, state = "active")
                  tkentryconfigure(MainPlotMenu5, 15, state = "active")
                }
                else {
                  tkmessageBox(message = "Non-Numeric values in the data. Please retry data upload.", 
                    type = "ok")
                  MGvar$activedata <<- as.matrix(0)
                  tclvalue(labelText) <- paste("No Active Dataset")
                  MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                  tkentryconfigure(dataMenu, 7, state = "disabled")
                }
                tkdestroy(namingtt)
                if (nrow(MGvar$activedata) > 50) {
                  LargeDataOps()
                }
            }
            tkplace(tkbutton(namingtt, text = "OK", width = 10, 
                command = function() OnOk.Name()), relx = 0.38, 
                rely = 0.925, `in` = namingtt)
            tkbind(namingtt, "<Return>", OnOk.Name)
            tkentryconfigure(EditDataMenu, 0, state = "active")
            MGvar$tShepx <<- as.vector(0)
            ClearRemindex()
        }
    }
    LoadDistSet <- function() {
        fileName <- tclvalue(tkgetOpenFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            loadeddata = read.table(fileName)
            MGvar$activedata <<- loadeddata
            MGvar$distmat <<- loadeddata
            namingtt <- tktoplevel()
            tkwm.resizable(namingtt, "0", "0")
            tkwm.deiconify(namingtt)
            tkwm.title(namingtt, "New Active Dissimilarity Matrix Options")
            tkwm.geometry(namingtt, "350x300")
            Loadcanvas = tkcanvas(namingtt, width = "1128", height = "756", 
                bg = col.sec)
            tkplace(Loadcanvas, relx = 0, rely = 0, relwidth = 1, 
                relheight = 1, `in` = namingtt)
            frameNaming <- tkwidget(namingtt, "TitleFrame", text = "Dataset Name", 
                background = "white")
            tkplace(frameNaming, relx = 0.02, rely = 0.02, relwidth = 0.96, 
                relheight = 0.25, `in` = namingtt)
            tkplace(tklabel(frameNaming, text = "Enter the name of your DataSet", 
                background = "white"), relx = 0.08, rely = 0.4, 
                `in` = frameNaming)
            ChangingName = tclVar("")
            NewNameBox = tkentry(namingtt, width = 15, textvariable = ChangingName)
            tkplace(NewNameBox, relx = 0.65, rely = 0.4, `in` = frameNaming)
            fontsmall <- tkfont.create(family = "times", size = 9)
            frameScale <- tkwidget(namingtt, "TitleFrame", text = "Dataset Scale", 
                background = "white")
            tkplace(frameScale, relx = 0.02, rely = 0.29, relwidth = 0.96, 
                relheight = 0.25, `in` = namingtt)
            tkplace(tklabel(frameScale, text = "Scaling Data will scale the columns of your data between 0 and 1.", 
                font = fontsmall), relx = 0.04, rely = 0.25, 
                `in` = frameScale)
            tkplace(tklabel(frameScale, text = "Scale your active data?", 
                background = "white"), relx = 0.08, rely = 0.6, 
                `in` = frameScale)
            ScDat.val <- tclVar(0)
            ScDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ScDat.CB, variable = ScDat.val)
            tkplace(ScDat.CB, relx = 0.75, rely = 0.6, `in` = frameScale)
            frameCol <- tkwidget(namingtt, "TitleFrame", text = "Dataset Colours", 
                background = "white")
            tkplace(frameCol, relx = 0.02, rely = 0.56, relwidth = 0.96, 
                relheight = 0.27, `in` = namingtt)
            tkplace(tklabel(frameCol, text = "Does the data contain a column of object category information?", 
                font = fontsmall), relx = 0.02, rely = 0.25, 
                `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Yes", background = "white"), 
                relx = 0.1, rely = 0.6, `in` = frameCol)
            ColDat <- tclVar(0)
            ColDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ColDat.CB, variable = ColDat)
            tkplace(ColDat.CB, relx = 0.25, rely = 0.6, `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Which Column?", 
                background = "white"), relx = 0.4, rely = 0.6, 
                `in` = frameCol)
            ColColumn <- tclVar("First")
            ColColumn.ComboBox <- tkwidget(namingtt, "ComboBox", 
                editable = FALSE, values = c("First", "Last"), 
                width = 6, textvariable = ColColumn)
            tkplace(ColColumn.ComboBox, relx = 0.7, rely = 0.6, 
                `in` = frameCol)
            OnOk.Name <- function() {
                MGvar$ClasTabCols <<- c()
                NName = as.character(tclvalue(ChangingName))
                MGvar$datatitle <<- NName
                tclvalue(labelText) <- paste("Active Dataset is", 
                  MGvar$datatitle)
                tkentryconfigure(dataMenu, 7, state = "disabled")
                YesCol <- as.character(tclvalue(ColDat))
                if (YesCol == "1") {
                  if (as.character(tclvalue(ColColumn)) == "First") {
                    MGvar$ClasVec <<- MGvar$distmat[, 1]
                    MGvar$activedata <<- MGvar$activedata[, -1]
                    MGvar$distmat <<- MGvar$distmat[, -1]
                  }
                  if (as.character(tclvalue(ColColumn)) == "Last") {
                    MGvar$ClasVec <<- MGvar$distmat[, ncol(MGvar$activedata)]
                    MGvar$activedata <<- MGvar$activedata[, -ncol(MGvar$activedata)]
                    MGvar$distmat <<- MGvar$distmat[, -ncol(MGvar$distmat)]
                  }
                  MGvar$ClasTab <<- table(MGvar$ClasVec)
                  potcols <- c(brewer.pal(9, "Set1"), brewer.pal(12, 
                    "Paired"))
                  MGvar$distmat <<- as.matrix(MGvar$distmat)
                  for (i in 1:nrow(MGvar$distmat)) {
                    for (j in 1:length(MGvar$ClasTab)) {
                      if (MGvar$ClasVec[i] == names(MGvar$ClasTab)[j]) {
                        if (j == length(potcols)) {
                          num = length(potcols)
                        }
                        if (j != length(potcols)) {
                          num = j%%length(potcols)
                        }
                        MGvar$MDSmat.Cols[i] <<- potcols[num]
                      }
                    }
                  }
                  for (i in 1:length(MGvar$ClasTab)) {
                    MGvar$ClasTabCols <<- c(MGvar$ClasTabCols, 
                      potcols[i])
                    names(MGvar$ClasTabCols)[i] <<- names(MGvar$ClasTab)[i]
                  }
                  MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
                  tkentryconfigure(dataMenu, 7, state = "active")
                }
                else {
                  FirstPointColInitialise()
                  PointColInitialise()
                }
                Scl <- as.character(tclvalue(ScDat.val))
                if (Scl == "1") {
                  MGvar$activedata <<- stand.1.range(MGvar$activedata)
                  MGvar$distmat <<- stand.1.range(MGvar$distmat)
                }
                if (is.numeric(as.matrix(MGvar$activedata))) {
                  if (nrow(MGvar$activedata) == ncol(MGvar$activedata)) {
                    if (all(diag(as.matrix(MGvar$activedata)) == 
                      0)) {
                      MGvar$originaldistmat <<- MGvar$distmat
                      MGvar$activedata <<- MGvar$distmat
                      MGvar$distmat.T1 <<- MGvar$distmat
                      MGvar$distmat.T2 <<- MGvar$distmat
                      MGvar$distmat.T3 <<- MGvar$distmat
                      MGvar$distmat.T4 <<- MGvar$distmat
                      MGvar$distmat.T5 <<- MGvar$distmat
                      MGvar$maxdims <<- nrow(MGvar$activedata) - 
                        1
                      MGvar$removedpointsactivedata <<- MGvar$activedata
                      tkentryconfigure(functionMenu, 5, state = "disabled")
                      tkentryconfigure(MainPlotMenu1, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu2, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu3, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu4, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu5, 15, state = "disabled")
                      MGvar$dMeas <<- NA
                      MGvar$dMeas.T1.temp <<- NA
                      MGvar$dMeas.T2.temp <<- NA
                      MGvar$dMeas.T3.temp <<- NA
                      MGvar$dMeas.T4.temp <<- NA
                      MGvar$dMeas.T5.temp <<- NA
                      distmatcheck()
                    }
                    else {
                      tkmessageBox(message = "The Dissimilarity Matrix requires zero's along the diagonal. Please retry data upload", 
                        type = "ok")
                      MGvar$activedata <<- as.matrix(0)
                      MGvar$distmat <<- as.matrix(0)
                      tclvalue(labelText) <- paste("No Active Dataset")
                      MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                      tkentryconfigure(dataMenu, 7, state = "disabled")
                    }
                  }
                  else {
                    tkmessageBox(message = "The Dissimilarity Matrix must be square! Please retry data upload.", 
                      type = "ok")
                    MGvar$activedata <<- as.matrix(0)
                    MGvar$distmat <<- as.matrix(0)
                    tclvalue(labelText) <- paste("No Active Dataset")
                    MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                    tkentryconfigure(dataMenu, 7, state = "disabled")
                  }
                }
                else {
                  tkmessageBox(message = "Non-Numeric values in the data. Please retry data upload.", 
                    type = "ok")
                  MGvar$activedata <<- as.matrix(0)
                  MGvar$distmat <<- as.matrix(0)
                  tclvalue(labelText) <- paste("No Active Dataset")
                  MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                  tkentryconfigure(dataMenu, 7, state = "disabled")
                }
                if (nrow(MGvar$activedata) > 50) {
                  LargeDataOps()
                }
                tkdestroy(namingtt)
            }
            tkplace(tkbutton(namingtt, text = "OK", width = 10, 
                command = function() OnOk.Name()), relx = 0.38, 
                rely = 0.88, `in` = namingtt)
            tkbind(namingtt, "<Return>", OnOk.Name)
            tkentryconfigure(EditDataMenu, 0, state = "disabled")
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 2, state = "disabled")
            tkentryconfigure(EditDataMenu, 3, state = "disabled")
            MGvar$tShepx <<- as.vector(0)
            ClearRemindex()
        }
    }
    LoadSimSet <- function() {
        fileName <- tclvalue(tkgetOpenFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            loadeddata = read.table(fileName)
            MGvar$activedata <<- loadeddata
            MGvar$simmat <<- loadeddata
            namingtt <- tktoplevel()
            tkwm.resizable(namingtt, "0", "0")
            tkwm.deiconify(namingtt)
            tkwm.title(namingtt, "New Active Similarity Matrix Options")
            tkwm.geometry(namingtt, "350x300")
            Loadcanvas = tkcanvas(namingtt, width = "1128", height = "756", 
                bg = col.sec)
            tkplace(Loadcanvas, relx = 0, rely = 0, relwidth = 1, 
                relheight = 1, `in` = namingtt)
            frameNaming <- tkwidget(namingtt, "TitleFrame", text = "Dataset Name", 
                background = "white")
            tkplace(frameNaming, relx = 0.02, rely = 0.02, relwidth = 0.96, 
                relheight = 0.25, `in` = namingtt)
            tkplace(tklabel(frameNaming, text = "Enter the name of your DataSet", 
                background = "white"), relx = 0.08, rely = 0.4, 
                `in` = frameNaming)
            ChangingName = tclVar("")
            NewNameBox = tkentry(namingtt, width = 15, textvariable = ChangingName)
            tkplace(NewNameBox, relx = 0.65, rely = 0.4, `in` = frameNaming)
            fontsmall <- tkfont.create(family = "times", size = 9)
            frameScale <- tkwidget(namingtt, "TitleFrame", text = "Dataset Scale", 
                background = "white")
            tkplace(frameScale, relx = 0.02, rely = 0.29, relwidth = 0.96, 
                relheight = 0.25, `in` = namingtt)
            tkplace(tklabel(frameScale, text = "Scaling Data will scale the columns of your data between 0 and 1.", 
                font = fontsmall), relx = 0.04, rely = 0.25, 
                `in` = frameScale)
            tkplace(tklabel(frameScale, text = "Scale your active data?", 
                background = "white"), relx = 0.08, rely = 0.6, 
                `in` = frameScale)
            ScDat.val <- tclVar(0)
            ScDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ScDat.CB, variable = ScDat.val)
            tkplace(ScDat.CB, relx = 0.75, rely = 0.6, `in` = frameScale)
            frameCol <- tkwidget(namingtt, "TitleFrame", text = "Dataset Colours", 
                background = "white")
            tkplace(frameCol, relx = 0.02, rely = 0.56, relwidth = 0.96, 
                relheight = 0.27, `in` = namingtt)
            tkplace(tklabel(frameCol, text = "Does tdahe data contain a column of object category information?", 
                font = fontsmall), relx = 0.02, rely = 0.25, 
                `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Yes", background = "white"), 
                relx = 0.1, rely = 0.6, `in` = frameCol)
            ColDat <- tclVar(0)
            ColDat.CB <- tk2checkbutton(namingtt)
            tkconfigure(ColDat.CB, variable = ColDat)
            tkplace(ColDat.CB, relx = 0.25, rely = 0.6, `in` = frameCol)
            tkplace(tklabel(frameCol, text = "Which Column?", 
                background = "white"), relx = 0.4, rely = 0.6, 
                `in` = frameCol)
            ColColumn <- tclVar("First")
            ColColumn.ComboBox <- tkwidget(namingtt, "ComboBox", 
                editable = FALSE, values = c("First", "Last"), 
                width = 6, textvariable = ColColumn)
            tkplace(ColColumn.ComboBox, relx = 0.7, rely = 0.6, 
                `in` = frameCol)
            OnOk.Name <- function() {
                MGvar$ClasTabCols <<- c()
                tkentryconfigure(dataMenu, 7, state = "disabled")
                YesCol <- as.character(tclvalue(ColDat))
                if (YesCol == "1") {
                  if (as.character(tclvalue(ColColumn)) == "First") {
                    MGvar$ClasVec <<- MGvar$simmat[, 1]
                    MGvar$activedata <<- MGvar$activedata[, -1]
                    MGvar$simmat <<- MGvar$simmat[, -1]
                  }
                  if (as.character(tclvalue(ColColumn)) == "Last") {
                    MGvar$ClasVec <<- MGvar$simmat[, ncol(MGvar$activedata)]
                    MGvar$activedata <<- MGvar$activedata[, -ncol(MGvar$activedata)]
                    MGvar$simmat <<- MGvar$simmat[, -ncol(MGvar$simmat)]
                  }
                  MGvar$ClasTab <<- table(MGvar$ClasVec)
                  potcols <- c(brewer.pal(9, "Set1"), brewer.pal(12, 
                    "Paired"))
                  MGvar$simmat <<- as.matrix(MGvar$simmat)
                  for (i in 1:nrow(MGvar$simmat)) {
                    for (j in 1:length(MGvar$ClasTab)) {
                      if (MGvar$ClasVec[i] == names(MGvar$ClasTab)[j]) {
                        if (j == length(potcols)) {
                          num = length(potcols)
                        }
                        if (j != length(potcols)) {
                          num = j%%length(potcols)
                        }
                        MGvar$MDSmat.Cols[i] <<- potcols[num]
                      }
                    }
                  }
                  for (i in 1:length(MGvar$ClasTab)) {
                    MGvar$ClasTabCols <<- c(MGvar$ClasTabCols, 
                      potcols[i])
                    names(MGvar$ClasTabCols)[i] <<- names(MGvar$ClasTab)[i]
                  }
                  MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
                  MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
                  tkentryconfigure(dataMenu, 7, state = "active")
                }
                else {
                  FirstPointColInitialise()
                  PointColInitialise()
                }
                NName = as.character(tclvalue(ChangingName))
                MGvar$datatitle <<- NName
                tclvalue(labelText) <- paste("Active Dataset is", 
                  MGvar$datatitle)
                Scl <- as.character(tclvalue(ScDat.val))
                if (Scl == "1") {
                  MGvar$activedata <<- stand.1.range(MGvar$activedata)
                  MGvar$simmat <<- stand.1.range(MGvar$distmat)
                }
                if (is.numeric(as.matrix(MGvar$activedata))) {
                  MGvar$distmat <<- 1 - MGvar$simmat
                  MGvar$distmat <<- as.matrix(MGvar$distmat)
                  diag(MGvar$distmat) <<- 0
                  if (nrow(MGvar$activedata) == ncol(MGvar$activedata)) {
                    if (all(diag(as.matrix(MGvar$activedata)) == 
                      0)) {
                      MGvar$originaldistmat <<- MGvar$distmat
                      MGvar$activedata <<- MGvar$distmat
                      MGvar$distmat.T1 <<- MGvar$distmat
                      MGvar$distmat.T2 <<- MGvar$distmat
                      MGvar$distmat.T3 <<- MGvar$distmat
                      MGvar$distmat.T4 <<- MGvar$distmat
                      MGvar$distmat.T5 <<- MGvar$distmat
                      MGvar$maxdims <<- nrow(MGvar$activedata) - 
                        1
                      MGvar$removedpointsactivedata <<- MGvar$activedata
                      tkentryconfigure(functionMenu, 5, state = "disabled")
                      tkentryconfigure(MainPlotMenu1, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu2, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu3, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu4, 15, state = "disabled")
                      tkentryconfigure(MainPlotMenu5, 15, state = "disabled")
                      MGvar$dMeas <<- NA
                      MGvar$dMeas.T1.temp <<- NA
                      MGvar$dMeas.T2.temp <<- NA
                      MGvar$dMeas.T3.temp <<- NA
                      MGvar$dMeas.T4.temp <<- NA
                      MGvar$dMeas.T5.temp <<- NA
                      distmatcheck()
                    }
                    else {
                      tkmessageBox(message = "The Similarity Matrix requires zero's along the diagonal. Please retry data upload", 
                        type = "ok")
                      MGvar$activedata <<- as.matrix(0)
                      MGvar$distmat <<- as.matrix(0)
                      tclvalue(labelText) <- paste("No Active Dataset")
                      MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                      MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                      tkentryconfigure(dataMenu, 7, state = "disabled")
                    }
                  }
                  else {
                    tkmessageBox(message = "The Similarity Matrix must be square! Please retry data upload.", 
                      type = "ok")
                    MGvar$activedata <<- as.matrix(0)
                    MGvar$distmat <<- as.matrix(0)
                    tclvalue(labelText) <- paste("No Active Dataset")
                    MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                    MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                    tkentryconfigure(dataMenu, 7, state = "disabled")
                  }
                }
                else {
                  tkmessageBox(message = "Non-Numeric values in the data. Please retry data upload.", 
                    type = "ok")
                  MGvar$activedata <<- as.matrix(0)
                  MGvar$distmat <<- as.matrix(0)
                  tclvalue(labelText) <- paste("No Active Dataset")
                  MGvar$MDSmat.Cols.T1 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T2 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T3 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T4 <<- as.vector(0)
                  MGvar$MDSmat.Cols.T5 <<- as.vector(0)
                  tkentryconfigure(dataMenu, 7, state = "disabled")
                }
                if (nrow(MGvar$activedata) > 50) {
                  LargeDataOps()
                }
                tkdestroy(namingtt)
            }
            tkplace(tkbutton(namingtt, text = "OK", width = 10, 
                command = function() OnOk.Name()), relx = 0.38, 
                rely = 0.88, `in` = namingtt)
            tkbind(namingtt, "<Return>", OnOk.Name)
            tkentryconfigure(EditDataMenu, 0, state = "disabled")
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 2, state = "active")
            tkentryconfigure(EditDataMenu, 3, state = "disabled")
            MGvar$tShepx <<- as.vector(0)
            ClearRemindex()
        }
    }
    LoadCorMat <- function() {
        fileName <- tclvalue(tkgetOpenFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            tkmessageBox(message = paste("The file selected was", 
                fileName))
            loadeddata = read.table(fileName)
            MGvar$cormat <<- loadeddata
            namingtt <- tktoplevel()
            tkwm.resizable(namingtt, "0", "0")
            tkwm.deiconify(namingtt)
            tkwm.title(namingtt, "New Active Correlation Matrix")
            tkwm.geometry(namingtt, "400x150")
            Loadcanvas = tkcanvas(namingtt, width = "1128", height = "756", 
                bg = col.sec)
            tkplace(Loadcanvas, relx = 0, rely = 0, relwidth = 1, 
                relheight = 1, `in` = namingtt)
            frameNaming <- tkwidget(namingtt, "TitleFrame", text = "New Correlation Matrix Options", 
                background = "white")
            tkplace(frameNaming, relx = 0.05, rely = 0.02, relwidth = 0.9, 
                relheight = 0.9, `in` = namingtt)
            tkplace(tklabel(frameNaming, text = "Enter the name of your DataSet", 
                background = "white"), relx = 0.05, rely = 0.3, 
                `in` = frameNaming)
            ChangingName = tclVar("")
            NewNameBox = tkentry(namingtt, width = 15, textvariable = ChangingName)
            tkplace(NewNameBox, relx = 0.6, rely = 0.3, `in` = frameNaming)
            OnOk.Name <- function() {
                NName = as.character(tclvalue(ChangingName))
                MGvar$datatitle <<- NName
                tclvalue(labelText) <- paste("Active Dataset is", 
                  MGvar$datatitle)
                tkdestroy(namingtt)
            }
            tkplace(tkbutton(namingtt, text = "OK", width = 10, 
                command = function() OnOk.Name()), relx = 0.4, 
                rely = 0.7, `in` = frameNaming)
            tkbind(namingtt, "<Return>", OnOk.Name)
            tkentryconfigure(EditDataMenu, 3, state = "active")
            MGvar$maxdims <<- nrow(MGvar$cormat) - 1
            MGvar$tShepx <<- as.vector(0)
            ClearRemindex()
        }
    }
    SaveDataSet = function() {
        fileName <- tclvalue(tkgetSaveFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            tkmessageBox(message = paste("The file selected was", 
                fileName))
            write.table(MGvar$activedata, fileName)
        }
    }
    SaveDistmat <- function() {
        fileName <- tclvalue(tkgetSaveFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            tkmessageBox(message = paste("The file selected was", 
                fileName))
            write.table(MGvar$distmat, fileName)
        }
    }
    SaveSimmat <- function() {
        fileName <- tclvalue(tkgetSaveFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            tkmessageBox(message = paste("The file selected was", 
                fileName))
            write.table(MGvar$simmat, fileName)
        }
    }
    SaveCormat <- function() {
        fileName <- tclvalue(tkgetSaveFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            tkmessageBox(message = paste("The file selected was", 
                fileName))
            write.table(MGvar$cormat, fileName)
        }
    }
    SaveMDScoords <- function() {
        fileName <- tclvalue(tkgetSaveFile())
        if (!nchar(fileName)) 
            tkmessageBox(message = "No file was selected!")
        else {
            tkmessageBox(message = paste("The file selected was", 
                fileName))
            write.table(MGvar$MDSmat, fileName)
        }
    }
    ShowActiveData = function() {
        if (nrow(MGvar$activedata) == 1 && ncol(MGvar$activedata) == 
            1) {
            tkmessageBox(message = "No Data! please Load Dataset", 
                icon = "error")
        }
        else {
            tempactivedata <- my.tk2edit(MGvar$activedata)
            for (i in 1:nrow(MGvar$activedata)) {
                for (j in 1:ncol(MGvar$activedata)) {
                  MGvar$activedata[i, j] <<- as.numeric(tempactivedata[i, 
                    j])
                }
            }
            MGvar$maxdims <<- nrow(MGvar$activedata) - 1
        }
    }
    ShowActiveDistmat = function() {
        if (nrow(MGvar$distmat) == 1 && ncol(MGvar$distmat) == 
            1) {
            tkmessageBox(message = "No Distance Matrix! Please Load or Calculate Distance Matrix", 
                icon = "error")
        }
        else {
            tempdistmat <- my.tk2edit(MGvar$distmat)
            for (i in 1:nrow(MGvar$distmat)) {
                for (j in 1:ncol(MGvar$distmat)) {
                  MGvar$distmat[i, j] <<- as.numeric(tempdistmat[i, 
                    j])
                }
            }
            ActivedistMat()
            MGvar$maxdims <<- nrow(MGvar$distmat) - 1
        }
    }
    ShowActiveSimmat = function() {
        if (nrow(MGvar$distmat) == 1 && ncol(MGvar$distmat) == 
            1) {
            tkmessageBox(message = "No Distance Matrix! Please Load or Calculate Distance Matrix", 
                icon = "error")
        }
        else {
            tempsimmat <- my.tk2edit(MGvar$simmat)
            for (i in 1:nrow(MGvar$simmat)) {
                for (j in 1:ncol(MGvar$simmat)) {
                  MGvar$simmat[i, j] <<- as.numeric(tempsimmat[i, 
                    j])
                }
            }
            ActivesimMat()
            MGvar$maxdims <<- nrow(MGvar$simmat) - 1
        }
    }
    ShowActiveCormat = function() {
        if (nrow(MGvar$distmat) == 1 && ncol(MGvar$distmat) == 
            1) {
            tkmessageBox(message = "No Distance Matrix! Please Load or Calculate Distance Matrix", 
                icon = "error")
        }
        else {
            tempcormat <- my.tk2edit(MGvar$cormat)
            for (i in 1:nrow(MGvar$cormat)) {
                for (j in 1:ncol(MGvar$cormat)) {
                  MGvar$cormat[i, j] <<- as.numeric(tempcormat[i, 
                    j])
                }
            }
            ActiveCorMat()
            MGvar$maxdims <<- nrow(MGvar$cormat) - 1
        }
    }
    ShowActiveMDSmat = function() {
        if (nrow(MGvar$MDSmat) == 1 && ncol(MGvar$MDSmat) == 
            1) {
            tkmessageBox(message = "No MDS process has been performed!", 
                icon = "error")
        }
        else {
            MGvar$MDSmat <<- my.tk2edit(MGvar$MDSmat)
            ActiveCoordMat()
            MGvar$maxdims <<- nrow(MGvar$MDSmat) - 1
            tabplot()
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
        }
    }
    Euc.dist = function(data) {
        rows = nrow(data)
        cols = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                dist = 0
                for (k in 1:cols) {
                  dist = dist + (data[i, k] - data[j, k])^2
                }
                MGvar$distmat[i, j] = sqrt(dist)
            }
        }
        return(MGvar$distmat)
    }
    Wtd.Euc <- function(data, weights) {
        rows = nrow(data)
        cols = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                dist = 0
                for (k in 1:cols) {
                  dist = dist + weights[k] * (data[i, k] - data[j, 
                    k])^2
                }
                MGvar$distmat[i, j] = sqrt(dist)
            }
        }
        return(MGvar$distmat)
    }
    Mahalanobis.dist <- function(X, sigma) {
        n <- nrow(X)
        D <- matrix(0, nrow = n, ncol = n)
        for (i in 1:(n - 1)) for (j in (i + 1):n) D[i, j] <- sqrt(t(X[i, 
            ] - X[j, ]) %*% solve(sigma) %*% (X[i, ] - X[j, ]))
        D + t(D)
    }
    CityBlock.dist <- function(data) {
        distm = dist(data, method = "manhattan")
        return(as.matrix(distm))
    }
    Minkowski.dist <- function(data) {
        distm = dist(data, method = "minkowski")
        return(as.matrix(distm))
    }
    Canberra.dist <- function(data) {
        distm = dist(data, method = "canberra")
        return(as.matrix(distm))
    }
    Div.dist <- function(data) {
        rows = nrow(data)
        cols = ncol(data)
        p = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                dist = 0
                for (k in 1:cols) {
                  dist = dist + ((data[i, k] - data[j, k])^2)/((data[i, 
                    k] + data[j, k])^2)
                }
                MGvar$distmat[i, j] = dist/p
            }
        }
        return(MGvar$distmat)
    }
    BC.dist <- function(data) {
        rows = nrow(data)
        cols = ncol(data)
        p = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                top = 0
                bottom = 0
                for (k in 1:cols) {
                  top = top + abs(data[i, k] - data[j, k])
                  bottom = bottom + (data[i, k] + data[j, k])
                }
                MGvar$distmat[i, j] = (1/p) * (top/bottom)
            }
        }
        return(MGvar$distmat)
    }
    Soergel.dist <- function(data) {
        rows = nrow(data)
        cols = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                top = 0
                bottom = 0
                for (k in 1:cols) {
                  top = top + abs(data[i, k] - data[j, k])
                  bottom = bottom + max(data[i, k], data[j, k])
                }
                MGvar$distmat[i, j] = (top/bottom)
            }
        }
        return(MGvar$distmat)
    }
    Bhatt.dist <- function(data) {
        rows = nrow(data)
        cols = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                dist = 0
                for (k in 1:cols) {
                  dist = dist + (sqrt(data[i, k]) - sqrt(data[j, 
                    k]))^2
                }
                MGvar$distmat[i, j] = sqrt(dist)
            }
        }
        return(MGvar$distmat)
    }
    WH.dist <- function(data) {
        rows = nrow(data)
        cols = ncol(data)
        p = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                dist = 0
                for (k in 1:cols) {
                  dist = dist + (1 - (min(data[i, k], data[j, 
                    k])/max(data[i, k], data[j, k])))
                }
                MGvar$distmat[i, j] = (1/p) * dist
            }
        }
        return(MGvar$distmat)
    }
    AngSep.dist <- function(data) {
        rows = nrow(data)
        cols = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                top = 0
                bottom1 = 0
                bottom2 = 0
                for (k in 1:cols) {
                  top = top + data[i, k] * data[j, k]
                  bottom1 = bottom1 + (data[i, k])^2
                  bottom2 = bottom2 + (data[j, k])^2
                }
                MGvar$distmat[i, j] = 1 - (top/sqrt(bottom1 * 
                  bottom2))
            }
        }
        return(MGvar$distmat)
    }
    Corr.dist = function(data) {
        rows = nrow(data)
        cols = ncol(data)
        MGvar$distmat = matrix(0, nrow = rows, ncol = rows)
        for (i in 1:rows) {
            for (j in 1:rows) {
                top = 0
                bottom1 = 0
                bottom2 = 0
                for (k in 1:cols) {
                  top = top + (data[i, k] - mean(data[i, ])) * 
                    (data[j, k] - mean(data[j, ]))
                  bottom1 = bottom1 + (data[i, k] - mean(data[i, 
                    ]))^2
                  bottom2 = bottom2 + (data[j, k] - mean(data[j, 
                    ]))^2
                }
                MGvar$distmat[i, j] = 1 - (top/sqrt(bottom1 * 
                  bottom2))
            }
        }
        return(MGvar$distmat)
    }
    Corr.dist2 = function(data) {
        return(1 - cor(t(data)))
    }
    PlotBlank <- function() {
        zeromat <- as.matrix(0)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkrreplot(img, function() plotting2D(zeromat))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkrreplot(img2, function() plotting2D(zeromat))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkrreplot(img3, function() plotting2D(zeromat))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkrreplot(img4, function() plotting2D(zeromat))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkrreplot(img5, function() plotting2D(zeromat))
        }
    }
    tabplot <- function() {
        if (MGvar$EnActivePlot.switch.T1 == "on" || MGvar$EnActivePlot.switch.T2 == 
            "on" || MGvar$EnActivePlot.switch.T3 == "on" || MGvar$EnActivePlot.switch.T4 == 
            "on" || MGvar$EnActivePlot.switch.T5 == "on") {
            EC.height <- as.numeric(tclvalue(tkwinfo("height", 
                MGcomp$EActive)))
            EC.width <- as.numeric(tclvalue(tkwinfo("width", 
                MGcomp$EActive)))
            WidthScale = 650/1.5
            EShscale <- EC.width/WidthScale
            dimrat = 1.5/1.4
            ESvscale <- EShscale/dimrat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkrreplot(img, function() plotting2D(MGvar$MDSmat.T1, 
                title = MGvar$activeplot.title.T1, Measure = MGvar$dMeas.T1, 
                showtitle = MGvar$activeplot.title.show.T1, showmeas = MGvar$activeplot.distmeas.T1, 
                xlabel = MGvar$activeplot.xlab.T1, ylabel = MGvar$activeplot.ylab.T1, 
                bgcol = MGvar$activeplot.bg.T1, pointcex = MGvar$activeplot.cex.T1, 
                showlabs = MGvar$activeplot.labs.T1, showpoints = MGvar$activeplot.showpoints.T1, 
                pointcol = MGvar$activeplot.pointcol.T1, pointshape = MGvar$activeplot.type.T1, 
                ymeas = MGvar$activeplot.yaxt.T1, xmeas = MGvar$activeplot.xaxt.T1, 
                axcol = MGvar$activeplot.axescol.T1, indexLabeled = MGvar$indexLabeled.T1, 
                zoomedcoords = MGvar$newCoords.T1, showreg = MGvar$activeplot.showreg.T1, 
                regcol = MGvar$activeplot.regcol.T1, showleg = MGvar$activeplot.showleg.T1, 
                Zrat = MGvar$zoominrat.T1, Mvup = MGvar$moveup.T1, 
                Mvdn = MGvar$movedown.T1, Mvlt = MGvar$moveleft.T1, 
                Mvrt = MGvar$moveright.T1, PTcolsindex = MGvar$MDSmat.Cols.T1, 
                showdist = MGvar$activeplot.showdist.T1, distcol = MGvar$activeplot.distcol.T1))
            if (MGvar$EnActivePlot.switch.T1 == "on") {
                tkrreplot(MGcomp$POimg, function() plotting2D(MGvar$MDSmat.T1, 
                  title = MGvar$activeplot.title.T1, Measure = MGvar$dMeas.T1, 
                  showtitle = MGvar$activeplot.title.show.T1, 
                  showmeas = MGvar$activeplot.distmeas.T1, xlabel = MGvar$activeplot.xlab.T1, 
                  ylabel = MGvar$activeplot.ylab.T1, bgcol = MGvar$activeplot.bg.T1, 
                  pointcex = MGvar$activeplot.cex.T1, showlabs = MGvar$activeplot.labs.T1, 
                  showpoints = MGvar$activeplot.showpoints.T1, 
                  pointcol = MGvar$activeplot.pointcol.T1, pointshape = MGvar$activeplot.type.T1, 
                  ymeas = MGvar$activeplot.yaxt.T1, xmeas = MGvar$activeplot.xaxt.T1, 
                  axcol = MGvar$activeplot.axescol.T1, indexLabeled = MGvar$indexLabeled.T1, 
                  zoomedcoords = MGvar$newCoords.T1, showreg = MGvar$activeplot.showreg.T1, 
                  regcol = MGvar$activeplot.regcol.T1, showleg = MGvar$activeplot.showleg.T1, 
                  Zrat = MGvar$zoominrat.T1, Mvup = MGvar$moveup.T1, 
                  Mvdn = MGvar$movedown.T1, Mvlt = MGvar$moveleft.T1, 
                  Mvrt = MGvar$moveright.T1, PTcolsindex = MGvar$MDSmat.Cols.T1, 
                  showdist = MGvar$activeplot.showdist.T1, distcol = MGvar$activeplot.distcol.T1), 
                  hscale = EShscale, vscale = ESvscale)
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkrreplot(img2, function() plotting2D(MGvar$MDSmat.T2, 
                title = MGvar$activeplot.title.T2, Measure = MGvar$dMeas.T2, 
                showtitle = MGvar$activeplot.title.show.T2, showmeas = MGvar$activeplot.distmeas.T2, 
                xlabel = MGvar$activeplot.xlab.T2, ylabel = MGvar$activeplot.ylab.T2, 
                bgcol = MGvar$activeplot.bg.T2, pointcex = MGvar$activeplot.cex.T2, 
                showlabs = MGvar$activeplot.labs.T2, showpoints = MGvar$activeplot.showpoints.T2, 
                pointcol = MGvar$activeplot.pointcol.T2, pointshape = MGvar$activeplot.type.T2, 
                ymeas = MGvar$activeplot.yaxt.T2, xmeas = MGvar$activeplot.xaxt.T2, 
                axcol = MGvar$activeplot.axescol.T2, indexLabeled = MGvar$indexLabeled.T2, 
                zoomedcoords = MGvar$newCoords.T2, showreg = MGvar$activeplot.showreg.T2, 
                regcol = MGvar$activeplot.regcol.T2, showleg = MGvar$activeplot.showleg.T2, 
                Zrat = MGvar$zoominrat.T2, Mvup = MGvar$moveup.T2, 
                Mvdn = MGvar$movedown.T2, Mvlt = MGvar$moveleft.T2, 
                Mvrt = MGvar$moveright.T2, PTcolsindex = MGvar$MDSmat.Cols.T2, 
                showdist = MGvar$activeplot.showdist.T2, distcol = MGvar$activeplot.distcol.T2))
            if (MGvar$EnActivePlot.switch.T2 == "on") {
                tkrreplot(MGcomp$POimg, function() plotting2D(MGvar$MDSmat.T2, 
                  title = MGvar$activeplot.title.T2, Measure = MGvar$dMeas.T2, 
                  showtitle = MGvar$activeplot.title.show.T2, 
                  showmeas = MGvar$activeplot.distmeas.T2, xlabel = MGvar$activeplot.xlab.T2, 
                  ylabel = MGvar$activeplot.ylab.T2, bgcol = MGvar$activeplot.bg.T2, 
                  pointcex = MGvar$activeplot.cex.T2, showlabs = MGvar$activeplot.labs.T2, 
                  showpoints = MGvar$activeplot.showpoints.T2, 
                  pointcol = MGvar$activeplot.pointcol.T2, pointshape = MGvar$activeplot.type.T2, 
                  ymeas = MGvar$activeplot.yaxt.T2, xmeas = MGvar$activeplot.xaxt.T2, 
                  axcol = MGvar$activeplot.axescol.T2, indexLabeled = MGvar$indexLabeled.T2, 
                  zoomedcoords = MGvar$newCoords.T2, showreg = MGvar$activeplot.showreg.T2, 
                  regcol = MGvar$activeplot.regcol.T2, showleg = MGvar$activeplot.showleg.T2, 
                  Zrat = MGvar$zoominrat.T2, Mvup = MGvar$moveup.T2, 
                  Mvdn = MGvar$movedown.T2, Mvlt = MGvar$moveleft.T2, 
                  Mvrt = MGvar$moveright.T2, PTcolsindex = MGvar$MDSmat.Cols.T2, 
                  showdist = MGvar$activeplot.showdist.T2, distcol = MGvar$activeplot.distcol.T2), 
                  hscale = EShscale, vscale = ESvscale)
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkrreplot(img3, function() plotting2D(MGvar$MDSmat.T3, 
                title = MGvar$activeplot.title.T3, Measure = MGvar$dMeas.T3, 
                showtitle = MGvar$activeplot.title.show.T3, showmeas = MGvar$activeplot.distmeas.T3, 
                xlabel = MGvar$activeplot.xlab.T3, ylabel = MGvar$activeplot.ylab.T3, 
                bgcol = MGvar$activeplot.bg.T3, pointcex = MGvar$activeplot.cex.T3, 
                showlabs = MGvar$activeplot.labs.T3, showpoints = MGvar$activeplot.showpoints.T3, 
                pointcol = MGvar$activeplot.pointcol.T3, pointshape = MGvar$activeplot.type.T3, 
                ymeas = MGvar$activeplot.yaxt.T3, xmeas = MGvar$activeplot.xaxt.T3, 
                axcol = MGvar$activeplot.axescol.T3, indexLabeled = MGvar$indexLabeled.T3, 
                zoomedcoords = MGvar$newCoords.T3, showreg = MGvar$activeplot.showreg.T3, 
                regcol = MGvar$activeplot.regcol.T3, showleg = MGvar$activeplot.showleg.T3, 
                Zrat = MGvar$zoominrat.T3, Mvup = MGvar$moveup.T3, 
                Mvdn = MGvar$movedown.T3, Mvlt = MGvar$moveleft.T3, 
                Mvrt = MGvar$moveright.T3, PTcolsindex = MGvar$MDSmat.Cols.T3, 
                showdist = MGvar$activeplot.showdist.T3, distcol = MGvar$activeplot.distcol.T3))
            if (MGvar$EnActivePlot.switch.T3 == "on") {
                tkrreplot(MGcomp$POimg, function() plotting2D(MGvar$MDSmat.T3, 
                  title = MGvar$activeplot.title.T3, Measure = MGvar$dMeas.T3, 
                  showtitle = MGvar$activeplot.title.show.T3, 
                  showmeas = MGvar$activeplot.distmeas.T3, xlabel = MGvar$activeplot.xlab.T3, 
                  ylabel = MGvar$activeplot.ylab.T3, bgcol = MGvar$activeplot.bg.T3, 
                  pointcex = MGvar$activeplot.cex.T3, showlabs = MGvar$activeplot.labs.T3, 
                  showpoints = MGvar$activeplot.showpoints.T3, 
                  pointcol = MGvar$activeplot.pointcol.T3, pointshape = MGvar$activeplot.type.T3, 
                  ymeas = MGvar$activeplot.yaxt.T3, xmeas = MGvar$activeplot.xaxt.T3, 
                  axcol = MGvar$activeplot.axescol.T3, indexLabeled = MGvar$indexLabeled.T3, 
                  zoomedcoords = MGvar$newCoords.T3, showreg = MGvar$activeplot.showreg.T3, 
                  regcol = MGvar$activeplot.regcol.T3, showleg = MGvar$activeplot.showleg.T3, 
                  Zrat = MGvar$zoominrat.T3, Mvup = MGvar$moveup.T3, 
                  Mvdn = MGvar$movedown.T3, Mvlt = MGvar$moveleft.T3, 
                  Mvrt = MGvar$moveright.T3, PTcolsindex = MGvar$MDSmat.Cols.T3, 
                  showdist = MGvar$activeplot.showdist.T3, distcol = MGvar$activeplot.distcol.T3), 
                  hscale = EShscale, vscale = ESvscale)
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkrreplot(img4, function() plotting2D(MGvar$MDSmat.T4, 
                title = MGvar$activeplot.title.T4, Measure = MGvar$dMeas.T4, 
                showtitle = MGvar$activeplot.title.show.T4, showmeas = MGvar$activeplot.distmeas.T4, 
                xlabel = MGvar$activeplot.xlab.T4, ylabel = MGvar$activeplot.ylab.T4, 
                bgcol = MGvar$activeplot.bg.T4, pointcex = MGvar$activeplot.cex.T4, 
                showlabs = MGvar$activeplot.labs.T4, showpoints = MGvar$activeplot.showpoints.T4, 
                pointcol = MGvar$activeplot.pointcol.T4, pointshape = MGvar$activeplot.type.T4, 
                ymeas = MGvar$activeplot.yaxt.T4, xmeas = MGvar$activeplot.xaxt.T4, 
                axcol = MGvar$activeplot.axescol.T4, indexLabeled = MGvar$indexLabeled.T4, 
                zoomedcoords = MGvar$newCoords.T4, showreg = MGvar$activeplot.showreg.T4, 
                regcol = MGvar$activeplot.regcol.T4, showleg = MGvar$activeplot.showleg.T4, 
                Zrat = MGvar$zoominrat.T4, Mvup = MGvar$moveup.T4, 
                Mvdn = MGvar$movedown.T4, Mvlt = MGvar$moveleft.T4, 
                Mvrt = MGvar$moveright.T4, PTcolsindex = MGvar$MDSmat.Cols.T4, 
                showdist = MGvar$activeplot.showdist.T4, distcol = MGvar$activeplot.distcol.T4))
            if (MGvar$EnActivePlot.switch.T4 == "on") {
                tkrreplot(MGcomp$POimg, function() plotting2D(MGvar$MDSmat.T4, 
                  title = MGvar$activeplot.title.T4, Measure = MGvar$dMeas.T4, 
                  showtitle = MGvar$activeplot.title.show.T4, 
                  showmeas = MGvar$activeplot.distmeas.T4, xlabel = MGvar$activeplot.xlab.T4, 
                  ylabel = MGvar$activeplot.ylab.T4, bgcol = MGvar$activeplot.bg.T4, 
                  pointcex = MGvar$activeplot.cex.T4, showlabs = MGvar$activeplot.labs.T4, 
                  showpoints = MGvar$activeplot.showpoints.T4, 
                  pointcol = MGvar$activeplot.pointcol.T4, pointshape = MGvar$activeplot.type.T4, 
                  ymeas = MGvar$activeplot.yaxt.T4, xmeas = MGvar$activeplot.xaxt.T4, 
                  axcol = MGvar$activeplot.axescol.T4, indexLabeled = MGvar$indexLabeled.T4, 
                  zoomedcoords = MGvar$newCoords.T4, showreg = MGvar$activeplot.showreg.T4, 
                  regcol = MGvar$activeplot.regcol.T4, showleg = MGvar$activeplot.showleg.T4, 
                  Zrat = MGvar$zoominrat.T4, Mvup = MGvar$moveup.T4, 
                  Mvdn = MGvar$movedown.T4, Mvlt = MGvar$moveleft.T4, 
                  Mvrt = MGvar$moveright.T4, PTcolsindex = MGvar$MDSmat.Cols.T4, 
                  showdist = MGvar$activeplot.showdist.T4, distcol = MGvar$activeplot.distcol.T4), 
                  hscale = EShscale, vscale = ESvscale)
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkrreplot(img5, function() plotting2D(MGvar$MDSmat.T5, 
                title = MGvar$activeplot.title.T5, Measure = MGvar$dMeas.T5, 
                showtitle = MGvar$activeplot.title.show.T5, showmeas = MGvar$activeplot.distmeas.T5, 
                xlabel = MGvar$activeplot.xlab.T5, ylabel = MGvar$activeplot.ylab.T5, 
                bgcol = MGvar$activeplot.bg.T5, pointcex = MGvar$activeplot.cex.T5, 
                showlabs = MGvar$activeplot.labs.T5, showpoints = MGvar$activeplot.showpoints.T5, 
                pointcol = MGvar$activeplot.pointcol.T5, pointshape = MGvar$activeplot.type.T5, 
                ymeas = MGvar$activeplot.yaxt.T5, xmeas = MGvar$activeplot.xaxt.T5, 
                axcol = MGvar$activeplot.axescol.T5, indexLabeled = MGvar$indexLabeled.T5, 
                zoomedcoords = MGvar$newCoords.T5, showreg = MGvar$activeplot.showreg.T5, 
                regcol = MGvar$activeplot.regcol.T5, showleg = MGvar$activeplot.showleg.T5, 
                Zrat = MGvar$zoominrat.T5, Mvup = MGvar$moveup.T5, 
                Mvdn = MGvar$movedown.T5, Mvlt = MGvar$moveleft.T5, 
                Mvrt = MGvar$moveright.T5, PTcolsindex = MGvar$MDSmat.Cols.T5, 
                showdist = MGvar$activeplot.showdist.T5, distcol = MGvar$activeplot.distcol.T5))
            if (MGvar$EnActivePlot.switch.T5 == "on") {
                tkrreplot(MGcomp$POimg, function() plotting2D(MGvar$MDSmat.T5, 
                  title = MGvar$activeplot.title.T5, Measure = MGvar$dMeas.T5, 
                  showtitle = MGvar$activeplot.title.show.T5, 
                  showmeas = MGvar$activeplot.distmeas.T5, xlabel = MGvar$activeplot.xlab.T5, 
                  ylabel = MGvar$activeplot.ylab.T5, bgcol = MGvar$activeplot.bg.T5, 
                  pointcex = MGvar$activeplot.cex.T5, showlabs = MGvar$activeplot.labs.T5, 
                  showpoints = MGvar$activeplot.showpoints.T5, 
                  pointcol = MGvar$activeplot.pointcol.T5, pointshape = MGvar$activeplot.type.T5, 
                  ymeas = MGvar$activeplot.yaxt.T5, xmeas = MGvar$activeplot.xaxt.T5, 
                  axcol = MGvar$activeplot.axescol.T5, indexLabeled = MGvar$indexLabeled.T5, 
                  zoomedcoords = MGvar$newCoords.T5, showreg = MGvar$activeplot.showreg.T5, 
                  regcol = MGvar$activeplot.regcol.T5, showleg = MGvar$activeplot.showleg.T5, 
                  Zrat = MGvar$zoominrat.T5, Mvup = MGvar$moveup.T5, 
                  Mvdn = MGvar$movedown.T5, Mvlt = MGvar$moveleft.T5, 
                  Mvrt = MGvar$moveright.T5, PTcolsindex = MGvar$MDSmat.Cols.T5, 
                  showdist = MGvar$activeplot.showdist.T5, distcol = MGvar$activeplot.distcol.T5), 
                  hscale = EShscale, vscale = ESvscale)
            }
        }
    }
    ActivedMeas <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$dMeas.T1 <<- MGvar$dMeas
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$dMeas.T2 <<- MGvar$dMeas
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$dMeas.T3 <<- MGvar$dMeas
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$dMeas.T4 <<- MGvar$dMeas
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$dMeas.T5 <<- MGvar$dMeas
        }
    }
    ActiveCoordMat <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDSmat.T1 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDSmat.T2 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDSmat.T3 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDSmat.T4 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDSmat.T5 <<- MGvar$MDSmat
        }
    }
    ActivedistMat <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$distmat.T1 <<- MGvar$distmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$distmat.T2 <<- MGvar$distmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$distmat.T3 <<- MGvar$distmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$distmat.T4 <<- MGvar$distmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$distmat.T5 <<- MGvar$distmat
        }
    }
    ActivesimMat <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$simmat.T1 <<- MGvar$simmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$simmat.T2 <<- MGvar$simmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$simmat.T3 <<- MGvar$simmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$simmat.T4 <<- MGvar$simmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$simmat.T5 <<- MGvar$simmat
        }
    }
    ActiveCorMat <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$cormat.T1 <<- MGvar$cormat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$cormat.T2 <<- MGvar$cormat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$cormat.T3 <<- MGvar$cormat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$cormat.T4 <<- MGvar$cormat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$cormat.T5 <<- MGvar$cormat
        }
    }
    BlankShep <- function() {
        AShepdistmat <<- as.matrix(0)
        AShepMDSmat <<- as.matrix(0)
        blankdist <- as.matrix(0)
        blankconf <- as.matrix(0)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            AShepdistmat.T1 <<- as.matrix(0)
            AShepMDSmat.T1 <<- as.matrix(0)
            tkrreplot(imgshep, function() plotShepard(blankdist, 
                blankconf, Tab = "Tab1"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            AShepdistmat.T2 <<- as.matrix(0)
            AShepMDSmat.T2 <<- as.matrix(0)
            tkrreplot(imgshep, function() plotShepard(blankdist, 
                blankconf, Tab = "Tab2"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            AShepdistmat.T3 <<- as.matrix(0)
            AShepMDSmat.T3 <<- as.matrix(0)
            tkrreplot(imgshep, function() plotShepard(blankdist, 
                blankconf, Tab = "Tab3"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            AShepdistmat.T4 <<- as.matrix(0)
            AShepMDSmat.T4 <<- as.matrix(0)
            tkrreplot(imgshep, function() plotShepard(blankdist, 
                blankconf, Tab = "Tab4"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            AShepdistmat.T5 <<- as.matrix(0)
            AShepMDSmat.T5 <<- as.matrix(0)
            tkrreplot(imgshep, function() plotShepard(blankdist, 
                blankconf, Tab = "Tab5"))
        }
    }
    ActiveShep <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            AShepdistmat.T1 <<- MGvar$distmat
            AShepMDSmat.T1 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            AShepdistmat.T2 <<- MGvar$distmat
            AShepMDSmat.T2 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            AShepdistmat.T3 <<- MGvar$distmat
            AShepMDSmat.T3 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            AShepdistmat.T4 <<- MGvar$distmat
            AShepMDSmat.T4 <<- MGvar$MDSmat
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            AShepdistmat.T5 <<- MGvar$distmat
            AShepMDSmat.T5 <<- MGvar$MDSmat
        }
    }
    BlankScree <- function() {
        blankstress <- as.vector(0)
        MGvar$scree.stress <<- as.vector(0)
        MGvar$screepoints.current <<- as.vector(0)
        MGvar$screepoints.best <<- as.vector(0)
        MGvar$Opt.dim <<- 1
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$scree.stress.T1 <<- as.vector(0)
            MGvar$screepoints.current.T1 <<- as.vector(0)
            MGvar$screepoints.best.T1 <<- as.vector(0)
            MGvar$Opt.dim.T1 <<- 1
            tkrreplot(imgscree, function() plotScree(blankstress, 
                Tab = "Tab1"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$scree.stress.T2 <<- as.vector(0)
            MGvar$screepoints.current.T2 <<- as.vector(0)
            MGvar$screepoints.best.T2 <<- as.vector(0)
            MGvar$Opt.dim.T2 <<- 1
            tkrreplot(imgscree, function() plotScree(blankstress, 
                Tab = "Tab2"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$scree.stress.T3 <<- as.vector(0)
            MGvar$screepoints.current.T3 <<- as.vector(0)
            MGvar$screepoints.best.T3 <<- as.vector(0)
            MGvar$Opt.dim.T3 <<- 1
            tkrreplot(imgscree, function() plotScree(blankstress, 
                Tab = "Tab3"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$scree.stress.T4 <<- as.vector(0)
            MGvar$screepoints.current.T4 <<- as.vector(0)
            MGvar$screepoints.best.T4 <<- as.vector(0)
            MGvar$Opt.dim.T4 <<- 1
            tkrreplot(imgscree, function() plotScree(blankstress, 
                Tab = "Tab4"))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$scree.stress.T5 <<- as.vector(0)
            MGvar$screepoints.current.T5 <<- as.vector(0)
            MGvar$screepoints.best.T5 <<- as.vector(0)
            MGvar$Opt.dim.T4 <<- 1
            tkrreplot(imgscree, function() plotScree(blankstress, 
                Tab = "Tab5"))
        }
    }
    ActiveScree <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$scree.stress.T1 <<- MGvar$scree.stress
            MGvar$screepoints.current.T1 <<- MGvar$screepoints.current
            MGvar$screepoints.best.T1 <<- MGvar$screepoints.best
            MGvar$Opt.dim.T1 <<- MGvar$Opt.dim
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$scree.stress.T2 <<- MGvar$scree.stress
            MGvar$screepoints.current.T2 <<- MGvar$screepoints.current
            MGvar$screepoints.best.T2 <<- MGvar$screepoints.best
            MGvar$Opt.dim.T2 <<- MGvar$Opt.dim
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$scree.stress.T3 <<- MGvar$scree.stress
            MGvar$screepoints.current.T3 <<- MGvar$screepoints.current
            MGvar$screepoints.best.T3 <<- MGvar$screepoints.best
            MGvar$Opt.dim.T3 <<- MGvar$Opt.dim
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$scree.stress.T4 <<- MGvar$scree.stress
            MGvar$screepoints.current.T4 <<- MGvar$screepoints.current
            MGvar$screepoints.best.T4 <<- MGvar$screepoints.best
            MGvar$Opt.dim.T4 <<- MGvar$Opt.dim
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$scree.stress.T5 <<- MGvar$scree.stress
            MGvar$screepoints.current.T5 <<- MGvar$screepoints.current
            MGvar$screepoints.best.T5 <<- MGvar$screepoints.best
            MGvar$Opt.dim.T5 <<- MGvar$Opt.dim
        }
    }
    ActiveType <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDStype.T1 <<- MGvar$MDStype
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDStype.T2 <<- MGvar$MDStype
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDStype.T3 <<- MGvar$MDStype
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDStype.T4 <<- MGvar$MDStype
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDStype.T5 <<- MGvar$MDStype
        }
    }
    ActiveTitle <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$activeplot.title.T1 <<- MGvar$activeplot.title
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$activeplot.title.T2 <<- MGvar$activeplot.title
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$activeplot.title.T3 <<- MGvar$activeplot.title
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$activeplot.title.T4 <<- MGvar$activeplot.title
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$activeplot.title.T5 <<- MGvar$activeplot.title
        }
    }
    ActiveZoomSwitch <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$Tab1.zoomedswitch <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$Tab2.zoomedswitch <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$Tab3.zoomedswitch <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$Tab4.zoomedswitch <<- "off"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$Tab5.zoomedswitch <<- "off"
        }
    }
    ActiveWatchCursor <- function() {
        tkconfigure(mytt, cursor = "watch")
        tkconfigure(myPlottingNB, cursor = "watch")
        tkconfigure(img, cursor = "watch")
        tkconfigure(img2, cursor = "watch")
        tkconfigure(img3, cursor = "watch")
        tkconfigure(img4, cursor = "watch")
        tkconfigure(img5, cursor = "watch")
        tkconfigure(imgscree, cursor = "watch")
        tkconfigure(imgshep, cursor = "watch")
    }
    ActiveArrowCursor <- function() {
        tkconfigure(mytt, cursor = "arrow")
        tkconfigure(myPlottingNB, cursor = "hand2")
        tkconfigure(img, cursor = "arrow")
        tkconfigure(img2, cursor = "arrow")
        tkconfigure(img3, cursor = "arrow")
        tkconfigure(img4, cursor = "arrow")
        tkconfigure(img5, cursor = "arrow")
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkconfigure(img, cursor = "hand2")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkconfigure(img2, cursor = "hand2")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkconfigure(img3, cursor = "hand2")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkconfigure(img4, cursor = "hand2")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkconfigure(img5, cursor = "hand2")
        }
        tkconfigure(imgscree, cursor = "arrow")
        tkconfigure(imgshep, cursor = "hand2")
    }
    ActiveStress <- function() {
        if (MGvar$MDStype != "INDSCAL" && MGvar$MDStype != "Gifi") {
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$MDSStress.T1 <<- round(MGvar$MDSStress, 
                  3)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$MDSStress.T2 <<- round(MGvar$MDSStress, 
                  3)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$MDSStress.T3 <<- round(MGvar$MDSStress, 
                  3)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$MDSStress.T4 <<- round(MGvar$MDSStress, 
                  3)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$MDSStress.T5 <<- round(MGvar$MDSStress, 
                  3)
            }
        }
        else {
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$MDSStress.T1 <<- paste("*")
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$MDSStress.T2 <<- paste("*")
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$MDSStress.T3 <<- paste("*")
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$MDSStress.T4 <<- paste("*")
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$MDSStress.T5 <<- paste("*")
            }
        }
    }
    ActiveTab <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tk2notetab.select(myPlottingNB, "Plot1")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tk2notetab.select(myPlottingNB, "Plot2")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tk2notetab.select(myPlottingNB, "Plot3")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tk2notetab.select(myPlottingNB, "Plot4")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tk2notetab.select(myPlottingNB, "Plot5")
        }
    }
    ActiveIndex <- function() {
        if (MGvar$GenSet.ClearIL == "yes") {
            MGvar$indexLabeled <<- c()
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$indexLabeled.T1 <<- c()
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$indexLabeled.T2 <<- c()
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$indexLabeled.T3 <<- c()
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$indexLabeled.T4 <<- c()
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$indexLabeled.T5 <<- c()
            }
        }
        if (MGvar$GenSet.KeepShep == "no") {
            MGvar$Shep.indexLabeled <<- c()
        }
    }
    ActiveLabs <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$activeplot.xlab.T1 <<- MGvar$activeplot.xlab
            MGvar$activeplot.ylab.T1 <<- MGvar$activeplot.ylab
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$activeplot.xlab.T2 <<- MGvar$activeplot.xlab
            MGvar$activeplot.ylab.T2 <<- MGvar$activeplot.ylab
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$activeplot.xlab.T3 <<- MGvar$activeplot.xlab
            MGvar$activeplot.ylab.T3 <<- MGvar$activeplot.ylab
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$activeplot.xlab.T4 <<- MGvar$activeplot.xlab
            MGvar$activeplot.ylab.T4 <<- MGvar$activeplot.ylab
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$activeplot.xlab.T5 <<- MGvar$activeplot.xlab
            MGvar$activeplot.ylab.T5 <<- MGvar$activeplot.ylab
        }
    }
    ActiveTabDims <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$TabDims.T1 <<- paste(MGvar$PlottingDimX, "&", 
                MGvar$PlottingDimY, sep = "")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$TabDims.T2 <<- paste(MGvar$PlottingDimX, "&", 
                MGvar$PlottingDimY, sep = "")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$TabDims.T3 <<- paste(MGvar$PlottingDimX, "&", 
                MGvar$PlottingDimY, sep = "")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$TabDims.T4 <<- paste(MGvar$PlottingDimX, "&", 
                MGvar$PlottingDimY, sep = "")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$TabDims.T5 <<- paste(MGvar$PlottingDimX, "&", 
                MGvar$PlottingDimY, sep = "")
        }
    }
    ActiveMDSdims <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDS.dimensions.T1 <<- MGvar$MDS.dimensions
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDS.dimensions.T2 <<- MGvar$MDS.dimensions
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDS.dimensions.T3 <<- MGvar$MDS.dimensions
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDS.dimensions.T4 <<- MGvar$MDS.dimensions
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDS.dimensions.T5 <<- MGvar$MDS.dimensions
        }
    }
    ActivePos <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$zoominrat.T1 <<- 1
            MGvar$moveup.T1 <<- 1
            MGvar$movedown.T1 <<- 1
            MGvar$moveleft.T1 <<- 1
            MGvar$moveright.T1 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$zoominrat.T2 <<- 1
            MGvar$moveup.T2 <<- 1
            MGvar$movedown.T2 <<- 1
            MGvar$moveleft.T2 <<- 1
            MGvar$moveright.T2 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$zoominrat.T3 <<- 1
            MGvar$moveup.T3 <<- 1
            MGvar$movedown.T3 <<- 1
            MGvar$moveleft.T3 <<- 1
            MGvar$moveright.T3 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$zoominrat.T4 <<- 1
            MGvar$moveup.T4 <<- 1
            MGvar$movedown.T4 <<- 1
            MGvar$moveleft.T4 <<- 1
            MGvar$moveright.T4 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$zoominrat.T5 <<- 1
            MGvar$moveup.T5 <<- 1
            MGvar$movedown.T5 <<- 1
            MGvar$moveleft.T5 <<- 1
            MGvar$moveright.T5 <<- 1
        }
    }
    ActiveIterTab <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDS.iter.T1 <<- length(MGvar$stressitervec)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDS.iter.T2 <<- length(MGvar$stressitervec)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDS.iter.T3 <<- length(MGvar$stressitervec)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDS.iter.T4 <<- length(MGvar$stressitervec)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDS.iter.T5 <<- length(MGvar$stressitervec)
        }
    }
    ActiveTolTab <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDS.tol.T1 <<- MGvar$MDS.tol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDS.tol.T2 <<- MGvar$MDS.tol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDS.tol.T3 <<- MGvar$MDS.tol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDS.tol.T4 <<- MGvar$MDS.tol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDS.tol.T5 <<- MGvar$MDS.tol
        }
    }
    ActiveDistFunc <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$DistFunc.T1 <<- MGvar$DistFunc
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$DistFunc.T2 <<- MGvar$DistFunc
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$DistFunc.T3 <<- MGvar$DistFunc
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$DistFunc.T4 <<- MGvar$DistFunc
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$DistFunc.T5 <<- MGvar$DistFunc
        }
    }
    ActiveIsMetric <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$is.Metric.T1 <<- MGvar$is.Metric
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$is.Metric.T2 <<- MGvar$is.Metric
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$is.Metric.T3 <<- MGvar$is.Metric
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$is.Metric.T4 <<- MGvar$is.Metric
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$is.Metric.T5 <<- MGvar$is.Metric
        }
    }
    ClasScal = function(data = MGvar$activedata, fromdistmat = FALSE) {
        DeactivateAll()
        if (nrow(data) == 1 && ncol(data) == 1) {
            tkmessageBox(message = "Invalid Data!", icon = "error")
        }
        else {
            MGvar$is.Metric <<- TRUE
            ActiveIsMetric()
            MGvar$distmat <<- MGvar$originaldistmat
            ActivedistMat()
            ActiveWatchCursor()
            distnameassign()
            ClearRemindex()
            MGvar$tempdims <<- MGvar$MDS.dimensions
            ActivePos()
            MGvar$MDStype <<- "ClasScal"
            CS = cmdscale(MGvar$distmat, k = MGvar$MDS.dimensions)
            MGvar$MDSStress <<- StressCalculation(CS, MGvar$distmat, 
                isMet = TRUE)
            if (MGvar$MDS.dimensions == 1) {
                UniCoords = matrix(0, nrow = nrow(MGvar$distmat), 
                  ncol = 2)
                for (i in 1:nrow(UniCoords)) {
                  UniCoords[i, ] = c(CS[i], 0)
                }
                MGvar$MDSmat <<- UniCoords
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$activeplot.title <<- paste("Unidimensional Classic Scaling on", 
                  MGvar$datatitle, "data")
            }
            if (MGvar$MDS.dimensions == 2) {
                MGvar$PlottingDimX <<- 1
                MGvar$PlottingDimY <<- 2
                MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                MGvar$MDSmat <<- CS
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$activeplot.title <<- paste("Classic Scaling on", 
                  MGvar$datatitle, "data")
            }
            if (MGvar$MDS.dimensions > 3) {
                MGvar$MDSmat.LD <<- CS
                LargeDimensions()
                MGvar$activeplot.title <<- paste("Classic Scaling on", 
                  MGvar$datatitle, "data: Dim", MGvar$PlottingDimX, 
                  "vs. Dim", MGvar$PlottingDimY)
                MGvar$sthreeDplot.title <<- paste("3D Classic Scaling on", 
                  MGvar$datatitle, "data")
                MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                ActiveLabs()
            }
            if (MGvar$MDS.dimensions == 3) {
                MGvar$MDSmat <<- CS
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$sthreeDplot.title <<- paste("3D Classic Scaling on", 
                  MGvar$datatitle, "data")
                threedimplotting()
            }
            if (MGvar$tempdims != 3) {
                ActivedMeas()
                ActiveMDSdims()
                ActiveTab()
                ActiveStress()
                ActiveTabDims()
                ActiveCoordMat()
                ActiveTitle()
                ActiveIndex()
                ActiveType()
                ActiveZoomSwitch()
                MGvar$activeplot.cex <<- 0.6
                tableupdate()
                if (MGvar$GenSet.ClearCols == "yes") {
                  PointColInitialise()
                }
                if (MGvar$GenSet.CalcMDS == "yes") {
                  tabplot()
                }
                else {
                  PlotBlank()
                }
            }
            tkfocus(mytt)
            if (MGvar$GenSet.CalcShep == "yes") {
                MGvar$shep.firstrun <<- "yes"
                AShepdistmat <<- MGvar$distmat
                AShepMDSmat <<- MGvar$MDSmat
                ActiveShep()
                tkrreplot(imgshep)
            }
            else {
                BlankShep()
            }
            tkwm.state(mytt, "normal")
            tkwm.focusmodel(mytt)
            if (MGvar$GenSet.CalcScree == "yes") {
                ndim = MGvar$MDS.dimensions
                bestchange = 0
                cap = min(MGvar$maxdims, 12)
                for (i in 1:cap) {
                  nextCS <- cmdscale(MGvar$distmat, k = i)
                  MGvar$scree.stress[i] <<- StressCalculation(nextCS, 
                    MGvar$distmat, isMet = TRUE)
                }
                for (i in 2:(cap - 2)) {
                  m1 = -(MGvar$scree.stress[i - 1] - MGvar$scree.stress[i])
                  m2 = -(MGvar$scree.stress[i] - MGvar$scree.stress[i + 
                    1])
                  tanthet = abs((m1 - m2)/(1 + m1 * m2))
                  newthet = atan(tanthet)
                  ultchange = MGvar$scree.stress[i - 1] - MGvar$scree.stress[i + 
                    1]
                  if (newthet > bestchange && i < 7) {
                    bestchange = newthet
                    MGvar$Opt.dim <<- i
                  }
                }
                for (i in 1:(cap - 2)) {
                  if (MGvar$Opt.dim == i) {
                    MGvar$screepoints.best[i] <<- MGvar$scree.stress[i]
                  }
                  else {
                    MGvar$screepoints.best[i] <<- ""
                  }
                  if (ndim == i) {
                    MGvar$screepoints.current[i] <<- MGvar$scree.stress[i]
                  }
                  else {
                    MGvar$screepoints.current[i] <<- ""
                  }
                }
                ActiveScree()
                tkrreplot(imgscree)
            }
            else {
                BlankScree()
            }
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 4, state = "active")
        }
        ActivateAll()
        ResetUndo()
    }
    smacof.UCT <- function(DeltaMat, Wmat, X.ini, tol = MGvar$MDS.tol, 
        iter = MGvar$MDS.iter.max) {
        start.time <- proc.time()
        n <- nrow(DeltaMat)
        termin <- "Maximum iterations not reached"
        if (is.null(Wmat)) {
            Wmat <- matrix(1, nrow = n, ncol = n)
            Wmat[lower.tri(Wmat)] <- 0
        }
        if (is.null(rownames(DeltaMat))) 
            rownames(DeltaMat) <- 1:n
        eta.sq.delta <- sum(DeltaMat * DeltaMat * Wmat)
        X <- X.ini
        V <- -(Wmat + t(Wmat)) + diag(apply((Wmat + t(Wmat)), 
            1, sum))
        Vplus <- solve(V + matrix(1, nrow = n, ncol = n))
        W.delta.mat <- DeltaMat * Wmat
        itnum <- 0
        stressval <- 0
        tel <- 0
        stress.begin <- Inf
        if (!MGvar$fromscree) {
            MGvar$stressitervec <<- as.vector(0)
        }
        if (MGvar$GenSet.UpConf == "yes") {
            MGvar$MDSmat <<- X
            rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
            ActiveCoordMat()
            tabplot()
        }
        repeat {
            if (MGvar$endprocess == "yes") {
                break
            }
            tel <- tel + 1
            if (tel > MGvar$MDS.iter.max) {
                termin <- "Maximum number iterations reached. Increase maximum iterations!"
                tkmessageBox(message = termin)
                MGvar$breakscree <<- TRUE
                break
            }
            D.z <- Euclid.dist(X)
            B.z.mat <- ifelse(upper.tri(W.delta.mat), W.delta.mat/D.z, 
                W.delta.mat)
            B.z.mat <- B.z.mat + t(B.z.mat)
            B.z.mat <- -B.z.mat + diag(apply(B.z.mat, 2, sum))
            eta.sq.x <- sum(diag(t(X) %*% V %*% X))
            rho.x <- sum(diag(t(X) %*% B.z.mat %*% X))
            sigma.r <- eta.sq.delta + eta.sq.x - 2 * rho.x
            X.up <- Vplus %*% B.z.mat %*% X
            X <- X.up
            if (MGvar$GenSet.UpConf == "yes") {
                MGvar$MDSmat <<- X
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                ActiveCoordMat()
                if (MGvar$MDS.dimensions == 2) {
                  tabplot()
                }
                if (MGvar$MDS.dimensions == 3) {
                  tkrreplot(img3Dstat)
                }
            }
            if (MGvar$GenSet.UpShep == "yes" && MGvar$GenSet.CalcShep == 
                "yes") {
                MGvar$MDSmat <<- X
                if (tk2notetab.text(mySecondaryNB) == "Shepard") {
                  tkrreplot(imgshep)
                }
                if (MGvar$EnShep.switch == "on") {
                  Replot.imgEShep()
                }
            }
            if (MGvar$GenSet.UpStress == "yes" && MGvar$GenSet.CalcStress == 
                "yes") {
                iter.time <- proc.time()
                MGvar$proctime <<- (iter.time - start.time)[3]
                if (tk2notetab.text(mySecondaryNB) == "Stress Plot") {
                  tkrreplot(imgstress)
                }
                if (tk2notetab.text(mySecondaryNB) == "Log Stress") {
                  tkrreplot(imgstress2)
                }
                if (MGvar$EnStress2.switch == "on") {
                  tkrreplot(MGcomp$imgEStress2, hscale = MGvar$ESt2hscale, 
                    vscale = MGvar$ESt2vscale)
                }
                if (MGvar$EnStress.switch == "on") {
                  tkrreplot(MGcomp$imgEStress, hscale = MGvar$ESt1hscale, 
                    vscale = MGvar$ESt1vscale)
                }
            }
            if ((stress.begin - sigma.r) < MGvar$MDS.tol) 
                break
            if (MGvar$GenSet.UpProg == "yes") {
                tkconfigure(MDSGUI.ProgressBar, value = (tel/MGvar$MDS.iter.max) * 
                  100)
            }
            stress.begin <- sigma.r
            itnum <- append(itnum, tel)
            MGvar$stressval <<- append(MGvar$stressval, sigma.r)
            if (!MGvar$fromscree) {
                MGvar$stressitervec[tel] <<- StressCalculation(MGvar$MDSmat, 
                  MGvar$distmat, isMet = TRUE)
            }
        }
        if (MGvar$endprocess == "no") {
            list(conf = X.up, stress = sigma.r, termin = termin)
        }
    }
    Euclid.dist <- function(X) {
        n <- nrow(X)
        B.mat <- X %*% t(X)
        d.squared <- matrix(diag(B.mat), nrow = n, ncol = n) + 
            t(matrix(diag(B.mat), nrow = n, ncol = n)) - 2 * 
            B.mat
        sqrt(d.squared)
    }
    stand.1.range <- function(mat) {
        mins <- apply(mat, 2, min)
        maks <- apply(mat, 2, max)
        sweep(sweep(mat, 2, mins, "-"), 2, maks - mins, "/")
    }
    datsetup.smacof <- function(datmat, iter = 1000, dim = 2) {
        datmat = as.matrix(datmat)
        datmat <- stand.1.range(datmat)
        delta.mat <- Euclid.dist(datmat)
        x.ini <- cmdscale(delta.mat, dim)
        x.ini <- stand.1.range(matrix(x.ini, ncol = dim))
        smacof.UCT(DeltaMat = delta.mat, Wmat = NULL, X.ini = x.ini, 
            iter = iter)
    }
    MyMetricSmacofSym <- function(data = MGvar$activedata, fromdistmat = FALSE, 
        initial = TRUE) {
        DeactivateAll()
        if (nrow(data) == 1 && ncol(data) == 1) {
            tkmessageBox(message = "Invalid Data!", icon = "error")
        }
        else {
            MGvar$is.Metric <<- TRUE
            ActiveIsMetric()
            if (initial) {
                ClearRemindex()
                MGvar$distmat <<- MGvar$originaldistmat
                ActivedistMat()
            }
            distnameassign()
            if (tclvalue(MGvar$MDSops.startconfig) == "ClasScal") {
                startconfig = cmdscale(MGvar$distmat, k = MGvar$MDS.dimensions)
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "RandConf") {
                startconfig = runif(MGvar$MDS.dimensions * nrow(data))
                startconfig <- stand.1.range(matrix(startconfig, 
                  ncol = MGvar$MDS.dimensions))
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "ExistConf") {
                if (MGvar$startconfig.UseEC.plot == "Plot1") {
                  startconfig = MGvar$MDSmat.T1
                }
                if (MGvar$startconfig.UseEC.plot == "Plot2") {
                  startconfig = MGvar$MDSmat.T2
                }
                if (MGvar$startconfig.UseEC.plot == "Plot3") {
                  startconfig = MGvar$MDSmat.T3
                }
                if (MGvar$startconfig.UseEC.plot == "Plot4") {
                  startconfig = MGvar$MDSmat.T4
                }
                if (MGvar$startconfig.UseEC.plot == "Plot5") {
                  startconfig = MGvar$MDSmat.T5
                }
            }
            if (MGvar$CSinitConf == "no") {
                startconfig = MGvar$MDSmat
                MGvar$CSinitConf <<- "yes"
            }
            startconfig <- stand.1.range(matrix(startconfig, 
                ncol = MGvar$MDS.dimensions))
            ActivePos()
            MGvar$tempdims <<- MGvar$MDS.dimensions
            ActiveWatchCursor()
            MGvar$MDStype <<- "M.Smac.Sym"
            SymMetSmac = smacof.UCT(DeltaMat = MGvar$distmat, 
                Wmat = NULL, X.ini = startconfig, iter = MGvar$MDS.iter.max)
            if (MGvar$endprocess == "no") {
                MGvar$MDSStress <<- StressCalculation(SymMetSmac$conf, 
                  MGvar$distmat, isMet = TRUE)
                if (MGvar$MDS.dimensions == 1) {
                  UniCoords = matrix(0, nrow = nrow(MGvar$distmat), 
                    ncol = 2)
                  for (i in 1:nrow(UniCoords)) {
                    UniCoords[i, ] = c(SymMetSmac$conf[i], 0)
                  }
                  MGvar$MDSmat <<- UniCoords
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$activeplot.title <<- paste("Unidimensional Metric SMACOF on", 
                    MGvar$datatitle, "data")
                }
                if (MGvar$MDS.dimensions == 2) {
                  MGvar$MDSmat <<- SymMetSmac$conf
                  MGvar$PlottingDimX <<- 1
                  MGvar$PlottingDimY <<- 2
                  MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                  MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$activeplot.title <<- paste("Metric SMACOF on", 
                    MGvar$datatitle, "data")
                }
                if (MGvar$MDS.dimensions > 3) {
                  MGvar$MDSmat.LD <<- SymMetSmac$conf
                  LargeDimensions()
                  MGvar$activeplot.title <<- paste("Metric SMACOF on", 
                    MGvar$datatitle, "data: Dim", MGvar$PlottingDimX, 
                    "vs. Dim", MGvar$PlottingDimY)
                  MGvar$sthreeDplot.title <<- paste("3D Metric SMACOF on", 
                    MGvar$datatitle, "data")
                  MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                  MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                  ActiveLabs()
                }
                if (MGvar$MDS.dimensions == 3) {
                  MGvar$MDSmat <<- SymMetSmac$conf
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$sthreeDplot.title <<- paste("3D Metric SMACOF on", 
                    MGvar$datatitle, "data")
                  threedimplotting()
                }
                if (MGvar$tempdims != 3) {
                  ActiveMDSdims()
                  ActivedMeas()
                  ActiveTab()
                  ActiveStress()
                  ActiveTabDims()
                  ActiveCoordMat()
                  ActiveTitle()
                  ActiveIndex()
                  ActiveType()
                  ActiveZoomSwitch()
                  ActiveIterTab()
                  ActiveTolTab()
                  MGvar$activeplot.cex <<- 0.6
                  tableupdate()
                  if (MGvar$GenSet.ClearCols == "yes") {
                    PointColInitialise()
                  }
                  if (MGvar$GenSet.CalcMDS == "yes") {
                    tabplot()
                  }
                  else {
                    PlotBlank()
                  }
                }
                if (MGvar$GenSet.CalcShep == "yes") {
                  MGvar$shep.firstrun <<- "yes"
                  AShepdistmat <<- MGvar$distmat
                  AShepMDSmat <<- MGvar$MDSmat
                  ActiveShep()
                  tkrreplot(imgshep)
                }
                else {
                  BlankShep()
                }
                tkfocus(mytt)
                if (MGvar$GenSet.CalcScree == "yes") {
                  ndim = MGvar$MDS.dimensions
                  bestchange = 0
                  cap = min(MGvar$maxdims, 12)
                  MGvar$fromscree <<- TRUE
                  if (MGvar$GenSet.UpShep == "yes") {
                    MGvar$GenSet.UpShep <<- "no"
                    sheptracker = 1
                  }
                  else {
                    sheptracker = 0
                  }
                  if (MGvar$GenSet.UpConf == "yes") {
                    MGvar$GenSet.UpConf <<- "no"
                    conftracker = 1
                  }
                  else {
                    conftracker = 0
                  }
                  if (MGvar$GenSet.UpStress == "yes") {
                    MGvar$GenSet.UpStress <<- "no"
                    stresstracker = 1
                  }
                  else {
                    stresstracker = 0
                  }
                  for (i in 1:cap) {
                    if (MGvar$endprocess == "yes") {
                      break
                      MGvar$breakscree <<- TRUE
                    }
                    if (MGvar$breakscree) {
                      break
                    }
                    newSmac <- smacof.UCT(DeltaMat = MGvar$distmat, 
                      Wmat = NULL, X.ini = cmdscale(MGvar$distmat, 
                        k = i), iter = MGvar$MDS.iter.max)$conf
                    MGvar$scree.stress[i] <<- StressCalculation(newSmac, 
                      MGvar$distmat, isMet = TRUE)
                  }
                  if (sheptracker == 1) {
                    MGvar$GenSet.UpShep <<- "yes"
                    sheptracker = 0
                  }
                  if (conftracker == 1) {
                    MGvar$GenSet.UpConf <<- "yes"
                    conftracker = 0
                  }
                  if (stresstracker == 1) {
                    MGvar$GenSet.UpStress <<- "yes"
                    stresstracker = 0
                  }
                  if (!MGvar$breakscree) {
                    for (i in 2:(cap - 2)) {
                      m1 = -(MGvar$scree.stress[i - 1] - MGvar$scree.stress[i])
                      m2 = -(MGvar$scree.stress[i] - MGvar$scree.stress[i + 
                        1])
                      tanthet = abs((m1 - m2)/(1 + m1 * m2))
                      newthet = atan(tanthet)
                      if (newthet > bestchange && i < 7) {
                        bestchange = newthet
                        MGvar$Opt.dim <<- i
                      }
                    }
                    for (i in 1:(cap - 2)) {
                      if (MGvar$Opt.dim == i) {
                        MGvar$screepoints.best[i] <<- MGvar$scree.stress[i]
                      }
                      else {
                        MGvar$screepoints.best[i] <<- ""
                      }
                      if (ndim == i) {
                        MGvar$screepoints.current[i] <<- MGvar$scree.stress[i]
                      }
                      else {
                        MGvar$screepoints.current[i] <<- ""
                      }
                    }
                    ActiveScree()
                    tkrreplot(imgscree)
                  }
                  MGvar$breakscree <<- FALSE
                  MGvar$fromscree <<- FALSE
                }
                else {
                  BlankScree()
                }
            }
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 4, state = "active")
            MGvar$endprocess <<- "no"
        }
        ActivateAll()
        ResetUndo()
    }
    my.SMACOF.ordinal <- function(Delta.mat, X.ini, W.mat = NULL, 
        iter = MGvar$MDS.iter.max, eps = MGvar$MDS.tol) {
        start.time <- proc.time()
        n <- nrow(Delta.mat)
        if (is.null(row.names(Delta.mat))) 
            row.names(Delta.mat) <- 1:n
        Z <- X.ini
        dimnames(Z) <- list(dimnames(Delta.mat)[[1]], NULL)
        if (is.null(W.mat)) 
            W.mat <- matrix(1, nrow = n, ncol = n)
        D.Z.mat <- Euclid.afstand(Z)
        matrix.combined <- matrix(0, nrow = n * (n - 1)/2, ncol = 4)
        dimnames(matrix.combined) <- list(NULL, c("ID", "deltaij", 
            "dij", "dhatij"))
        matrix.combined[, 1] <- 1:(n * (n - 1)/2)
        matrix.combined[, 2] <- proximitymat2vector(Delta.mat)$vec
        matrix.combined[, 3] <- proximitymat2vector(D.Z.mat)$vec
        matrix.combined.ordered <- matrix.combined[order(matrix.combined[, 
            2]), ]
        matrix.combined.ordered[, 4] <- UpDownBlocks.algorithm(x = matrix.combined.ordered[, 
            3])$d.hat.vec
        matrix.combined.back <- matrix.combined.ordered[order(matrix.combined.ordered[, 
            1]), ]
        D.hat.mat <- vector2proximitymat(matrix.combined.back[, 
            4], n = n)
        D.hat.mat <- D.hat.mat * sqrt(n * (n - 1))/sqrt(sum(W.mat * 
            D.hat.mat * D.hat.mat))
        eta.sq.d.hat <- 0.5 * (sum(D.hat.mat * D.hat.mat * W.mat))
        eta.sq.X <- 0.5 * (sum(D.Z.mat * D.Z.mat * W.mat))
        rho.d.hat.X <- 0.5 * sum(D.hat.mat * D.Z.mat * W.mat)
        sigma.r.d.hat.X.ini <- eta.sq.d.hat + eta.sq.X - 2 * 
            rho.d.hat.X
        V <- -W.mat
        diag(V) <- 0
        diag(V) <- -apply(V, 1, sum)
        Vplus <- solve(V + matrix(1, nrow = n, ncol = n)) - matrix(1, 
            nrow = n, ncol = n)/(n * n)
        W.d.hat.mat <- D.hat.mat * W.mat
        B.z.mat <- ifelse(zapsmall(W.d.hat.mat) == 0, 0, W.d.hat.mat/D.Z.mat)
        diag(B.z.mat) <- 0
        B.z.mat <- -B.z.mat + diag(apply(B.z.mat, 2, sum))
        X.update <- Vplus %*% B.z.mat %*% Z
        tel <- 0
        stress.begin <- sigma.r.d.hat.X.ini
        performance <- data.frame(IterNo = tel, Stress = stress.begin)
        itnum <- 0
        stressval <- 0
        tel <- 0
        stress.begin <- Inf
        if (!MGvar$fromscree) {
            MGvar$stressitervec <<- as.vector(0)
        }
        if (MGvar$GenSet.UpConf == "yes") {
            MGvar$MDSmat <<- Z
            rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
            ActiveCoordMat()
            tabplot()
        }
        repeat {
            if (MGvar$endprocess == "yes") {
                break
            }
            tel <- tel + 1
            if (tel > MGvar$MDS.iter.max) {
                termin <- "Maximum number iterations reached. Increase maximum iterations!"
                tkmessageBox(message = termin)
                MGvar$breakscree <<- TRUE
                break
            }
            Z <- X.update
            D.Z.mat <- Euclid.afstand(Z)
            matrix.combined[, 3] <- proximitymat2vector(D.Z.mat)$vec
            matrix.combined.ordered <- matrix.combined[order(matrix.combined[, 
                2]), ]
            matrix.combined.ordered[, 4] <- UpDownBlocks.algorithm(x = matrix.combined.ordered[, 
                3])$d.hat.vec
            matrix.combined.back <- matrix.combined.ordered[order(matrix.combined.ordered[, 
                1]), ]
            D.hat.mat <- vector2proximitymat(matrix.combined.back[, 
                4], n = n)
            D.hat.mat <- D.hat.mat * sqrt(n * (n - 1))/sqrt(sum(W.mat * 
                D.hat.mat * D.hat.mat))
            eta.sq.d.hat <- 0.5 * (sum(D.hat.mat * D.hat.mat * 
                W.mat))
            eta.sq.X <- 0.5 * (sum(D.Z.mat * D.Z.mat * W.mat))
            rho.d.hat.X <- 0.5 * sum(D.hat.mat * D.Z.mat * W.mat)
            sigma.r.d.hat.X.new <- eta.sq.d.hat + eta.sq.X - 
                2 * rho.d.hat.X
            if (stress.begin - sigma.r.d.hat.X.new < eps) 
                break
            stress.begin <- sigma.r.d.hat.X.new
            W.d.hat.mat <- D.hat.mat * W.mat
            B.z.mat <- ifelse(zapsmall(W.d.hat.mat) == 0, 0, 
                W.d.hat.mat/D.Z.mat)
            diag(B.z.mat) <- 0
            B.z.mat <- -B.z.mat + diag(apply(B.z.mat, 2, sum))
            X.update <- Vplus %*% B.z.mat %*% Z
            performance <- rbind(performance, c(tel, stress.begin))
            if (MGvar$GenSet.UpConf == "yes") {
                MGvar$MDSmat <<- X.update
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                ActiveCoordMat()
                if (MGvar$MDS.dimensions == 2) {
                  tabplot()
                }
                if (MGvar$MDS.dimensions == 3) {
                  tkrreplot(img3Dstat)
                }
            }
            if (MGvar$GenSet.UpShep == "yes" && MGvar$GenSet.CalcShep == 
                "yes") {
                MGvar$MDSmat <<- X.update
                if (tk2notetab.text(mySecondaryNB) == "Shepard") {
                  tkrreplot(imgshep)
                }
                if (MGvar$EnShep.switch == "on") {
                  Replot.imgEShep()
                }
            }
            if (MGvar$GenSet.UpStress == "yes" && MGvar$GenSet.CalcStress == 
                "yes") {
                iter.time <- proc.time()
                MGvar$proctime <<- (iter.time - start.time)[3]
                if (tk2notetab.text(mySecondaryNB) == "Stress Plot") {
                  tkrreplot(imgstress)
                }
                if (tk2notetab.text(mySecondaryNB) == "Log Stress") {
                  tkrreplot(imgstress2)
                }
                if (MGvar$EnStress2.switch == "on") {
                  tkrreplot(MGcomp$imgEStress2, hscale = MGvar$ESt2hscale, 
                    vscale = MGvar$ESt2vscale)
                }
                if (MGvar$EnStress.switch == "on") {
                  tkrreplot(MGcomp$imgEStress, hscale = MGvar$ESt1hscale, 
                    vscale = MGvar$ESt1vscale)
                }
            }
            if (MGvar$GenSet.UpProg == "yes") {
                tkconfigure(MDSGUI.ProgressBar, value = (tel/MGvar$MDS.iter.max) * 
                  100)
            }
            row.names(Z) <- row.names(Delta.mat)
            stress.begin <- sigma.r.d.hat.X.new
            itnum <- append(itnum, tel)
            MGvar$stressval <<- append(MGvar$stressval, sigma.r.d.hat.X.new)
            if (!MGvar$fromscree) {
                MGvar$stressitervec[tel] <<- StressCalculation(MGvar$MDSmat, 
                  MGvar$distmat, isMet = FALSE)
            }
        }
        if (MGvar$endprocess == "no") {
            list(conf = Z, iter = tel, stress = sigma.r.d.hat.X.new)
        }
    }
    UpDownBlocks.algorithm <- function(x, w = rep(1, length(x)), 
        block = weighted.mean) {
        isUpSatisfied <- function(x, i) (i == length(x)) || (x[i] <= 
            x[i + 1])
        isDownSatisfied <- function(x, i) (i == 1) || (x[i - 
            1] <= x[i])
        mergeBlockup <- function(blocklist, blockvalues, x, w, 
            i, block) {
            n <- length(blockvalues)
            nn <- 1:n
            ii <- which(i + 1 != nn)
            blocklist[i, ] <- c(blocklist[i, 1], blocklist[i + 
                1, 2])
            indi <- blocklist[i, 1]:blocklist[i + 1, 2]
            blockvalues[i] <- block(x[indi], w[indi])
            blocklist <- blocklist[ii, ]
            if (length(ii) == 1) 
                dim(blocklist) <- c(1, 2)
            blockvalues <- blockvalues[ii]
            list(v = blockvalues, l = blocklist)
        }
        putBack <- function(n, blocklist, blockvalues) {
            x <- rep(0, n)
            nb <- length(blockvalues)
            for (i in 1:nb) {
                x[blocklist[i, 1]:blocklist[i, 2]] <- blockvalues[i]
            }
            return(x)
        }
        nblock <- n <- length(x)
        blocklist <- matrix(1:n, nrow = n, ncol = 2)
        blockvalues <- x
        active <- 1
        repeat {
            if (!isUpSatisfied(blockvalues, active)) {
                blockmerge <- mergeBlockup(blocklist, blockvalues, 
                  x, w, active, block)
                blockvalues <- blockmerge$v
                blocklist <- blockmerge$l
                nblock <- nblock - 1
                while (!isDownSatisfied(blockvalues, active)) {
                  blockmerge <- mergeBlockup(blocklist, blockvalues, 
                    x, w, active - 1, block)
                  blockvalues <- blockmerge$v
                  blocklist <- blockmerge$l
                  nblock <- nblock - 1
                  active <- active - 1
                }
            }
            else if (active == nblock) 
                break()
            else active <- active + 1
        }
        d.hat <- putBack(n, blocklist, blockvalues)
        list(d.vec = x, d.hat.vec = d.hat)
    }
    Euclid.afstand <- function(dat) {
        X = as.matrix(dist(dat))
        return(X)
    }
    proximitymat2vector <- function(dat) {
        V = dat[lower.tri(dat)]
        list(vec = V)
    }
    vector2proximitymat <- function(vec, n) {
        mat <- matrix(0, n, n)
        mat[lower.tri(mat)] <- vec
        mat <- mat + t(mat)
        return(mat)
    }
    MyNonMetricSmacofSym <- function(data = MGvar$activedata, 
        fromdistmat = FALSE, initial = TRUE) {
        DeactivateAll()
        if (nrow(data) == 1 && ncol(data) == 1) {
            tkmessageBox(message = "Invalid Data!", icon = "error")
        }
        else {
            MGvar$is.Metric <<- FALSE
            ActiveIsMetric()
            if (initial) {
                ClearRemindex()
                MGvar$distmat <<- MGvar$originaldistmat
                ActivedistMat()
            }
            distnameassign()
            if (tclvalue(MGvar$MDSops.startconfig) == "ClasScal") {
                startconfig = cmdscale(MGvar$distmat, k = MGvar$MDS.dimensions)
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "RandConf") {
                startconfig = runif(MGvar$MDS.dimensions * nrow(data))
                startconfig <- stand.1.range(matrix(startconfig, 
                  ncol = MGvar$MDS.dimensions))
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "ExistConf") {
                if (MGvar$startconfig.UseEC.plot == "Plot1") {
                  startconfig = MGvar$MDSmat.T1
                }
                if (MGvar$startconfig.UseEC.plot == "Plot2") {
                  startconfig = MGvar$MDSmat.T2
                }
                if (MGvar$startconfig.UseEC.plot == "Plot3") {
                  startconfig = MGvar$MDSmat.T3
                }
                if (MGvar$startconfig.UseEC.plot == "Plot4") {
                  startconfig = MGvar$MDSmat.T4
                }
                if (MGvar$startconfig.UseEC.plot == "Plot5") {
                  startconfig = MGvar$MDSmat.T5
                }
            }
            if (MGvar$CSinitConf == "no") {
                startconfig = MGvar$MDSmat
                MGvar$CSinitConf <<- "yes"
            }
            startconfig <- stand.1.range(matrix(startconfig, 
                ncol = MGvar$MDS.dimensions))
            ActivePos()
            MGvar$tempdims <<- MGvar$MDS.dimensions
            ActiveWatchCursor()
            MGvar$MDStype <<- "NM.Smac.Sym"
            SymNMetSmac = my.SMACOF.ordinal(Delta.mat = MGvar$distmat, 
                X.ini = startconfig, iter = MGvar$MDS.iter.max, 
                eps = MGvar$MDS.tol)
            if (MGvar$endprocess == "no") {
                MGvar$MDSStress <<- StressCalculation(SymNMetSmac$conf, 
                  MGvar$distmat, isMet = FALSE)
                if (MGvar$MDS.dimensions == 1) {
                  UniCoords = matrix(0, nrow = nrow(MGvar$distmat), 
                    ncol = 2)
                  for (i in 1:nrow(UniCoords)) {
                    UniCoords[i, ] = c(SymNMetSmac$conf[i], 0)
                  }
                  MGvar$MDSmat <<- UniCoords
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$activeplot.title <<- paste("Unidimensional Non-Metric SMACOF on", 
                    MGvar$datatitle, "data")
                }
                if (MGvar$MDS.dimensions == 2) {
                  MGvar$MDSmat <<- SymNMetSmac$conf
                  MGvar$PlottingDimX <<- 1
                  MGvar$PlottingDimY <<- 2
                  MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                  MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$activeplot.title <<- paste("Non-Metric SMACOF on", 
                    MGvar$datatitle, "data")
                }
                if (MGvar$MDS.dimensions > 3) {
                  MGvar$MDSmat.LD <<- SymNMetSmac$conf
                  LargeDimensions()
                  MGvar$activeplot.title <<- paste("Non-Metric SMACOF on", 
                    MGvar$datatitle, "data: Dim", MGvar$PlottingDimX, 
                    "vs. Dim", MGvar$PlottingDimY)
                  MGvar$sthreeDplot.title <<- paste("3D Non-Metric SMACOF on", 
                    MGvar$datatitle, "data")
                  MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                  MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                  ActiveLabs()
                }
                if (MGvar$MDS.dimensions == 3) {
                  MGvar$MDSmat <<- SymNMetSmac$conf
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$sthreeDplot.title <<- paste("3D non-Metric SMACOF on", 
                    MGvar$datatitle, "data")
                  threedimplotting()
                }
                if (MGvar$tempdims != 3) {
                  ActiveMDSdims()
                  ActivedMeas()
                  ActiveTab()
                  ActiveStress()
                  ActiveTabDims()
                  ActiveCoordMat()
                  ActiveTitle()
                  ActiveIndex()
                  ActiveType()
                  ActiveZoomSwitch()
                  ActiveIterTab()
                  ActiveTolTab()
                  MGvar$activeplot.cex <<- 0.6
                  tableupdate()
                  if (MGvar$GenSet.ClearCols == "yes") {
                    PointColInitialise()
                  }
                  if (MGvar$GenSet.CalcMDS == "yes") {
                    tabplot()
                  }
                  else {
                    PlotBlank()
                  }
                }
                if (MGvar$GenSet.CalcShep == "yes") {
                  MGvar$shep.firstrun <<- "yes"
                  AShepdistmat <<- MGvar$distmat
                  AShepMDSmat <<- MGvar$MDSmat
                  ActiveShep()
                  tkrreplot(imgshep)
                }
                else {
                  BlankShep()
                }
                tkfocus(mytt)
                if (MGvar$GenSet.CalcScree == "yes") {
                  ndim = MGvar$MDS.dimensions
                  bestchange = 0
                  cap = min(MGvar$maxdims, 12)
                  MGvar$fromscree <<- TRUE
                  if (MGvar$GenSet.UpShep == "yes") {
                    MGvar$GenSet.UpShep <<- "no"
                    sheptracker = 1
                  }
                  else {
                    sheptracker = 0
                  }
                  if (MGvar$GenSet.UpConf == "yes") {
                    MGvar$GenSet.UpConf <<- "no"
                    conftracker = 1
                  }
                  else {
                    conftracker = 0
                  }
                  if (MGvar$GenSet.UpStress == "yes") {
                    MGvar$GenSet.UpStress <<- "no"
                    stresstracker = 1
                  }
                  else {
                    stresstracker = 0
                  }
                  for (i in 1:cap) {
                    if (MGvar$endprocess == "yes") {
                      break
                      MGvar$breakscree <<- TRUE
                    }
                    if (MGvar$breakscree) {
                      break
                    }
                    newSmac <- my.SMACOF.ordinal(Delta.mat = MGvar$distmat, 
                      X.ini = cmdscale(MGvar$distmat, k = i), 
                      iter = MGvar$MDS.iter.max)$conf
                    MGvar$scree.stress[i] <<- StressCalculation(newSmac, 
                      MGvar$distmat, isMet = TRUE)
                  }
                  if (sheptracker == 1) {
                    MGvar$GenSet.UpShep <<- "yes"
                    sheptracker = 0
                  }
                  if (conftracker == 1) {
                    MGvar$GenSet.UpConf <<- "yes"
                    conftracker = 0
                  }
                  if (stresstracker == 1) {
                    MGvar$GenSet.UpStress <<- "yes"
                    stresstracker = 0
                  }
                  if (!MGvar$breakscree) {
                    for (i in 2:(cap - 2)) {
                      m1 = -(MGvar$scree.stress[i - 1] - MGvar$scree.stress[i])
                      m2 = -(MGvar$scree.stress[i] - MGvar$scree.stress[i + 
                        1])
                      tanthet = abs((m1 - m2)/(1 + m1 * m2))
                      newthet = atan(tanthet)
                      if (newthet > bestchange && i < 7) {
                        bestchange = newthet
                        MGvar$Opt.dim <<- i
                      }
                    }
                    for (i in 1:(cap - 2)) {
                      if (MGvar$Opt.dim == i) {
                        MGvar$screepoints.best[i] <<- MGvar$scree.stress[i]
                      }
                      else {
                        MGvar$screepoints.best[i] <<- ""
                      }
                      if (ndim == i) {
                        MGvar$screepoints.current[i] <<- MGvar$scree.stress[i]
                      }
                      else {
                        MGvar$screepoints.current[i] <<- ""
                      }
                    }
                    ActiveScree()
                    tkrreplot(imgscree)
                  }
                  MGvar$breakscree <<- FALSE
                  MGvar$fromscree <<- FALSE
                }
                else {
                  BlankScree()
                }
            }
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 4, state = "active")
            MGvar$endprocess <<- "no"
        }
        ActivateAll()
        ResetUndo()
    }
    MyKruskalsMDS <- function(data = MGvar$activedata, fromdistmat = FALSE, 
        initial = TRUE) {
        DeactivateAll()
        if (nrow(data) == 1 && ncol(data) == 1) {
            tkmessageBox(message = "Invalid Data!", icon = "error")
        }
        else {
            MGvar$is.Metric <<- FALSE
            ActiveIsMetric()
            if (initial) {
                ClearRemindex()
                MGvar$distmat <<- MGvar$originaldistmat
                ActivedistMat()
            }
            distnameassign()
            if (tclvalue(MGvar$MDSops.startconfig) == "ClasScal") {
                startconfig = cmdscale(MGvar$distmat, k = MGvar$MDS.dimensions)
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "RandConf") {
                startconfig = runif(MGvar$MDS.dimensions * nrow(data))
                startconfig <- stand.1.range(matrix(startconfig, 
                  ncol = MGvar$MDS.dimensions))
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "ExistConf") {
                if (MGvar$startconfig.UseEC.plot == "Plot1") {
                  startconfig = MGvar$MDSmat.T1
                }
                if (MGvar$startconfig.UseEC.plot == "Plot2") {
                  startconfig = MGvar$MDSmat.T2
                }
                if (MGvar$startconfig.UseEC.plot == "Plot3") {
                  startconfig = MGvar$MDSmat.T3
                }
                if (MGvar$startconfig.UseEC.plot == "Plot4") {
                  startconfig = MGvar$MDSmat.T4
                }
                if (MGvar$startconfig.UseEC.plot == "Plot5") {
                  startconfig = MGvar$MDSmat.T5
                }
            }
            if (MGvar$CSinitConf == "no") {
                startconfig = MGvar$MDSmat
                MGvar$CSinitConf <<- "yes"
            }
            ActivePos()
            MGvar$tempdims <<- MGvar$MDS.dimensions
            ActiveWatchCursor()
            MGvar$MDStype <<- "Kruskal"
            myKrusk = isoMDS(MGvar$distmat, y = startconfig, 
                k = MGvar$MDS.dimensions, tol = MGvar$MDS.tol)
            MGvar$MDSStress <<- StressCalculation(myKrusk$points, 
                MGvar$distmat, isMet = FALSE)
            if (MGvar$MDS.dimensions == 1) {
                UniCoords = matrix(0, nrow = nrow(MGvar$distmat), 
                  ncol = 2)
                for (i in 1:nrow(UniCoords)) {
                  UniCoords[i, ] = c(myKrusk$points[i], 0)
                }
                MGvar$MDSmat <<- UniCoords
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$activeplot.title <<- paste("Unidimensional Kruskals Non-Metric MDS on", 
                  MGvar$datatitle, "data")
            }
            if (MGvar$MDS.dimensions == 2) {
                MGvar$MDSmat <<- myKrusk$points
                MGvar$PlottingDimX <<- 1
                MGvar$PlottingDimY <<- 2
                MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$activeplot.title <<- paste("Kruskals Non-Metric MDS on", 
                  MGvar$datatitle, "data")
            }
            if (MGvar$MDS.dimensions > 3) {
                MGvar$MDSmat.LD <<- myKrusk$points
                LargeDimensions()
                MGvar$activeplot.title <<- paste("Kruskals Non-Metric MDS on", 
                  MGvar$datatitle, "data: Dim", MGvar$PlottingDimX, 
                  "vs. Dim", MGvar$PlottingDimY)
                MGvar$sthreeDplot.title <<- paste("3D Kruskals Non-Metric MDS on", 
                  MGvar$datatitle, "data")
                MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                ActiveLabs()
            }
            if (MGvar$MDS.dimensions == 3) {
                MGvar$MDSmat <<- myKrusk$points
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$sthreeDplot.title <<- paste("3D Kruskals Non-Metric MDS on", 
                  MGvar$datatitle, "data")
                threedimplotting()
            }
            if (MGvar$tempdims != 3) {
                ActiveMDSdims()
                ActivedMeas()
                ActiveTab()
                ActiveStress()
                ActiveTabDims()
                ActiveCoordMat()
                ActiveTitle()
                ActiveIndex()
                ActiveType()
                ActiveZoomSwitch()
                ActiveTolTab()
                MGvar$activeplot.cex <<- 0.6
                tableupdate()
                if (MGvar$GenSet.ClearCols == "yes") {
                  PointColInitialise()
                }
                if (MGvar$GenSet.CalcMDS == "yes") {
                  tabplot()
                }
                else {
                  PlotBlank()
                }
            }
            if (MGvar$GenSet.CalcShep == "yes") {
                MGvar$shep.firstrun <<- "yes"
                AShepdistmat <<- MGvar$distmat
                AShepMDSmat <<- MGvar$MDSmat
                ActiveShep()
                tkrreplot(imgshep)
            }
            else {
                BlankShep()
            }
            tkfocus(mytt)
            if (MGvar$GenSet.CalcScree == "yes") {
                ndim = MGvar$MDS.dimensions
                bestchange = 0
                cap = min(MGvar$maxdims, 12)
                for (i in 1:cap) {
                  nextKrusk <- isoMDS(MGvar$distmat, k = i)$points
                  MGvar$scree.stress[i] <<- StressCalculation(nextKrusk, 
                    MGvar$distmat, isMet = FALSE)
                }
                for (i in 2:(cap - 2)) {
                  m1 = -(MGvar$scree.stress[i - 1] - MGvar$scree.stress[i])
                  m2 = -(MGvar$scree.stress[i] - MGvar$scree.stress[i + 
                    1])
                  tanthet = abs((m1 - m2)/(1 + m1 * m2))
                  newthet = atan(tanthet)
                  if (newthet > bestchange && i < 7) {
                    bestchange = newthet
                    MGvar$Opt.dim <<- i
                  }
                }
                for (i in 1:(cap - 2)) {
                  if (MGvar$Opt.dim == i) {
                    MGvar$screepoints.best[i] <<- MGvar$scree.stress[i]
                  }
                  else {
                    MGvar$screepoints.best[i] <<- ""
                  }
                  if (ndim == i) {
                    MGvar$screepoints.current[i] <<- MGvar$scree.stress[i]
                  }
                  else {
                    MGvar$screepoints.current[i] <<- ""
                  }
                }
                ActiveScree()
                tkrreplot(imgscree)
            }
            else {
                BlankScree()
            }
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 4, state = "active")
        }
        ActivateAll()
        ResetUndo()
    }
    SammonMapping <- function(data = MGvar$activedata, fromdistmat = FALSE, 
        initial = TRUE) {
        DeactivateAll()
        if (nrow(data) == 1 && ncol(data) == 1) {
            tkmessageBox(message = "Invalid Data!", icon = "error")
        }
        else {
            MGvar$is.Metric <<- FALSE
            ActiveIsMetric()
            if (initial) {
                ClearRemindex()
                MGvar$distmat <<- MGvar$originaldistmat
                ActivedistMat()
            }
            distnameassign()
            if (tclvalue(MGvar$MDSops.startconfig) == "ClasScal") {
                startconfig = cmdscale(MGvar$distmat, k = MGvar$MDS.dimensions)
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "RandConf") {
                startconfig = runif(MGvar$MDS.dimensions * nrow(data))
                startconfig <- stand.1.range(matrix(startconfig, 
                  ncol = MGvar$MDS.dimensions))
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "ExistConf") {
                if (MGvar$startconfig.UseEC.plot == "Plot1") {
                  startconfig = MGvar$MDSmat.T1
                }
                if (MGvar$startconfig.UseEC.plot == "Plot2") {
                  startconfig = MGvar$MDSmat.T2
                }
                if (MGvar$startconfig.UseEC.plot == "Plot3") {
                  startconfig = MGvar$MDSmat.T3
                }
                if (MGvar$startconfig.UseEC.plot == "Plot4") {
                  startconfig = MGvar$MDSmat.T4
                }
                if (MGvar$startconfig.UseEC.plot == "Plot5") {
                  startconfig = MGvar$MDSmat.T5
                }
            }
            if (MGvar$CSinitConf == "no") {
                startconfig = MGvar$MDSmat
                MGvar$CSinitConf <<- "yes"
            }
            ActivePos()
            MGvar$tempdims <<- MGvar$MDS.dimensions
            ActiveWatchCursor()
            MGvar$MDStype <<- "Sammon"
            mySam = sammon(MGvar$distmat, k = MGvar$MDS.dimensions, 
                y = startconfig, tol = MGvar$MDS.tol)
            MGvar$MDSStress <<- StressCalculation(mySam$points, 
                MGvar$distmat, isMet = FALSE)
            if (MGvar$MDS.dimensions == 1) {
                UniCoords = matrix(0, nrow = nrow(MGvar$distmat), 
                  ncol = 2)
                for (i in 1:nrow(UniCoords)) {
                  UniCoords[i, ] = c(mySam$points[i], 0)
                }
                MGvar$MDSmat <<- UniCoords
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$activeplot.title <<- paste("Unidimensional Sammon Mapping on", 
                  MGvar$datatitle, "data")
            }
            if (MGvar$MDS.dimensions == 2) {
                MGvar$MDSmat <<- mySam$points
                MGvar$PlottingDimX <<- 1
                MGvar$PlottingDimY <<- 2
                MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$activeplot.title <<- paste("Sammon Mapping on", 
                  MGvar$datatitle, "data")
            }
            if (MGvar$MDS.dimensions > 3) {
                MGvar$MDSmat.LD <<- mySam$points
                LargeDimensions()
                MGvar$activeplot.title <<- paste("Sammon Mapping on", 
                  MGvar$datatitle, "data: Dim", MGvar$PlottingDimX, 
                  "vs. Dim", MGvar$PlottingDimY)
                MGvar$sthreeDplot.title <<- paste("3D Sammon Mapping on", 
                  MGvar$datatitle, "data")
                MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                ActiveLabs()
            }
            if (MGvar$MDS.dimensions == 3) {
                MGvar$MDSmat <<- mySam$points
                rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                MGvar$sthreeDplot.title <<- paste("3D Sammon Mapping on", 
                  MGvar$datatitle, "data")
                threedimplotting()
            }
            if (MGvar$tempdims != 3) {
                ActiveMDSdims()
                ActivedMeas()
                ActiveTab()
                ActiveStress()
                ActiveTabDims()
                ActiveCoordMat()
                ActiveTitle()
                ActiveIndex()
                ActiveType()
                ActiveZoomSwitch()
                ActiveTolTab()
                MGvar$activeplot.cex <<- 0.6
                tableupdate()
                if (MGvar$GenSet.ClearCols == "yes") {
                  PointColInitialise()
                }
                if (MGvar$GenSet.CalcMDS == "yes") {
                  tabplot()
                }
                else {
                  PlotBlank()
                }
            }
            if (MGvar$GenSet.CalcShep == "yes") {
                MGvar$shep.firstrun <<- "yes"
                AShepdistmat <<- MGvar$distmat
                AShepMDSmat <<- MGvar$MDSmat
                ActiveShep()
                tkrreplot(imgshep)
            }
            else {
                BlankShep()
            }
            tkfocus(mytt)
            if (MGvar$GenSet.CalcScree == "yes") {
                ndim = MGvar$MDS.dimensions
                bestchange = 0
                cap = min(MGvar$maxdims, 12)
                for (i in 1:(cap)) {
                  nextSam <- sammon(MGvar$distmat, k = i, niter = MGvar$MDS.iter.max)$points
                  MGvar$scree.stress[i] <<- StressCalculation(nextSam, 
                    MGvar$distmat, isMet = FALSE)
                }
                for (i in 2:(cap - 2)) {
                  m1 = -(MGvar$scree.stress[i - 1] - MGvar$scree.stress[i])
                  m2 = -(MGvar$scree.stress[i] - MGvar$scree.stress[i + 
                    1])
                  tanthet = abs((m1 - m2)/(1 + m1 * m2))
                  newthet = atan(tanthet)
                  if (newthet > bestchange && i < 7) {
                    bestchange = newthet
                    MGvar$Opt.dim <<- i
                  }
                }
                for (i in 1:(cap - 2)) {
                  if (MGvar$Opt.dim == i) {
                    MGvar$screepoints.best[i] <<- MGvar$scree.stress[i]
                  }
                  else {
                    MGvar$screepoints.best[i] <<- ""
                  }
                  if (ndim == i) {
                    MGvar$screepoints.current[i] <<- MGvar$scree.stress[i]
                  }
                  else {
                    MGvar$screepoints.current[i] <<- ""
                  }
                }
                ActiveScree()
                tkrreplot(imgscree)
            }
            else {
                BlankScree()
            }
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 4, state = "active")
        }
        ActivateAll()
        ResetUndo()
    }
    MetLSC <- function(Delta, start.config = cmdscale(Delta), 
        dims = MGvar$MDS.dimensions) {
        stress.func <- function(y.vec, Delta.mat, f.func = function(x) {
            x
        }, dim = dims) {
            Y <- matrix(y.vec, ncol = dim)
            D.mat <- as.matrix(dist(Y))
            diff.sq <- (f.func(Delta.mat) - D.mat)^2
            diff.sq <- diff.sq[lower.tri(diff.sq)]
            sum(diff.sq)/sum(D.mat[lower.tri(D.mat)]^2)
        }
        Delta <- as.matrix(Delta)
        y <- as.vector(start.config)
        Y <- matrix(optim(y, stress.func, Delta.mat = Delta)$par, 
            ncol = dims)
        list(conf = Y)
    }
    MetricLeastSquaresScaling <- function(data = MGvar$activedata, 
        fromdistmat = FALSE, initial = TRUE) {
        DeactivateAll()
        if (nrow(data) == 1 && ncol(data) == 1) {
            tkmessageBox(message = "Invalid Data!", icon = "error")
        }
        else {
            MGvar$is.Metric <<- TRUE
            ActiveIsMetric()
            if (initial) {
                ClearRemindex()
                MGvar$distmat <<- MGvar$originaldistmat
                ActivedistMat()
            }
            distnameassign()
            if (tclvalue(MGvar$MDSops.startconfig) == "ClasScal") {
                startconfig = cmdscale(MGvar$distmat, k = MGvar$MDS.dimensions)
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "RandConf") {
                startconfig = runif(MGvar$MDS.dimensions * nrow(data))
                startconfig <- stand.1.range(matrix(startconfig, 
                  ncol = MGvar$MDS.dimensions))
            }
            if (tclvalue(MGvar$MDSops.startconfig) == "ExistConf") {
                if (MGvar$startconfig.UseEC.plot == "Plot1") {
                  startconfig = MGvar$MDSmat.T1
                }
                if (MGvar$startconfig.UseEC.plot == "Plot2") {
                  startconfig = MGvar$MDSmat.T2
                }
                if (MGvar$startconfig.UseEC.plot == "Plot3") {
                  startconfig = MGvar$MDSmat.T3
                }
                if (MGvar$startconfig.UseEC.plot == "Plot4") {
                  startconfig = MGvar$MDSmat.T4
                }
                if (MGvar$startconfig.UseEC.plot == "Plot5") {
                  startconfig = MGvar$MDSmat.T5
                }
            }
            if (MGvar$CSinitConf == "no") {
                startconfig = MGvar$MDSmat
                MGvar$CSinitConf <<- "yes"
            }
            ActivePos()
            MGvar$tempdims <<- MGvar$MDS.dimensions
            ActiveWatchCursor()
            MGvar$MDStype <<- "Met.LeastSqs."
            MetLSCres <- MetLSC(Delta = MGvar$distmat, start.config = startconfig)
            if (MGvar$endprocess == "no") {
                MGvar$MDSStress <<- StressCalculation(MetLSCres$conf, 
                  MGvar$distmat, isMet = TRUE)
                if (MGvar$MDS.dimensions == 1) {
                  UniCoords = matrix(0, nrow = nrow(MGvar$distmat), 
                    ncol = 2)
                  for (i in 1:nrow(UniCoords)) {
                    UniCoords[i, ] = c(MetLSCres$conf[i], 0)
                  }
                  MGvar$MDSmat <<- UniCoords
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$activeplot.title <<- paste("Unidimensional Metric Least Squares Scaling on", 
                    MGvar$datatitle, "data")
                }
                if (MGvar$MDS.dimensions == 2) {
                  MGvar$MDSmat <<- MetLSCres$conf
                  MGvar$PlottingDimX <<- 1
                  MGvar$PlottingDimY <<- 2
                  MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                  MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$activeplot.title <<- paste("Metric Least Squares Scaling on", 
                    MGvar$datatitle, "data")
                }
                if (MGvar$MDS.dimensions > 3) {
                  MGvar$MDSmat.LD <<- MetLSCres$conf
                  LargeDimensions()
                  MGvar$activeplot.title <<- paste("Metric Least Squares Scaling on", 
                    MGvar$datatitle, "data: Dim", MGvar$PlottingDimX, 
                    "vs. Dim", MGvar$PlottingDimY)
                  MGvar$sthreeDplot.title <<- paste("3D Metric SMACOF on", 
                    MGvar$datatitle, "data")
                  MGvar$activeplot.xlab <<- paste("Dim", MGvar$PlottingDimX)
                  MGvar$activeplot.ylab <<- paste("Dim", MGvar$PlottingDimY)
                  ActiveLabs()
                }
                if (MGvar$MDS.dimensions == 3) {
                  MGvar$MDSmat <<- MetLSCres$conf
                  rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
                  MGvar$sthreeDplot.title <<- paste("3D Metric Least Squares Scaling on", 
                    MGvar$datatitle, "data")
                  threedimplotting()
                }
                if (MGvar$tempdims != 3) {
                  ActiveMDSdims()
                  ActivedMeas()
                  ActiveTab()
                  ActiveStress()
                  ActiveTabDims()
                  ActiveCoordMat()
                  ActiveTitle()
                  ActiveIndex()
                  ActiveType()
                  ActiveZoomSwitch()
                  ActiveIterTab()
                  ActiveTolTab()
                  MGvar$activeplot.cex <<- 0.6
                  tableupdate()
                  if (MGvar$GenSet.ClearCols == "yes") {
                    PointColInitialise()
                  }
                  if (MGvar$GenSet.CalcMDS == "yes") {
                    tabplot()
                  }
                  else {
                    PlotBlank()
                  }
                }
                if (MGvar$GenSet.CalcShep == "yes") {
                  MGvar$shep.firstrun <<- "yes"
                  AShepdistmat <<- MGvar$distmat
                  AShepMDSmat <<- MGvar$MDSmat
                  ActiveShep()
                  tkrreplot(imgshep)
                }
                else {
                  BlankShep()
                }
                tkfocus(mytt)
                if (MGvar$GenSet.CalcScree == "yes") {
                  ndim = MGvar$MDS.dimensions
                  bestchange = 0
                  cap = min(MGvar$maxdims, 12)
                  MGvar$fromscree <<- TRUE
                  if (MGvar$GenSet.UpShep == "yes") {
                    MGvar$GenSet.UpShep <<- "no"
                    sheptracker = 1
                  }
                  else {
                    sheptracker = 0
                  }
                  if (MGvar$GenSet.UpConf == "yes") {
                    MGvar$GenSet.UpConf <<- "no"
                    conftracker = 1
                  }
                  else {
                    conftracker = 0
                  }
                  if (MGvar$GenSet.UpStress == "yes") {
                    MGvar$GenSet.UpStress <<- "no"
                    stresstracker = 1
                  }
                  else {
                    stresstracker = 0
                  }
                  for (i in 1:cap) {
                    if (MGvar$endprocess == "yes") {
                      break
                      MGvar$breakscree <<- TRUE
                    }
                    if (MGvar$breakscree) {
                      break
                    }
                    newMetLSC <- MetLSC(Delta = MGvar$distmat, 
                      cmdscale(MGvar$distmat, k = i), dims = i)$conf
                    MGvar$scree.stress[i] <<- StressCalculation(newMetLSC, 
                      MGvar$distmat, isMet = TRUE)
                  }
                  if (sheptracker == 1) {
                    MGvar$GenSet.UpShep <<- "yes"
                    sheptracker = 0
                  }
                  if (conftracker == 1) {
                    MGvar$GenSet.UpConf <<- "yes"
                    conftracker = 0
                  }
                  if (stresstracker == 1) {
                    MGvar$GenSet.UpStress <<- "yes"
                    stresstracker = 0
                  }
                  if (!MGvar$breakscree) {
                    for (i in 2:(cap - 2)) {
                      m1 = -(MGvar$scree.stress[i - 1] - MGvar$scree.stress[i])
                      m2 = -(MGvar$scree.stress[i] - MGvar$scree.stress[i + 
                        1])
                      tanthet = abs((m1 - m2)/(1 + m1 * m2))
                      newthet = atan(tanthet)
                      if (newthet > bestchange && i < 7) {
                        bestchange = newthet
                        MGvar$Opt.dim <<- i
                      }
                    }
                    for (i in 1:(cap - 2)) {
                      if (MGvar$Opt.dim == i) {
                        MGvar$screepoints.best[i] <<- MGvar$scree.stress[i]
                      }
                      else {
                        MGvar$screepoints.best[i] <<- ""
                      }
                      if (ndim == i) {
                        MGvar$screepoints.current[i] <<- MGvar$scree.stress[i]
                      }
                      else {
                        MGvar$screepoints.current[i] <<- ""
                      }
                    }
                    ActiveScree()
                    tkrreplot(imgscree)
                  }
                  MGvar$breakscree <<- FALSE
                  MGvar$fromscree <<- FALSE
                }
                else {
                  BlankScree()
                }
            }
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            tkentryconfigure(EditDataMenu, 4, state = "active")
            MGvar$endprocess <<- "no"
        }
        ActivateAll()
        ResetUndo()
    }
    threedimplotting <- function() {
        threeD <- tktoplevel()
        tkwm.resizable(threeD, "0", "0")
        tkwm.deiconify(threeD)
        tkwm.title(threeD, "3D Plotting")
        tkwm.geometry(threeD, "420x250")
        TDcanvas <- tkcanvas(threeD, width = "420", height = "250", 
            bg = col.sec)
        tkplace(TDcanvas, `in` = threeD)
        frame3D <- tkwidget(threeD, "TitleFrame", text = "3D Plot Location", 
            background = "white")
        tkplace(frame3D, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            relheight = 0.8, `in` = threeD)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frame3D, text = "A three dimensional plot will now be created. This plot may been\ndone as a static image in the Static 3D Plot tab or as an interactive\nRGL plot in an external window (Or both!). Please make your\nchoice below.", 
            font = fontsmall), relx = 0.08, rely = 0.1, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "Static plot in Tab", 
            background = "white"), relx = 0.15, rely = 0.45, 
            `in` = frame3D)
        Cstat.CB <- tk2checkbutton(threeD)
        Cstat.val <- tclVar("1")
        tkconfigure(Cstat.CB, variable = Cstat.val)
        tkplace(Cstat.CB, relx = 0.75, rely = 0.45, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "External RGL plot", 
            background = "white"), relx = 0.15, rely = 0.6, `in` = frame3D)
        Crgl.CB <- tk2checkbutton(threeD)
        Crgl.val <- tclVar("1")
        tkconfigure(Crgl.CB, variable = Crgl.val)
        tkplace(Crgl.CB, relx = 0.75, rely = 0.6, `in` = frame3D)
        On.3dplot <- function() {
            StatChoice = as.character(tclvalue(Cstat.val))
            if (StatChoice == "1") {
                tk2notetab.select(myPlottingNB, "Static 3D Plot")
                tkrreplot(img3Dstat)
                MGvar$MDSStress.3S <<- round(MGvar$MDSStress, 
                  3)
                MGvar$sthreeDplot.xlab <<- paste("Dim", MGvar$PlottingDimX.3D)
                MGvar$sthreeDplot.ylab <<- paste("Dim", MGvar$PlottingDimY.3D)
                MGvar$sthreeDplot.zlab <<- paste("Dim", MGvar$PlottingDimZ.3D)
                MGvar$TabDims.3S <<- paste(MGvar$PlottingDimX.3D, 
                  ",", MGvar$PlottingDimY.3D, "&", MGvar$PlottingDimZ.3D, 
                  sep = "")
                MGvar$MDS.dimensions.3S <<- MGvar$MDS.dimensions
                MGvar$dMeas.3S <<- MGvar$dMeas
                MGvar$MDStype.3S <<- MGvar$MDStype
                tableupdate()
            }
            RGLChoice = as.character(tclvalue(Crgl.val))
            if (RGLChoice == "1") {
                plotting3D(MGvar$MDSmat)
                MGvar$MDSStress.3R <<- round(MGvar$MDSStress, 
                  3)
                MGvar$TabDims.3R <<- paste(MGvar$PlottingDimX.3D, 
                  ",", MGvar$PlottingDimY.3D, "&", MGvar$PlottingDimZ.3D, 
                  sep = "")
                MGvar$MDS.dimensions.3R <<- MGvar$MDS.dimensions
                MGvar$dMeas.3R <<- MGvar$dMeas
                MGvar$MDStype.3R <<- MGvar$MDStype
                tableupdate()
            }
            tkdestroy(threeD)
        }
        On.cancel <- function() {
            tkdestroy(threeD)
        }
        tkplace(tkbutton(threeD, text = "Plot", width = 15, command = function() On.3dplot()), 
            relx = 0.35, rely = 0.8, `in` = frame3D)
        tkplace(tkbutton(threeD, text = "Cancel", width = 15, 
            command = function() On.cancel()), relx = 0.365, 
            rely = 0.85, `in` = threeD)
        tkfocus(threeD)
        tkwait.window(threeD)
    }
    ActiveTabChanges <- function() {
        if (MGvar$ActivePlottingTab == "Tab1") {
            MGvar$indexLabeled.T1 <<- MGvar$indexLabeled
        }
        if (MGvar$ActivePlottingTab == "Tab2") {
            MGvar$indexLabeled.T2 <<- MGvar$indexLabeled
        }
        if (MGvar$ActivePlottingTab == "Tab3") {
            MGvar$indexLabeled.T3 <<- MGvar$indexLabeled
        }
        if (MGvar$ActivePlottingTab == "Tab4") {
            MGvar$indexLabeled.T4 <<- MGvar$indexLabeled
        }
        if (MGvar$ActivePlottingTab == "Tab5") {
            MGvar$indexLabeled.T5 <<- MGvar$indexLabeled
        }
        DoNothing = function() {
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tk2notetab.select(myPlottingNB, "Plot1")
            MGvar$MDSmat <<- MGvar$MDSmat.T1
            MGvar$distmat <<- MGvar$distmat.T1
            MGvar$scree.stress <<- MGvar$scree.stress.T1
            MGvar$screepoints.current <<- MGvar$screepoints.current.T1
            MGvar$screepoints.best <<- MGvar$screepoints.best.T1
            MGvar$activeplot.title <<- MGvar$activeplot.title.T1
            MGvar$indexLabeled <<- MGvar$indexLabeled.T1
            MGvar$MDSStress <<- MGvar$MDSStress.T1
            MGvar$Opt.dim <<- MGvar$Opt.dim.T1
            MGvar$dMeas <<- MGvar$dMeas.T1
            AShepdistmat <<- MGvar$distmat.T1
            AShepMDSmat <<- MGvar$MDSmat.T1
            tclvalue(PlotText) <- paste("Active Plot is Plot1")
            tkbind(img, "<Button-1>", OnPlotLeftClick)
            tkbind(img2, "<Button-1>", DoNothing)
            tkbind(img3, "<Button-1>", DoNothing)
            tkbind(img4, "<Button-1>", DoNothing)
            tkbind(img5, "<Button-1>", DoNothing)
            tkbind(img, "<B1-Motion>", BrushingPointMove)
            tkbind(img2, "<B1-Motion>", DoNothing)
            tkbind(img3, "<B1-Motion>", DoNothing)
            tkbind(img4, "<B1-Motion>", DoNothing)
            tkbind(img5, "<B1-Motion>", DoNothing)
            tkbind(img, "<ButtonRelease-1>", OnRelease.Main)
            tkbind(img2, "<ButtonRelease-1>", DoNothing)
            tkbind(img3, "<ButtonRelease-1>", DoNothing)
            tkbind(img4, "<ButtonRelease-1>", DoNothing)
            tkbind(img5, "<ButtonRelease-1>", DoNothing)
            tkconfigure(img, cursor = "hand2")
            tkconfigure(img2, cursor = "arrow")
            tkconfigure(img3, cursor = "arrow")
            tkconfigure(img4, cursor = "arrow")
            tkconfigure(img5, cursor = "arrow")
            MGvar$MDStype <<- MGvar$MDStype.T1
            MGvar$newCoords <<- MGvar$newCoords.T1
            MGvar$zoomedplot.title.show <<- MGvar$activeplot.title.show.T1
            MGvar$zoomedplot.distmeas <<- MGvar$activeplot.distmeas.T1
            MGvar$zoomedplot.xlab <<- MGvar$activeplot.xlab.T1
            MGvar$zoomedplot.ylab <<- MGvar$activeplot.ylab.T1
            MGvar$zoomedplot.bg <<- MGvar$activeplot.bg.T1
            MGvar$zoomedplot.cex <<- MGvar$activeplot.cex.T1
            MGvar$zoomedplot.labs <<- MGvar$activeplot.labs.T1
            MGvar$zoomedplot.showpoints <<- MGvar$activeplot.showpoints.T1
            MGvar$zoomedplot.pointcol <<- MGvar$activeplot.pointcol.T1
            MGvar$zoomedplot.type <<- MGvar$activeplot.type.T1
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.yaxt.T1
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.xaxt.T1
            MGvar$zoomedplot.axescol <<- MGvar$activeplot.axescol.T1
            MGvar$zoomedplot.showreg <<- MGvar$activeplot.showreg.T1
            MGvar$zoomedplot.regcol <<- MGvar$activeplot.regcol.T1
            MGvar$removedpointsactivedata <<- MGvar$removedpointsactivedata.T1
            MGvar$removedpoints <<- MGvar$removedpoints.T1
            MGvar$remindex <<- MGvar$remindex.T1
            MGvar$RemovedPoints <<- MGvar$RemovedPoints.T1
            MGvar$remAxcompindex <<- MGvar$remAxcompindex.T1
            MGvar$RemovedAxes <<- MGvar$RemovedAxes.T1
            MGvar$remPcompindex <<- MGvar$remPcompindex.T1
            MGvar$DistFunc <<- MGvar$DistFunc.T1
            MGvar$is.Metric <<- MGvar$is.Metric.T1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkconfigure(img2, state = "active")
            tk2notetab.select(myPlottingNB, "Plot2")
            MGvar$MDSmat <<- MGvar$MDSmat.T2
            MGvar$distmat <<- MGvar$distmat.T2
            MGvar$scree.stress <<- MGvar$scree.stress.T2
            MGvar$screepoints.current <<- MGvar$screepoints.current.T2
            MGvar$screepoints.best <<- MGvar$screepoints.best.T2
            MGvar$activeplot.title <<- MGvar$activeplot.title.T2
            MGvar$indexLabeled <<- MGvar$indexLabeled.T2
            MGvar$MDSStress <<- MGvar$MDSStress.T2
            MGvar$Opt.dim <<- MGvar$Opt.dim.T2
            MGvar$dMeas <<- MGvar$dMeas.T2
            AShepdistmat <<- MGvar$distmat.T2
            AShepMDSmat <<- MGvar$MDSmat.T2
            tclvalue(PlotText) <- paste("Active Plot is Plot2")
            tkbind(img, "<Button-1>", DoNothing)
            tkbind(img2, "<Button-1>", OnPlotLeftClick)
            tkbind(img3, "<Button-1>", DoNothing)
            tkbind(img4, "<Button-1>", DoNothing)
            tkbind(img5, "<Button-1>", DoNothing)
            tkbind(img, "<B1-Motion>", DoNothing)
            tkbind(img2, "<B1-Motion>", BrushingPointMove)
            tkbind(img3, "<B1-Motion>", DoNothing)
            tkbind(img4, "<B1-Motion>", DoNothing)
            tkbind(img5, "<B1-Motion>", DoNothing)
            tkbind(img, "<ButtonRelease-1>", DoNothing)
            tkbind(img2, "<ButtonRelease-1>", OnRelease.Main)
            tkbind(img3, "<ButtonRelease-1>", DoNothing)
            tkbind(img4, "<ButtonRelease-1>", DoNothing)
            tkbind(img5, "<ButtonRelease-1>", DoNothing)
            tkconfigure(img, cursor = "arrow")
            tkconfigure(img2, cursor = "hand2")
            tkconfigure(img3, cursor = "arrow")
            tkconfigure(img4, cursor = "arrow")
            tkconfigure(img5, cursor = "arrow")
            MGvar$MDStype <<- MGvar$MDStype.T2
            MGvar$newCoords <<- MGvar$newCoords.T2
            MGvar$zoomedplot.title.show <<- MGvar$activeplot.title.show.T2
            MGvar$zoomedplot.distmeas <<- MGvar$activeplot.distmeas.T2
            MGvar$zoomedplot.xlab <<- MGvar$activeplot.xlab.T2
            MGvar$zoomedplot.ylab <<- MGvar$activeplot.ylab.T2
            MGvar$zoomedplot.bg <<- MGvar$activeplot.bg.T2
            MGvar$zoomedplot.cex <<- MGvar$activeplot.cex.T2
            MGvar$zoomedplot.labs <<- MGvar$activeplot.labs.T2
            MGvar$zoomedplot.showpoints <<- MGvar$activeplot.showpoints.T2
            MGvar$zoomedplot.pointcol <<- MGvar$activeplot.pointcol.T2
            MGvar$zoomedplot.type <<- MGvar$activeplot.type.T2
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.yaxt.T2
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.xaxt.T2
            MGvar$zoomedplot.axescol <<- MGvar$activeplot.axescol.T2
            MGvar$zoomedplot.showreg <<- MGvar$activeplot.showreg.T2
            MGvar$zoomedplot.regcol <<- MGvar$activeplot.regcol.T2
            MGvar$removedpointsactivedata <<- MGvar$removedpointsactivedata.T2
            MGvar$removedpoints <<- MGvar$removedpoints.T2
            MGvar$remindex <<- MGvar$remindex.T2
            MGvar$RemovedPoints <<- MGvar$RemovedPoints.T2
            MGvar$remAxcompindex <<- MGvar$remAxcompindex.T2
            MGvar$RemovedAxes <<- MGvar$RemovedAxes.T2
            MGvar$remPcompindex <<- MGvar$remPcompindex.T2
            MGvar$DistFunc <<- MGvar$DistFunc.T2
            MGvar$is.Metric <<- MGvar$is.Metric.T2
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkconfigure(img3, state = "active")
            tk2notetab.select(myPlottingNB, "Plot3")
            MGvar$MDSmat <<- MGvar$MDSmat.T3
            MGvar$distmat <<- MGvar$distmat.T3
            MGvar$scree.stress <<- MGvar$scree.stress.T3
            MGvar$screepoints.current <<- MGvar$screepoints.current.T3
            MGvar$screepoints.best <<- MGvar$screepoints.best.T3
            MGvar$activeplot.title <<- MGvar$activeplot.title.T3
            MGvar$indexLabeled <<- MGvar$indexLabeled.T3
            MGvar$Opt.dim <<- MGvar$Opt.dim.T3
            MGvar$MDSStress <<- MGvar$MDSStress.T3
            MGvar$dMeas <<- MGvar$dMeas.T3
            AShepdistmat <<- MGvar$distmat.T3
            AShepMDSmat <<- MGvar$MDSmat.T3
            tclvalue(PlotText) <- paste("Active Plot is Plot3")
            tkbind(img, "<Button-1>", DoNothing)
            tkbind(img2, "<Button-1>", DoNothing)
            tkbind(img3, "<Button-1>", OnPlotLeftClick)
            tkbind(img4, "<Button-1>", DoNothing)
            tkbind(img5, "<Button-1>", DoNothing)
            tkbind(img, "<B1-Motion>", DoNothing)
            tkbind(img2, "<B1-Motion>", DoNothing)
            tkbind(img3, "<B1-Motion>", BrushingPointMove)
            tkbind(img4, "<B1-Motion>", DoNothing)
            tkbind(img5, "<B1-Motion>", DoNothing)
            tkbind(img, "<ButtonRelease-1>", DoNothing)
            tkbind(img2, "<ButtonRelease-1>", DoNothing)
            tkbind(img3, "<ButtonRelease-1>", OnRelease.Main)
            tkbind(img4, "<ButtonRelease-1>", DoNothing)
            tkbind(img5, "<ButtonRelease-1>", DoNothing)
            tkconfigure(img, cursor = "arrow")
            tkconfigure(img2, cursor = "arrow")
            tkconfigure(img3, cursor = "hand2")
            tkconfigure(img4, cursor = "arrow")
            tkconfigure(img5, cursor = "arrow")
            MGvar$MDStype <<- MGvar$MDStype.T3
            MGvar$newCoords <<- MGvar$newCoords.T3
            MGvar$zoomedplot.title.show <<- MGvar$activeplot.title.show.T3
            MGvar$zoomedplot.distmeas <<- MGvar$activeplot.distmeas.T3
            MGvar$zoomedplot.xlab <<- MGvar$activeplot.xlab.T3
            MGvar$zoomedplot.ylab <<- MGvar$activeplot.ylab.T3
            MGvar$zoomedplot.bg <<- MGvar$activeplot.bg.T3
            MGvar$zoomedplot.cex <<- MGvar$activeplot.cex.T3
            MGvar$zoomedplot.labs <<- MGvar$activeplot.labs.T3
            MGvar$zoomedplot.showpoints <<- MGvar$activeplot.showpoints.T3
            MGvar$zoomedplot.pointcol <<- MGvar$activeplot.pointcol.T3
            MGvar$zoomedplot.type <<- MGvar$activeplot.type.T3
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.yaxt.T3
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.xaxt.T3
            MGvar$zoomedplot.axescol <<- MGvar$activeplot.axescol.T3
            MGvar$zoomedplot.showreg <<- MGvar$activeplot.showreg.T3
            MGvar$zoomedplot.regcol <<- MGvar$activeplot.regcol.T3
            MGvar$removedpointsactivedata <<- MGvar$removedpointsactivedata.T3
            MGvar$removedpoints <<- MGvar$removedpoints.T3
            MGvar$remindex <<- MGvar$remindex.T3
            MGvar$RemovedPoints <<- MGvar$RemovedPoints.T3
            MGvar$remAxcompindex <<- MGvar$remAxcompindex.T3
            MGvar$RemovedAxes <<- MGvar$RemovedAxes.T3
            MGvar$remPcompindex <<- MGvar$remPcompindex.T3
            MGvar$DistFunc <<- MGvar$DistFunc.T3
            MGvar$is.Metric <<- MGvar$is.Metric.T3
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkconfigure(img4, state = "active")
            tk2notetab.select(myPlottingNB, "Plot4")
            MGvar$MDSmat <<- MGvar$MDSmat.T4
            MGvar$distmat <<- MGvar$distmat.T4
            MGvar$scree.stress <<- MGvar$scree.stress.T4
            MGvar$screepoints.current <<- MGvar$screepoints.current.T4
            MGvar$screepoints.best <<- MGvar$screepoints.best.T4
            MGvar$activeplot.title <<- MGvar$activeplot.title.T4
            MGvar$indexLabeled <<- MGvar$indexLabeled.T4
            MGvar$Opt.dim <<- MGvar$Opt.dim.T4
            MGvar$MDSStress <<- MGvar$MDSStress.T4
            MGvar$dMeas <<- MGvar$dMeas.T4
            AShepdistmat <<- MGvar$distmat.T4
            AShepMDSmat <<- MGvar$MDSmat.T4
            tclvalue(PlotText) <- paste("Active Plot is Plot4")
            tkbind(img, "<Button-1>", DoNothing)
            tkbind(img2, "<Button-1>", DoNothing)
            tkbind(img3, "<Button-1>", DoNothing)
            tkbind(img4, "<Button-1>", OnPlotLeftClick)
            tkbind(img5, "<Button-1>", DoNothing)
            tkbind(img, "<B1-Motion>", DoNothing)
            tkbind(img2, "<B1-Motion>", DoNothing)
            tkbind(img3, "<B1-Motion>", DoNothing)
            tkbind(img4, "<B1-Motion>", BrushingPointMove)
            tkbind(img5, "<B1-Motion>", DoNothing)
            tkbind(img, "<ButtonRelease-1>", DoNothing)
            tkbind(img2, "<ButtonRelease-1>", DoNothing)
            tkbind(img3, "<ButtonRelease-1>", DoNothing)
            tkbind(img4, "<ButtonRelease-1>", OnRelease.Main)
            tkbind(img5, "<ButtonRelease-1>", DoNothing)
            tkconfigure(img, cursor = "arrow")
            tkconfigure(img2, cursor = "arrow")
            tkconfigure(img3, cursor = "arrow")
            tkconfigure(img4, cursor = "hand2")
            tkconfigure(img5, cursor = "arrow")
            MGvar$MDStype <<- MGvar$MDStype.T4
            MGvar$newCoords <<- MGvar$newCoords.T4
            MGvar$zoomedplot.title.show <<- MGvar$activeplot.title.show.T4
            MGvar$zoomedplot.distmeas <<- MGvar$activeplot.distmeas.T4
            MGvar$zoomedplot.xlab <<- MGvar$activeplot.xlab.T4
            MGvar$zoomedplot.ylab <<- MGvar$activeplot.ylab.T4
            MGvar$zoomedplot.bg <<- MGvar$activeplot.bg.T4
            MGvar$zoomedplot.cex <<- MGvar$activeplot.cex.T4
            MGvar$zoomedplot.labs <<- MGvar$activeplot.labs.T4
            MGvar$zoomedplot.showpoints <<- MGvar$activeplot.showpoints.T4
            MGvar$zoomedplot.pointcol <<- MGvar$activeplot.pointcol.T4
            MGvar$zoomedplot.type <<- MGvar$activeplot.type.T4
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.yaxt.T4
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.xaxt.T4
            MGvar$zoomedplot.axescol <<- MGvar$activeplot.axescol.T4
            MGvar$zoomedplot.showreg <<- MGvar$activeplot.showreg.T4
            MGvar$zoomedplot.regcol <<- MGvar$activeplot.regcol.T4
            MGvar$removedpointsactivedata <<- MGvar$removedpointsactivedata.T4
            MGvar$removedpoints <<- MGvar$removedpoints.T4
            MGvar$remindex <<- MGvar$remindex.T4
            MGvar$RemovedPoints <<- MGvar$RemovedPoints.T4
            MGvar$remAxcompindex <<- MGvar$remAxcompindex.T4
            MGvar$RemovedAxes <<- MGvar$RemovedAxes.T4
            MGvar$remPcompindex <<- MGvar$remPcompindex.T4
            MGvar$DistFunc <<- MGvar$DistFunc.T4
            MGvar$is.Metric <<- MGvar$is.Metric.T4
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkconfigure(img5, state = "active")
            tk2notetab.select(myPlottingNB, "Plot5")
            MGvar$MDSmat <<- MGvar$MDSmat.T5
            MGvar$distmat <<- MGvar$distmat.T5
            MGvar$scree.stress <<- MGvar$scree.stress.T5
            MGvar$screepoints.current <<- MGvar$screepoints.current.T5
            MGvar$screepoints.best <<- MGvar$screepoints.best.T5
            MGvar$activeplot.title <<- MGvar$activeplot.title.T5
            MGvar$indexLabeled <<- MGvar$indexLabeled.T5
            MGvar$Opt.dim <<- MGvar$Opt.dim.T5
            MGvar$MDSStress <<- MGvar$MDSStress.T5
            MGvar$dMeas <<- MGvar$dMeas.T5
            AShepdistmat <<- MGvar$distmat.T5
            AShepMDSmat <<- MGvar$MDSmat.T5
            tclvalue(PlotText) <- paste("Active Plot is Plot5")
            tkbind(img, "<Button-1>", DoNothing)
            tkbind(img2, "<Button-1>", DoNothing)
            tkbind(img3, "<Button-1>", DoNothing)
            tkbind(img4, "<Button-1>", DoNothing)
            tkbind(img5, "<Button-1>", OnPlotLeftClick)
            tkbind(img, "<B1-Motion>", DoNothing)
            tkbind(img2, "<B1-Motion>", DoNothing)
            tkbind(img3, "<B1-Motion>", DoNothing)
            tkbind(img4, "<B1-Motion>", DoNothing)
            tkbind(img5, "<B1-Motion>", BrushingPointMove)
            tkbind(img, "<ButtonRelease-1>", DoNothing)
            tkbind(img2, "<ButtonRelease-1>", DoNothing)
            tkbind(img3, "<ButtonRelease-1>", DoNothing)
            tkbind(img4, "<ButtonRelease-1>", DoNothing)
            tkbind(img5, "<ButtonRelease-1>", OnRelease.Main)
            tkconfigure(img, cursor = "arrow")
            tkconfigure(img2, cursor = "arrow")
            tkconfigure(img3, cursor = "arrow")
            tkconfigure(img4, cursor = "arrow")
            tkconfigure(img5, cursor = "hand2")
            MGvar$MDStype <<- MGvar$MDStype.T5
            MGvar$newCoords <<- MGvar$newCoords.T5
            MGvar$zoomedplot.title.show <<- MGvar$activeplot.title.show.T5
            MGvar$zoomedplot.distmeas <<- MGvar$activeplot.distmeas.T5
            MGvar$zoomedplot.xlab <<- MGvar$activeplot.xlab.T5
            MGvar$zoomedplot.ylab <<- MGvar$activeplot.ylab.T5
            MGvar$zoomedplot.bg <<- MGvar$activeplot.bg.T5
            MGvar$zoomedplot.cex <<- MGvar$activeplot.cex.T5
            MGvar$zoomedplot.labs <<- MGvar$activeplot.labs.T5
            MGvar$zoomedplot.showpoints <<- MGvar$activeplot.showpoints.T5
            MGvar$zoomedplot.pointcol <<- MGvar$activeplot.pointcol.T5
            MGvar$zoomedplot.type <<- MGvar$activeplot.type.T5
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.yaxt.T5
            MGvar$zoomedplot.yaxt <<- MGvar$activeplot.xaxt.T5
            MGvar$zoomedplot.axescol <<- MGvar$activeplot.axescol.T5
            MGvar$zoomedplot.showreg <<- MGvar$activeplot.showreg.T5
            MGvar$zoomedplot.regcol <<- MGvar$activeplot.regcol.T5
            MGvar$removedpointsactivedata <<- MGvar$removedpointsactivedata.T5
            MGvar$removedpoints <<- MGvar$removedpoints.T5
            MGvar$remindex <<- MGvar$remindex.T5
            MGvar$RemovedPoints <<- MGvar$RemovedPoints.T5
            MGvar$remAxcompindex <<- MGvar$remAxcompindex.T5
            MGvar$RemovedAxes <<- MGvar$RemovedAxes.T5
            MGvar$remPcompindex <<- MGvar$remPcompindex.T5
            MGvar$DistFunc <<- MGvar$DistFunc.T5
            MGvar$is.Metric <<- MGvar$is.Metric.T5
        }
        MGvar$activeplot.title.show <<- "yes"
        MGvar$activeplot.distmeas <<- "yes"
        MGvar$activeplot.xlab <<- ""
        MGvar$activeplot.ylab <<- ""
        MGvar$activeplot.bg <<- "white"
        MGvar$activeplot.bg.temp <<- "white"
        MGvar$activeplot.cex <<- 0.6
        MGvar$activeplot.labs <<- "yes"
        MGvar$activeplot.showpoints <<- "no"
        MGvar$activeplot.pointcol <<- "black"
        MGvar$activeplot.pointcol.temp <<- "black"
        MGvar$activeplot.type <<- "1"
        MGvar$activeplot.yaxt <<- "n"
        MGvar$activeplot.yaxt <<- "n"
        MGvar$activeplot.axescol <<- "black"
        MGvar$activeplot.axescol.temp <<- "black"
        MGvar$activeplot.showreg.T5 <- "no"
        MGvar$activeplot.regcol <- "red"
        MGvar$activeplot.regcol.temp <- "red"
        MGvar$tShepx <<- as.vector(0)
        tabplot()
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
        tkrreplot(imgscree)
        tableupdate.Remp()
        tableupdate.RemAx()
        ResetUndo()
    }
    QuitGUI = function() {
    }
    RotateandReflect <- function() {
        RAR.Rotate.Val = tclVar(0)
        RAR.Rotate.Dir = tclVar("Clockwise")
        R.A.R = tktoplevel()
        tkwm.resizable(R.A.R, "0", "0")
        tkwm.deiconify(R.A.R)
        tkwm.title(R.A.R, "Plotting Point Options")
        tkwm.geometry(R.A.R, "350x300")
        RARcanvas = tkcanvas(R.A.R, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(RARcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = R.A.R)
        RARnotebook <- tk2notebook(R.A.R, tabs = NULL)
        RAR.1 <- tk2frame(RARnotebook)
        tkadd(RARnotebook, RAR.1, text = "Rotation & Reflection")
        frameRAR1 <- tkwidget(RAR.1, "TitleFrame", text = "Rotation", 
            background = "white")
        tkplace(frameRAR1, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.47, `in` = RAR.1)
        tkplace(tklabel(frameRAR1, text = "Direction", background = "white"), 
            relx = 0.1, rely = 0.2, `in` = frameRAR1)
        RAR.spinner.Dir <- tk2spinbox(R.A.R, values = c("Anti-Clockwise", 
            "Clockwise"), textvariable = RAR.Rotate.Dir, width = 15)
        tkplace(RAR.spinner.Dir, relx = 0.6, rely = 0.2, `in` = frameRAR1)
        tkplace(tklabel(frameRAR1, text = "Angle of Rotation", 
            background = "white"), relx = 0.1, rely = 0.6, `in` = frameRAR1)
        RAR.spinner.Deg <- tk2spinbox(R.A.R, from = 0, to = 180, 
            increment = 10, textvariable = RAR.Rotate.Val, width = 15)
        tkplace(RAR.spinner.Deg, relx = 0.6, rely = 0.6, `in` = frameRAR1)
        frameRAR2 <- tkwidget(RAR.1, "TitleFrame", text = "Reflection", 
            background = "white")
        tkplace(frameRAR2, relx = 0.02, relwidth = 0.96, rely = 0.51, 
            relheight = 0.47, `in` = RAR.1)
        tkplace(tklabel(frameRAR2, text = "Rotate about X-Axis", 
            background = "white"), relx = 0.1, rely = 0.2, `in` = frameRAR2)
        tkplace(tklabel(frameRAR2, text = "Rotate about Y-Axis", 
            background = "white"), relx = 0.1, rely = 0.6, `in` = frameRAR2)
        RAR.cb.xrot <- tk2checkbutton(R.A.R)
        RARcbVx <- tclVar("0")
        tkconfigure(RAR.cb.xrot, variable = RARcbVx)
        RAR.cb.yrot <- tk2checkbutton(R.A.R)
        RARcbVy <- tclVar("0")
        tkconfigure(RAR.cb.yrot, variable = RARcbVy)
        tkplace(RAR.cb.xrot, relx = 0.75, rely = 0.2, `in` = frameRAR2)
        tkplace(RAR.cb.yrot, relx = 0.75, rely = 0.6, `in` = frameRAR2)
        tkplace(RARnotebook, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = R.A.R)
        rotate2D = function(coords, theta, direction = c("C", 
            "AC")) {
            tkconfigure(mytt, cursor = "watch")
            direction = as.character(direction)
            if (direction == "C") {
                theta = -theta
            }
            if (direction == "AC") {
            }
            radians = theta/57.2957795130823
            rotmat = matrix(c(cos(radians), -sin(radians), sin(radians), 
                cos(radians)), nrow = 2)
            rotcoords = coords %*% rotmat
            dimnames(rotcoords) = dimnames(coords)
            MGvar$MDSmat <<- rotcoords
            ActiveCoordMat()
            tabplot()
            tkfocus(mytt)
            tkconfigure(mytt, cursor = "arrow")
        }
        reflect2D = function(coords, xaxis = c("yes", "no"), 
            yaxis = c("yes", "no")) {
            tkconfigure(mytt, cursor = "watch")
            refcoords = matrix(0, nrow(coords), ncol = 2)
            xaxis = as.character(xaxis)
            yaxis = as.character(yaxis)
            if (yaxis == "yes") {
                refcoords[, 1] = -coords[, 1]
            }
            else {
                refcoords[, 1] = coords[, 1]
            }
            if (xaxis == "yes") {
                refcoords[, 2] = -coords[, 2]
            }
            else {
                refcoords[, 2] = coords[, 2]
            }
            dimnames(refcoords) = dimnames(coords)
            MGvar$MDSmat <<- refcoords
            ActiveCoordMat()
            tabplot()
            tkfocus(mytt)
            tkconfigure(mytt, cursor = "arrow")
        }
        OnOK.RAR <- function() {
            RotDir = tclvalue(RAR.Rotate.Dir)
            RotDeg = tclvalue(RAR.Rotate.Val)
            XTVal = as.character(tclvalue(RARcbVx))
            YTVal = as.character(tclvalue(RARcbVy))
            RotDeg = as.numeric(RotDeg)
            RotDir = as.character(RotDir)
            if (RotDeg == 0 && XTVal == 0 && YTVal == 0) {
                tkdestroy(R.A.R)
            }
            else {
                if (nrow(MGvar$MDSmat) == 1 && ncol(MGvar$MDSmat) == 
                  1) {
                  tkmessageBox(message = "No Active Coordinates. Perform MDS procedure!", 
                    icon = "error")
                }
                else {
                  RAR.Option = tkmessageBox(message = "Are you sure you want to confirm these changes?", 
                    icon = "question", type = "yesno", default = "yes")
                  RAR.Option = as.character(RAR.Option)
                  if (RAR.Option == "yes") {
                    if (RotDir == "Clockwise") {
                      rotate2D(MGvar$MDSmat, RotDeg, direction = "C")
                    }
                    if (RotDir == "Anti-Clockwise") {
                      rotate2D(MGvar$MDSmat, RotDeg, direction = "AC")
                    }
                    if (XTVal == "0") {
                      XT = "no"
                    }
                    else {
                      XT = "yes"
                    }
                    if (YTVal == "0") {
                      YT = "no"
                    }
                    else {
                      YT = "yes"
                    }
                    reflect2D(MGvar$MDSmat, xaxis = XT, yaxis = YT)
                    tkdestroy(R.A.R)
                  }
                  if (RAR.Option == "no") {
                    tkmessageBox(message = "No changes have been made")
                  }
                }
            }
        }
        OnCancel.RAR <- function() {
            tkdestroy(R.A.R)
        }
        tkplace(tkbutton(R.A.R, text = "OK", width = 15, command = function() OnOK.RAR()), 
            relx = 0.15, rely = 0.88, `in` = R.A.R)
        tkplace(tkbutton(R.A.R, text = "Cancel", width = 15, 
            command = function() OnCancel.RAR()), relx = 0.55, 
            rely = 0.88, `in` = R.A.R)
        tkfocus(R.A.R)
        tkbind(R.A.R, "<Return>", OnOK.RAR)
        tkwait.window(R.A.R)
    }
    MDSOptions <- function() {
        MDS.dim = tclVar(MGvar$MDS.dimensions)
        MDSops = tktoplevel()
        tkwm.resizable(MDSops, "0", "0")
        tkwm.deiconify(MDSops)
        tkwm.title(MDSops, "MDS Options")
        tkwm.geometry(MDSops, "400x300")
        MDSOcanvas = tkcanvas(MDSops, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(MDSOcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = MDSops)
        MDSOpsNB <- tk2notebook(MDSops, tabs = NULL)
        Dimensions <- tk2frame(MDSOpsNB)
        tkadd(MDSOpsNB, Dimensions, text = "Dimensions")
        frameDims <- tkwidget(Dimensions, "TitleFrame", text = "Choice of Dimension", 
            background = "white")
        tkplace(frameDims, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = Dimensions)
        tkplace(tklabel(frameDims, text = "Choose Number of Dimensions for MDS", 
            background = "white"), relx = 0.18, rely = 0.1, `in` = frameDims)
        dims = seq(1:MGvar$maxdims)
        Dim.ComboBox <- tkwidget(MDSops, "ComboBox", editable = FALSE, 
            values = dims, width = 12, textvariable = MDS.dim)
        tkplace(Dim.ComboBox, relx = 0.36, rely = 0.3, `in` = frameDims)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frameDims, text = "Please note that only Dimensions 1,2 and 3 can be plotted.\nMDS output with 4 dimensions or more will only be in\ncoordinate format. Otherwise a subset of any 2 or 3 dimensions\nwill be available for plotting.", 
            font = fontsmall), relx = 0.04, rely = 0.485, `in` = frameDims)
        ChangeD <- function() {
            MGvar$temp.dims <<- as.numeric(tclvalue(MDS.dim))
            if (MGvar$temp.dims == 1) {
                tkmessageBox(message = paste("You have chosen to perform MDS in only 1 Dimension. Please be warned that while Unidimensional Scaling is possible it is NOT recommended. This is due to Unidimensional Scaling being particularly prone to finding local minima rather than global minima. It is highly recommended that the number of dimensions be increased!"))
            }
            MGvar$MDS.dimensions <<- MGvar$temp.dims
            tclvalue(dimText) <<- paste("p =", MGvar$MDS.dimensions)
        }
        tkplace(tkbutton(MDSops, text = "Change", width = 15, 
            command = function() ChangeD()), relx = 0.325, rely = 0.82, 
            `in` = frameDims)
        SConfig <- tkframe(MDSOpsNB)
        tkadd(MDSOpsNB, SConfig, text = "Starting Configurations")
        frameSConf <- tkwidget(SConfig, "TitleFrame", text = "Starting Configuration Options", 
            background = "white")
        tkplace(frameSConf, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = SConfig)
        tkplace(tklabel(frameSConf, text = "Many of the MDS procedures in the MDS-GUI make use of\na set of coordinates as the starting configuration. The default\nstarting configuration is in the form of the result of Classical\nScaling. Please choose your desired starting configuration.", 
            font = fontsmall), relx = 0.06, rely = 0.1, `in` = frameSConf)
        tkplace(tklabel(frameSConf, text = "Result of Classic Scaling", 
            background = "white"), relx = 0.18, rely = 0.42, 
            `in` = frameSConf)
        MDSops.UseCS.RB <- tk2radiobutton(MDSops)
        tkconfigure(MDSops.UseCS.RB, variable = MGvar$MDSops.startconfig, 
            value = "ClasScal")
        tkplace(MDSops.UseCS.RB, relx = 0.1, rely = 0.42, `in` = frameSConf)
        tkplace(tklabel(frameSConf, text = "Random Configuration", 
            background = "white"), relx = 0.18, rely = 0.54, 
            `in` = frameSConf)
        MDSops.UseRC.RB <- tk2radiobutton(MDSops)
        tkconfigure(MDSops.UseRC.RB, variable = MGvar$MDSops.startconfig, 
            value = "RandConf")
        tkplace(MDSops.UseRC.RB, relx = 0.1, rely = 0.54, `in` = frameSConf)
        tkplace(tklabel(frameSConf, text = "Use Existing Configuration in:", 
            background = "white"), relx = 0.18, rely = 0.66, 
            `in` = frameSConf)
        MDSops.UseEC.RB <- tk2radiobutton(MDSops)
        tkconfigure(MDSops.UseEC.RB, variable = MGvar$MDSops.startconfig, 
            value = "ExistConf")
        tkplace(MDSops.UseEC.RB, relx = 0.1, rely = 0.66, `in` = frameSConf)
        MDSops.UseEC.ComboBox <- tkwidget(MDSops, "ComboBox", 
            editable = FALSE, values = c("Plot1", "Plot2", "Plot3", 
                "Plot4", "Plot5"), width = 8, textvariable = MGvar$MDSops.UseEC.Plot)
        tkplace(MDSops.UseEC.ComboBox, relx = 0.7, rely = 0.66, 
            `in` = frameSConf)
        ChangeC <- function() {
            tclvalue(MGvar$MDSops.startconfig) <<- tclvalue(MGvar$MDSops.startconfig)
            if (as.character(tclvalue(MGvar$MDSops.startconfig)) == 
                "ExistConf") {
                MGvar$startconfig.UseEC.plot <<- tclvalue(MGvar$MDSops.UseEC.Plot)
            }
        }
        tkplace(tkbutton(MDSops, text = "Change", width = 15, 
            command = function() ChangeC()), relx = 0.325, rely = 0.82, 
            `in` = frameSConf)
        StressCalc <- tkframe(MDSOpsNB)
        tkadd(MDSOpsNB, StressCalc, text = "Stress")
        frameSC <- tkwidget(StressCalc, "TitleFrame", text = "Choice of Stress Calculation", 
            background = "white")
        tkplace(frameSC, relx = 0.05, relwidth = 0.9, rely = 0.02, 
            relheight = 0.9, `in` = StressCalc)
        tkplace(tklabel(frameSC, text = "A number of options exist for the indication of goodness\nof fit for an MDS configuration. This measure is referred\nto as 'Stress'. Select the Stress you would like to use.", 
            font = fontsmall), relx = 0.06, rely = 0.1, `in` = frameSC)
        tkplace(tklabel(frameSC, text = "Stress Calculation Method", 
            background = "white"), relx = 0.27, rely = 0.4, `in` = frameSC)
        SC.ComboBox <- tkwidget(MDSops, "ComboBox", editable = FALSE, 
            values = c("Norm. Raw Stress", "STRESS1", "STRESS2", 
                "Correlation Coefficient"), width = 20, height = 4, 
            textvariable = MGvar$MDSops.StressCalc)
        tkplace(SC.ComboBox, relx = 0.28, rely = 0.55, `in` = frameSC)
        ChangeSt <- function() {
            Calc = as.character(tclvalue(MGvar$MDSops.StressCalc))
            if (Calc == "Norm. Raw Stress") {
                tclvalue(MGvar$MDSops.StressCalc) <<- "Norm. Raw Stress"
                MGvar$StCalc <<- "NormRawStress"
            }
            if (Calc == "STRESS1") {
                tclvalue(MGvar$MDSops.StressCalc) <<- "STRESS1"
                MGvar$StCalc <<- "Stress-1"
            }
            if (Calc == "STRESS2") {
                tclvalue(MGvar$MDSops.StressCalc) <<- "STRESS2"
                MGvar$StCalc <<- "Stress-2"
            }
            if (Calc == "Correlation Coefficient") {
                tclvalue(MGvar$MDSops.StressCalc) <<- "Correlation Coefficient"
                MGvar$StCalc <<- "CorCoef"
            }
            StressUpdate()
            tableupdate()
        }
        tkplace(tkbutton(MDSops, text = "Change", width = 15, 
            command = function() ChangeSt()), relx = 0.325, rely = 0.82, 
            `in` = frameSC)
        tkplace(MDSOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = MDSops)
        OnOK.Ops <- function() {
            ChangeD()
            ChangeC()
            ChangeSt()
            tkdestroy(MDSops)
        }
        OnCancel.Ops <- function() {
            tkdestroy(MDSops)
        }
        tkplace(tkbutton(MDSops, text = "OK", width = 15, command = function() OnOK.Ops()), 
            relx = 0.15, rely = 0.9, `in` = MDSops)
        tkplace(tkbutton(MDSops, text = "Cancel", width = 15, 
            command = function() OnCancel.Ops()), relx = 0.55, 
            rely = 0.9, `in` = MDSops)
        tkfocus(MDSops)
        tkbind(MDSops, "<Return>", OnOK.Ops)
        tkwait.window(MDSops)
    }
    LargeDimensions <- function() {
        LargeDim <- tktoplevel()
        tkwm.resizable(LargeDim, "0", "0")
        tkwm.deiconify(LargeDim)
        tkwm.title(LargeDim, "Dimension Options")
        tkwm.geometry(LargeDim, "450x350")
        LDcanvas <- tkcanvas(LargeDim, width = "450", height = "350", 
            bg = col.sec)
        tkplace(LDcanvas, `in` = LargeDim)
        frameLD <- tkwidget(LargeDim, "TitleFrame", text = "Four or More Dimension Options", 
            background = "white")
        tkplace(frameLD, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            relheight = 0.8, `in` = LargeDim)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frameLD, text = "You are attempting to perform an MDS procedure with more than 3 dimensions.\nYou may now choose to plot any two of the resultant dimensions in 2D or 3D\nformat. Any resulting 2D configuration will be output to the Active Plot.", 
            font = fontsmall), relx = 0.02, rely = 0.08, `in` = frameLD)
        LargeDimNB <- tk2notebook(LargeDim, tabs = NULL)
        LD.t1 <- tk2frame(LargeDimNB)
        tkadd(LargeDimNB, LD.t1, text = "    2D    ")
        frame2D <- tkframe(LD.t1, relief = "groove", borderwidth = 2, 
            background = "white")
        tkplace(frame2D, relx = 0.005, relwidth = 0.99, rely = 0.005, 
            relheight = 0.99, `in` = LD.t1)
        tkplace(tklabel(frame2D, text = "Plotting X Dimension", 
            background = "white"), relx = 0.1, rely = 0.12, `in` = frame2D)
        tkplace(tklabel(frame2D, text = "Dimension", background = "white"), 
            relx = 0.6, rely = 0.12, `in` = frame2D)
        dimensions = c()
        for (i in 0:(MGvar$MDS.dimensions - 1)) {
            newdim = paste(MGvar$MDS.dimensions - i)
            dimensions = c(dimensions, newdim)
        }
        XDim.spin <- tk2spinbox(LargeDim, values = dimensions, 
            width = 2)
        XDim.val <- tclVar("1")
        tkconfigure(XDim.spin, textvariable = XDim.val)
        tkplace(XDim.spin, relx = 0.8, rely = 0.12, `in` = frame2D)
        tkplace(tklabel(frame2D, text = "Plotting Y Dimension", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frame2D)
        tkplace(tklabel(frame2D, text = "Dimension", background = "white"), 
            relx = 0.6, rely = 0.4, `in` = frame2D)
        YDim.spin <- tk2spinbox(LargeDim, values = dimensions, 
            width = 2)
        YDim.val <- tclVar("2")
        tkconfigure(YDim.spin, textvariable = YDim.val)
        tkplace(YDim.spin, relx = 0.8, rely = 0.4, `in` = frame2D)
        On.2DPlot <- function() {
            XDim = as.numeric(tclvalue(XDim.val))
            MGvar$PlottingDimX <<- XDim
            YDim = as.numeric(tclvalue(YDim.val))
            MGvar$PlottingDimY <<- YDim
            MGvar$MDSmat <<- matrix(0, nrow = nrow(MGvar$MDSmat.LD), 
                ncol = 2)
            MGvar$MDSmat[, 1] <<- MGvar$MDSmat.LD[, MGvar$PlottingDimX]
            MGvar$MDSmat[, 2] <<- MGvar$MDSmat.LD[, MGvar$PlottingDimY]
            rownames(MGvar$MDSmat) <<- rownames(MGvar$activedata)
            tkdestroy(LargeDim)
        }
        tkplace(tkbutton(LargeDim, text = "Plot 2D", width = 15, 
            command = function() On.2DPlot()), relx = 0.35, rely = 0.75, 
            `in` = frame2D)
        LD.t2 <- tk2frame(LargeDimNB)
        tkadd(LargeDimNB, LD.t2, text = "    3D    ")
        frame3D <- tkframe(LD.t2, relief = "groove", borderwidth = 2, 
            background = "white")
        tkplace(frame3D, relx = 0.005, relwidth = 0.99, rely = 0.005, 
            relheight = 0.99, `in` = LD.t2)
        tkplace(tklabel(frame3D, text = "Plotting X Dimension", 
            background = "white"), relx = 0.1, rely = 0.12, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "Dimension", background = "white"), 
            relx = 0.6, rely = 0.12, `in` = frame3D)
        XDim3D.spin <- tk2spinbox(LargeDim, values = dimensions, 
            width = 2)
        XDim3D.val <- tclVar("1")
        tkconfigure(XDim3D.spin, textvariable = XDim3D.val)
        tkplace(XDim3D.spin, relx = 0.8, rely = 0.12, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "Plotting Y Dimension", 
            background = "white"), relx = 0.1, rely = 0.3, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "Dimension", background = "white"), 
            relx = 0.6, rely = 0.3, `in` = frame3D)
        YDim3D.spin <- tk2spinbox(LargeDim, values = dimensions, 
            width = 2)
        YDim3D.val <- tclVar("2")
        tkconfigure(YDim3D.spin, textvariable = YDim3D.val)
        tkplace(YDim3D.spin, relx = 0.8, rely = 0.3, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "Plotting Z Dimension", 
            background = "white"), relx = 0.1, rely = 0.48, `in` = frame3D)
        tkplace(tklabel(frame3D, text = "Dimension", background = "white"), 
            relx = 0.6, rely = 0.48, `in` = frame3D)
        ZDim3D.spin <- tk2spinbox(LargeDim, values = dimensions, 
            width = 2)
        ZDim3D.val <- tclVar("3")
        tkconfigure(ZDim3D.spin, textvariable = ZDim3D.val)
        tkplace(ZDim3D.spin, relx = 0.8, rely = 0.48, `in` = frame3D)
        On.3DPlot <- function() {
            MGvar$tempdims <<- 3
            XDim3D = as.numeric(tclvalue(XDim3D.val))
            MGvar$PlottingDimX.3D <<- XDim3D
            YDim3D = as.numeric(tclvalue(YDim3D.val))
            MGvar$PlottingDimY.3D <<- YDim3D
            ZDim3D = as.numeric(tclvalue(ZDim3D.val))
            MGvar$PlottingDimZ.3D <<- ZDim3D
            MGvar$MDSmat <<- matrix(0, nrow = nrow(MGvar$MDSmat.LD), 
                ncol = 3)
            MGvar$MDSmat[, 1] <<- MGvar$MDSmat.LD[, MGvar$PlottingDimX.3D]
            MGvar$MDSmat[, 2] <<- MGvar$MDSmat.LD[, MGvar$PlottingDimY.3D]
            MGvar$MDSmat[, 3] <<- MGvar$MDSmat.LD[, MGvar$PlottingDimZ.3D]
            rownames(MGvar$MDSmat) <<- rownames(MGvar$distmat)
            threedimplotting()
            tkdestroy(LargeDim)
        }
        tkplace(tkbutton(LargeDim, text = "Plot 3D", width = 15, 
            command = function() On.3DPlot()), relx = 0.35, rely = 0.75, 
            `in` = frame3D)
        tkplace(LargeDimNB, relx = 0.03, relwidth = 0.94, rely = 0.28, 
            relheight = 0.69, `in` = frameLD)
        On.Coords <- function() {
            MGvar$MDSmat <<- MGvar$MDSmat.LD
            ShowActiveMDSmat()
            tkdestroy(LargeDim)
        }
        On.Cancel <- function() {
            tkdestroy(LargeDim)
        }
        tkplace(tkbutton(LargeDim, text = "Just Coordinates", 
            width = 15, command = function() On.Coords()), relx = 0.15, 
            rely = 0.87, `in` = LargeDim)
        tkplace(tkbutton(LargeDim, text = "Cancel", width = 15, 
            command = function() On.Cancel()), relx = 0.6, rely = 0.87, 
            `in` = LargeDim)
        tkfocus(LargeDim)
        tkwait.window(LargeDim)
    }
    ConfPlotOptions <- function() {
        CPlotOps = tktoplevel()
        tkwm.resizable(CPlotOps, "0", "0")
        tkwm.deiconify(CPlotOps)
        tkwm.title(CPlotOps, "MDS Configuration Plotting Options")
        tkwm.geometry(CPlotOps, "350x400")
        CPlotcanvas = tkcanvas(CPlotOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(CPlotcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = CPlotOps)
        ConfPlotcanvas = tkcanvas(CPlotOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(ConfPlotcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = CPlotOps)
        CPlotOpsNB <- tk2notebook(CPlotOps, tabs = NULL)
        ConfGen <- tk2frame(CPlotOpsNB)
        tkadd(CPlotOpsNB, ConfGen, text = "General")
        frameCGen <- tkwidget(ConfGen, "TitleFrame", text = "MDS Configuration General", 
            background = "white")
        tkplace(frameCGen, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ConfGen)
        tkplace(tklabel(frameCGen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameCGen)
        Conf.Main.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Main.CB, variable = MGvar$Conf.Main.val)
        tkplace(Conf.Main.CB, relx = 0.7, rely = 0.1, `in` = frameCGen)
        tkplace(tklabel(frameCGen, text = "Display Distance Measure", 
            background = "white"), relx = 0.1, rely = 0.22, `in` = frameCGen)
        Conf.Dist.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Dist.CB, variable = MGvar$Conf.Dist.val)
        tkplace(Conf.Dist.CB, relx = 0.7, rely = 0.22, `in` = frameCGen)
        tkplace(tklabel(frameCGen, text = "Display Legend", background = "white"), 
            relx = 0.1, rely = 0.34, `in` = frameCGen)
        Conf.Leg.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Leg.CB, variable = MGvar$Conf.Leg.val)
        tkplace(Conf.Leg.CB, relx = 0.7, rely = 0.34, `in` = frameCGen)
        tkplace(tklabel(frameCGen, text = "Provide Label for Y-Axis", 
            background = "white"), relx = 0.1, rely = 0.46, `in` = frameCGen)
        Conf.Ylab.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Ylab.CB, variable = MGvar$Conf.Ylab.val)
        tkplace(Conf.Ylab.CB, relx = 0.7, rely = 0.46, `in` = frameCGen)
        tkplace(tklabel(frameCGen, text = "Provide Label for X-Axis", 
            background = "white"), relx = 0.1, rely = 0.58, `in` = frameCGen)
        Conf.Xlab.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Xlab.CB, variable = MGvar$Conf.Xlab.val)
        tkplace(Conf.Xlab.CB, relx = 0.7, rely = 0.58, `in` = frameCGen)
        tkplace(tklabel(frameCGen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameCGen)
        ChangeColBG <- function() {
            MGvar$activeplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$activeplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$activeplot.bg.temp) > 0) 
                tkconfigure(ConfColBG, bg = MGvar$activeplot.bg.temp)
        }
        ConfColBG <- tkbutton(CPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$activeplot.bg, command = function() ChangeColBG())
        tkplace(ConfColBG, relx = 0.7, rely = 0.7, `in` = frameCGen)
        ChangeGen <- function() {
            MainC = as.character(tclvalue(MGvar$Conf.Main.val))
            if (MainC == "0") {
                MGvar$activeplot.title.show <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.title.show.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.title.show.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.title.show.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.title.show.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.title.show.T5 <<- "no"
                }
                MGvar$zoomedplot.title.show <<- "no"
                tclvalue(MGvar$Zoom.Main.val) <<- 0
                tclvalue(MGvar$Conf.Main.val) <<- 0
            }
            if (MainC == "1") {
                MGvar$activeplot.title.show <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.title.show.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.title.show.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.title.show.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.title.show.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.title.show.T5 <<- "yes"
                }
                MGvar$zoomedplot.title.show <<- "yes"
                tclvalue(MGvar$Zoom.Main.val) <<- 1
            }
            DistC = as.character(tclvalue(MGvar$Conf.Dist.val))
            if (DistC == "0") {
                MGvar$activeplot.distmeas <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.distmeas.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.distmeas.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.distmeas.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.distmeas.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.distmeas.T5 <<- "no"
                }
                MGvar$zoomedplot.distmeas <<- "no"
                tclvalue(MGvar$Zoom.Dist.val) <<- 0
                tclvalue(MGvar$Conf.Dist.val) <<- 0
            }
            if (DistC == "1") {
                MGvar$activeplot.distmeas <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.distmeas.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.distmeas.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.distmeas.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.distmeas.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.distmeas.T5 <<- "yes"
                }
                MGvar$zoomedplot.distmeas <<- "yes"
                tclvalue(MGvar$Zoom.Dist.val) <<- 1
                tclvalue(MGvar$Conf.Dist.val) <<- 1
            }
            LegC = as.character(tclvalue(MGvar$Conf.Leg.val))
            if (LegC == "1") {
                MGvar$activeplot.showleg <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showleg.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showleg.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showleg.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showleg.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showleg.T5 <<- "yes"
                }
                MGvar$zoomedplot.showleg <<- "yes"
                tclvalue(MGvar$Zoom.Leg.val) <<- 1
                tclvalue(MGvar$Conf.Leg.val) <<- 1
            }
            if (LegC == "0") {
                MGvar$activeplot.showleg <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showleg.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showleg.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showleg.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showleg.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showleg.T5 <<- "no"
                }
                MGvar$zoomedplot.showleg <<- "no"
                tclvalue(MGvar$Zoom.Leg.val) <<- 0
                tclvalue(MGvar$Conf.Leg.val) <<- 0
            }
            YlabC = as.character(tclvalue(MGvar$Conf.Ylab.val))
            if (YlabC == "0") {
                MGvar$activeplot.ylab <<- ""
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.ylab.T1 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.ylab.T2 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.ylab.T3 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.ylab.T4 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.ylab.T5 <<- ""
                }
                MGvar$zoomedplot.ylab <<- ""
                tclvalue(MGvar$Zoom.Ylab.val) <<- 0
            }
            if (YlabC == "1") {
                tclvalue(MGvar$Conf.Ylab.val) <<- 1
                tclvalue(MGvar$Zoom.Ylab.val) <<- 1
                CYlabtt = tktoplevel()
                tkwm.resizable(CYlabtt, "0", "0")
                tkwm.deiconify(CYlabtt)
                tkwm.title(CYlabtt, "Configuration Plot Y Label")
                tkwm.geometry(CYlabtt, "310x100")
                CYlabcanvas = tkcanvas(CYlabtt, width = "310", 
                  height = "100", bg = col.sec)
                tkplace(CYlabcanvas, `in` = CYlabtt)
                frameCY <- tkwidget(CYlabtt, "TitleFrame", text = "Y Label Input", 
                  background = "white")
                tkplace(frameCY, relx = 0.01, relwidth = 0.98, 
                  rely = 0.01, relheight = 0.98, `in` = CYlabtt)
                tkplace(tklabel(frameCY, text = "Enter your Label for Y-Axis", 
                  background = "white"), relx = 0.05, rely = 0.25, 
                  `in` = frameCY)
                YInput = tclVar("")
                Ytext = tkentry(CYlabtt, width = 15, textvariable = YInput)
                tkplace(Ytext, relx = 0.6, rely = 0.25, `in` = frameCY)
                On.Enter <- function() {
                  tkdestroy(CYlabtt)
                  YL = as.character(tclvalue(YInput))
                  MGvar$activeplot.ylab <<- YL
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                    MGvar$activeplot.ylab.T1 <<- YL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                    MGvar$activeplot.ylab.T2 <<- YL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                    MGvar$activeplot.ylab.T3 <<- YL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                    MGvar$activeplot.ylab.T4 <<- YL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                    MGvar$activeplot.ylab.T5 <<- YL
                  }
                  MGvar$zoomedplot.ylab <<- YL
                }
                tkplace(tkbutton(CYlabtt, width = 15, text = "Enter", 
                  command = function() On.Enter()), relx = 0.32, 
                  rely = 0.65, `in` = frameCY)
                tkbind(CYlabtt, "<Return>", On.Enter)
                tkwait.window(CYlabtt)
            }
            XlabC = as.character(tclvalue(MGvar$Conf.Xlab.val))
            if (XlabC == "0") {
                MGvar$activeplot.xlab <<- ""
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.xlab.T1 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.xlab.T2 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.xlab.T3 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.xlab.T4 <<- ""
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.xlab.T5 <<- ""
                }
                MGvar$zoomedplot.xlab <<- ""
                tclvalue(MGvar$Zoom.Xlab.val) <<- 0
            }
            if (XlabC == "1") {
                tclvalue(MGvar$Conf.Xlab.val) <<- 1
                tclvalue(MGvar$Zoom.Xlab.val) <<- 1
                CXlabtt = tktoplevel()
                tkwm.resizable(CXlabtt, "0", "0")
                tkwm.deiconify(CXlabtt)
                tkwm.title(CXlabtt, "Configuration Plot X Label")
                tkwm.geometry(CXlabtt, "310x100")
                CXlabcanvas = tkcanvas(CXlabtt, width = "310", 
                  height = "100", bg = col.sec)
                tkplace(CXlabcanvas, `in` = CXlabtt)
                frameCX <- tkwidget(CXlabtt, "TitleFrame", text = "X Label Input", 
                  background = "white")
                tkplace(frameCX, relx = 0.01, relwidth = 0.98, 
                  rely = 0.01, relheight = 0.98, `in` = CXlabtt)
                tkplace(tklabel(frameCX, text = "Enter your Label for X-Axis", 
                  background = "white"), relx = 0.05, rely = 0.25, 
                  `in` = frameCX)
                XInput = tclVar("")
                Xtext = tkentry(CXlabtt, width = 15, textvariable = XInput)
                tkplace(Xtext, relx = 0.6, rely = 0.25, `in` = frameCX)
                On.EnterX <- function() {
                  tkdestroy(CXlabtt)
                  XL = as.character(tclvalue(XInput))
                  MGvar$activeplot.xlab <<- XL
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                    MGvar$activeplot.xlab.T1 <<- XL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                    MGvar$activeplot.xlab.T2 <<- XL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                    MGvar$activeplot.xlab.T3 <<- XL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                    MGvar$activeplot.xlab.T4 <<- XL
                  }
                  if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                    MGvar$activeplot.xlab.T5 <<- XL
                  }
                  MGvar$zoomedplot.xlab <<- XL
                }
                tkplace(tkbutton(CXlabtt, width = 15, text = "Enter", 
                  command = function() On.EnterX()), relx = 0.32, 
                  rely = 0.65, `in` = frameCX)
                tkbind(CXlabtt, "<Return>", On.EnterX)
                tkwait.window(CXlabtt)
            }
            MGvar$activeplot.bg <<- MGvar$activeplot.bg.temp
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.bg.T1 <<- MGvar$activeplot.bg.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.bg.T2 <<- MGvar$activeplot.bg.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.bg.T3 <<- MGvar$activeplot.bg.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.bg.T4 <<- MGvar$activeplot.bg.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.bg.T5 <<- MGvar$activeplot.bg.temp
            }
            MGvar$zoomedplot.bg.temp <<- MGvar$activeplot.bg.temp
            MGvar$zoomedplot.bg <<- MGvar$zoomedplot.bg.temp
            tabplot()
        }
        tkplace(tkbutton(CPlotOps, text = "Change", width = 15, 
            command = function() ChangeGen()), relx = 0.325, 
            rely = 0.85, `in` = frameCGen)
        ConfPoints <- tk2frame(CPlotOpsNB)
        tkadd(CPlotOpsNB, ConfPoints, text = "Points")
        frameCPoints <- tkwidget(ConfPoints, "TitleFrame", text = "MDS Configuration Points", 
            background = "white")
        tkplace(frameCPoints, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ConfPoints)
        tkplace(tklabel(frameCPoints, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameCPoints)
        Conf.Points.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Points.CB, variable = MGvar$Conf.Points.val)
        tkplace(Conf.Points.CB, relx = 0.7, rely = 0.1, `in` = frameCPoints)
        tkplace(tklabel(frameCPoints, text = "Display Point Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameCPoints)
        Conf.Labels.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.Labels.CB, variable = MGvar$Conf.Labels.val)
        tkplace(Conf.Labels.CB, relx = 0.7, rely = 0.25, `in` = frameCPoints)
        tkplace(tklabel(frameCPoints, text = "Point Size", background = "white"), 
            relx = 0.1, rely = 0.4, `in` = frameCPoints)
        Conf.PS.var <- tclVar(MGvar$activeplot.cex)
        Conf.PS.spin <- tk2spinbox(CPlotOps, from = 0.1, to = 2, 
            increment = 0.1, width = 12)
        tkconfigure(Conf.PS.spin, textvariable = Conf.PS.var)
        tkplace(Conf.PS.spin, relx = 0.6, rely = 0.4, `in` = frameCPoints)
        tkplace(tklabel(frameCPoints, text = "Point Type", background = "white"), 
            relx = 0.1, rely = 0.55, `in` = frameCPoints)
        shapes <- c("Empty Boxes", "Filled Boxes", "Crosses", 
            "Empty Triangles", "Filled Triangles", "Filled Circles", 
            "Empty Circles")
        Conf.PT.ComboBox <- tkwidget(CPlotOps, "ComboBox", editable = FALSE, 
            values = shapes, width = 12, textvariable = MGvar$Conf.PT.var)
        tkplace(Conf.PT.ComboBox, relx = 0.6, rely = 0.55, `in` = frameCPoints)
        tkplace(tklabel(frameCPoints, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameCPoints)
        ChangeColPoints <- function() {
            MGvar$activeplot.pointcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$activeplot.pointcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$activeplot.pointcol.temp) > 0) 
                tkconfigure(ConfColPoints, bg = MGvar$activeplot.pointcol.temp)
            MGvar$activeplot.pointcol <<- MGvar$activeplot.pointcol.temp
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.pointcol.T1 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.pointcol.T2 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.pointcol.T3 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.pointcol.T4 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.pointcol.T5 <<- MGvar$activeplot.pointcol.temp
            }
            PointColInitialise()
        }
        ConfColPoints <- tkbutton(CPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$activeplot.pointcol, command = function() ChangeColPoints())
        tkplace(ConfColPoints, relx = 0.7, rely = 0.7, `in` = frameCPoints)
        ChangePoints <- function() {
            DispP = as.character(tclvalue(MGvar$Conf.Points.val))
            if (DispP == "0") {
                MGvar$activeplot.showpoints <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showpoints.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showpoints.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showpoints.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showpoints.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showpoints.T5 <<- "no"
                }
                MGvar$zoomedplot.showpoints <<- "no"
                tclvalue(MGvar$Zoom.Points.val) <<- 0
            }
            if (DispP == "1") {
                MGvar$activeplot.showpoints <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showpoints.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showpoints.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showpoints.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showpoints.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showpoints.T5 <<- "yes"
                }
                MGvar$zoomedplot.showpoints <<- "yes"
                tclvalue(MGvar$Zoom.Points.val) <<- 1
                tclvalue(MGvar$Conf.Points.val) <<- 1
            }
            DispL = as.character(tclvalue(MGvar$Conf.Labels.val))
            if (DispL == "0") {
                MGvar$activeplot.labs <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.labs.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.labs.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.labs.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.labs.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.labs.T5 <<- "no"
                }
                MGvar$zoomedplot.labs <<- "no"
                tclvalue(MGvar$Conf.Labels.val) <<- 0
                tclvalue(MGvar$Zoom.Labels.val) <<- 0
            }
            if (DispL == "1") {
                MGvar$activeplot.labs <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.labs.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.labs.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.labs.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.labs.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.labs.T5 <<- "yes"
                }
                MGvar$zoomedplot.labs <<- "yes"
                tclvalue(MGvar$Zoom.Labels.val) <<- 1
            }
            SizeP = as.numeric(tclvalue(Conf.PS.var))
            MGvar$activeplot.cex <<- SizeP
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.cex.T1 <<- SizeP
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.cex.T2 <<- SizeP
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.cex.T3 <<- SizeP
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.cex.T4 <<- SizeP
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.cex.T5 <<- SizeP
            }
            TypeP = as.character(tclvalue(MGvar$Conf.PT.var))
            if (TypeP == "Empty Circles") {
                MGvar$activeplot.type <<- 1
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 1
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 1
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 1
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 1
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 1
                }
                MGvar$zoomedplot.type <<- 1
                tclvalue(MGvar$Conf.PT.var) <<- "Empty Circles"
                tclvalue(MGvar$Zoom.PT.var) <<- "Empty Circles"
            }
            if (TypeP == "Filled Circles") {
                MGvar$activeplot.type <<- 16
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 16
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 16
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 16
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 16
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 16
                }
                MGvar$zoomedplot.type <<- 16
                tclvalue(MGvar$Conf.PT.var) <<- "Filled Circles"
                tclvalue(MGvar$Zoom.PT.var) <<- "Filled Circles"
            }
            if (TypeP == "Empty Boxes") {
                MGvar$activeplot.type <<- 22
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 22
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 22
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 22
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 22
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 22
                }
                MGvar$zoomedplot.type <<- 22
                tclvalue(MGvar$Conf.PT.var) <<- "Empty Boxes"
                tclvalue(MGvar$Zoom.PT.var) <<- "Empty Boxes"
            }
            if (TypeP == "Filled Boxes") {
                MGvar$activeplot.type <<- 15
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 15
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 15
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 15
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 15
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 15
                }
                MGvar$zoomedplot.type <<- 15
                tclvalue(MGvar$Conf.PT.var) <<- "Filled Boxes"
                tclvalue(MGvar$Zoom.PT.var) <<- "Filled Boces"
            }
            if (TypeP == "Crosses") {
                MGvar$activeplot.type <<- 4
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 4
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 4
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 4
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 4
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 4
                }
                MGvar$zoomedplot.type <<- 4
                tclvalue(MGvar$Conf.PT.var) <<- "Crosses"
                tclvalue(MGvar$Zoom.PT.var) <<- "Crosses"
            }
            if (TypeP == "Empty Triangles") {
                MGvar$activeplot.type <<- 24
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 24
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 24
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 24
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 24
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 24
                }
                MGvar$zoomedplot.type <<- 24
                tclvalue(MGvar$Conf.PT.var) <<- "Empty Triangles"
                tclvalue(MGvar$Zoom.PT.var) <<- "Empty Triangles"
            }
            if (TypeP == "Filled Triangles") {
                MGvar$activeplot.type <<- 17
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.type.T1 <<- 17
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.type.T2 <<- 17
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.type.T3 <<- 17
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.type.T4 <<- 17
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.type.T5 <<- 17
                }
                MGvar$zoomedplot.type <<- 17
                tclvalue(MGvar$Conf.PT.var) <<- "Filled Triangles"
                tclvalue(MGvar$Zoom.PT.var) <<- "Filled Triangles"
            }
            MGvar$activeplot.pointcol <<- MGvar$activeplot.pointcol.temp
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.pointcol.T1 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.pointcol.T2 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.pointcol.T3 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.pointcol.T4 <<- MGvar$activeplot.pointcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.pointcol.T5 <<- MGvar$activeplot.pointcol.temp
            }
            MGvar$zoomedplot.pointcol.temp <<- MGvar$activeplot.pointcol.temp
            MGvar$zoomedplot.pointcol <<- MGvar$zoomedplot.pointcol.temp
            tabplot()
        }
        tkplace(tkbutton(CPlotOps, text = "Change", width = 15, 
            command = function() ChangePoints()), relx = 0.325, 
            rely = 0.85, `in` = frameCPoints)
        ConfLines <- tk2frame(CPlotOpsNB)
        tkadd(CPlotOpsNB, ConfLines, text = "Lines")
        frameCLines <- tkwidget(ConfLines, "TitleFrame", text = "MDS Configuration Lines", 
            background = "white")
        tkplace(frameCLines, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ConfLines)
        tkplace(tklabel(frameCLines, text = "Display Regression Axes", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameCLines)
        Conf.RegLine.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.RegLine.CB, variable = MGvar$Conf.RegLine.val)
        tkplace(Conf.RegLine.CB, relx = 0.7, rely = 0.1, `in` = frameCLines)
        tkplace(tklabel(frameCLines, text = "Display Distance Lines Between \nPoints Corresponding to \nShepard Points", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameCLines)
        Conf.DistLine.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.DistLine.CB, variable = MGvar$Conf.DistLine.val)
        tkplace(Conf.DistLine.CB, relx = 0.7, rely = 0.45, `in` = frameCLines)
        ChangeColDist <- function() {
            MGvar$activeplot.distcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$activeplot.distcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$activeplot.distcol.temp) > 0) 
                tkconfigure(ConfColDist, bg = MGvar$activeplot.distcol.temp)
        }
        ConfColDist <- tkbutton(CPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$activeplot.distcol, command = function() ChangeColDist())
        ChangeLines <- function() {
            ShowR = as.character(tclvalue(MGvar$Conf.RegLine.val))
            if (ShowR == "1") {
                MGvar$activeplot.showreg <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showreg.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showreg.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showreg.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showreg.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showreg.T5 <<- "yes"
                }
                MGvar$zoomedplot.showreg <<- "yes"
                tclvalue(MGvar$Zoom.RegLine.val) <<- 1
                tclvalue(MGvar$Conf.RegLine.val) <<- 1
            }
            if (ShowR == "0") {
                MGvar$activeplot.showreg <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showreg.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showreg.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showreg.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showreg.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showreg.T5 <<- "no"
                }
                MGvar$zoomedplot.showreg <<- "no"
                tclvalue(MGvar$Zoom.RegLine.val) <<- 0
                tclvalue(MGvar$Conf.RegLine.val) <<- 0
            }
            MGvar$activeplot.regcol <<- MGvar$activeplot.regcol.temp
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.regcol.T1 <<- MGvar$activeplot.regcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.regcol.T2 <<- MGvar$activeplot.regcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.regcol.T3 <<- MGvar$activeplot.regcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.regcol.T4 <<- MGvar$activeplot.regcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.regcol.T5 <<- MGvar$activeplot.regcol.temp
            }
            MGvar$zoomedplot.regcol.temp <<- MGvar$activeplot.regcol.temp
            MGvar$zoomedplot.regcol <<- MGvar$zoomedplot.regcol.temp
            ShowD = as.character(tclvalue(MGvar$Conf.DistLine.val))
            if (ShowD == "1") {
                MGvar$activeplot.showdist <<- "yes"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showdist.T1 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showdist.T2 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showdist.T3 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showdist.T4 <<- "yes"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showdist.T5 <<- "yes"
                }
                tclvalue(MGvar$Conf.DistLine.val) <<- 1
            }
            if (ShowD == "0") {
                MGvar$activeplot.showdist <<- "no"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.showdist.T1 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.showdist.T2 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.showdist.T3 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.showdist.T4 <<- "no"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.showdist.T5 <<- "no"
                }
                tclvalue(MGvar$Conf.DistLine.val) <<- 0
            }
            MGvar$activeplot.distcol <<- MGvar$activeplot.distcol.temp
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.distcol.T1 <<- MGvar$activeplot.distcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.distcol.T2 <<- MGvar$activeplot.distcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.distcol.T3 <<- MGvar$activeplot.distcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.distcol.T4 <<- MGvar$activeplot.distcol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.distcol.T5 <<- MGvar$activeplot.distcol.temp
            }
            tabplot()
        }
        tkplace(tkbutton(CPlotOps, text = "Change", width = 15, 
            command = function() ChangeLines()), relx = 0.325, 
            rely = 0.85, `in` = frameCLines)
        ConfAxes <- tk2frame(CPlotOpsNB)
        tkadd(CPlotOpsNB, ConfAxes, text = "Axes")
        frameCAxes <- tkwidget(ConfAxes, "TitleFrame", text = "MDS Configuration Axes", 
            background = "white")
        tkplace(frameCAxes, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ConfAxes)
        tkplace(tklabel(frameCAxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameCAxes)
        Conf.AxesMeas.CB <- tk2checkbutton(CPlotOps)
        tkconfigure(Conf.AxesMeas.CB, variable = MGvar$Conf.AxesMeas.val)
        tkplace(Conf.AxesMeas.CB, relx = 0.7, rely = 0.1, `in` = frameCAxes)
        tkplace(tklabel(frameCAxes, text = "Axes Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameCAxes)
        ChangeColAxes <- function() {
            MGvar$activeplot.axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$activeplot.axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$activeplot.axescol.temp) > 0) 
                tkconfigure(ConfColAxes, bg = MGvar$activeplot.axescol.temp)
        }
        ConfColAxes <- tkbutton(CPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$activeplot.axescol, command = function() ChangeColAxes())
        tkplace(ConfColAxes, relx = 0.7, rely = 0.25, `in` = frameCAxes)
        ChangeAxes <- function() {
            ShowA = as.character(tclvalue(MGvar$Conf.AxesMeas.val))
            if (ShowA == "0") {
                MGvar$activeplot.yaxt <<- "n"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.yaxt.T1 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.yaxt.T2 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.yaxt.T3 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.yaxt.T4 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.yaxt.T5 <<- "n"
                }
                MGvar$activeplot.xaxt <<- "n"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.xaxt.T1 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.xaxt.T2 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.xaxt.T3 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.xaxt.T4 <<- "n"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.xaxt.T5 <<- "n"
                }
                MGvar$zoomedplot.yaxt <<- "n"
                MGvar$zoomedplot.xaxt <<- "n"
                tclvalue(MGvar$Zoom.AxesMeas.val) <<- 0
            }
            if (ShowA == "1") {
                MGvar$activeplot.yaxt <<- "s"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.yaxt.T1 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.yaxt.T2 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.yaxt.T3 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.yaxt.T4 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.yaxt.T5 <<- "s"
                }
                MGvar$activeplot.xaxt <<- "s"
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.xaxt.T1 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.xaxt.T2 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.xaxt.T3 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.xaxt.T4 <<- "s"
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.xaxt.T5 <<- "s"
                }
                MGvar$zoomedplot.yaxt <<- "s"
                MGvar$zoomedplot.xaxt <<- "s"
                tclvalue(MGvar$Conf.AxesMeas.val) <<- 1
                tclvalue(MGvar$Zoom.AxesMeas.val) <<- 1
            }
            MGvar$activeplot.axescol <<- MGvar$activeplot.axescol.temp
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.axescol.T1 <<- MGvar$activeplot.axescol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.axescol.T2 <<- MGvar$activeplot.axescol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.axescol.T3 <<- MGvar$activeplot.axescol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.axescol.T4 <<- MGvar$activeplot.axescol.temp
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.axescol.T5 <<- MGvar$activeplot.axescol.temp
            }
            MGvar$zoomedplot.axescol.temp <<- MGvar$activeplot.axescol
            MGvar$zoomedplot.axescol <<- MGvar$zoomedplot.axescol.temp
            tabplot()
        }
        tkplace(tkbutton(CPlotOps, text = "Change", width = 15, 
            command = function() ChangeAxes()), relx = 0.325, 
            rely = 0.85, `in` = frameCAxes)
        tkplace(CPlotOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = CPlotOps)
        OnOK.CP <- function() {
            ChangeGen()
            ChangePoints()
            ChangeLines()
            ChangeAxes()
            tkdestroy(CPlotOps)
        }
        OnCancel.CP <- function() {
            tkdestroy(CPlotOps)
        }
        OnDefault <- function() {
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$activeplot.title.show.T1 <<- "yes"
                MGvar$activeplot.distmeas.T1 <<- "yes"
                MGvar$activeplot.xlab.T1 <<- ""
                MGvar$activeplot.ylab.T1 <<- ""
                MGvar$activeplot.bg.T1 <<- "white"
                MGvar$activeplot.cex.T1 <<- 0.6
                MGvar$activeplot.labs.T1 <<- "yes"
                MGvar$activeplot.showpoints.T1 <<- "no"
                MGvar$activeplot.pointcol.T1 <<- "black"
                MGvar$activeplot.type.T1 <<- "1"
                MGvar$activeplot.xaxt.T1 <<- "n"
                MGvar$activeplot.yaxt.T1 <<- "n"
                MGvar$activeplot.axescol.T1 <<- "black"
                MGvar$activeplot.showreg.T1 <<- "no"
                MGvar$activeplot.regcol.T1 <<- "red"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$activeplot.title.show.T2 <<- "yes"
                MGvar$activeplot.distmeas.T2 <<- "yes"
                MGvar$activeplot.xlab.T2 <<- ""
                MGvar$activeplot.ylab.T2 <<- ""
                MGvar$activeplot.bg.T2 <<- "white"
                MGvar$activeplot.cex.T2 <<- 0.6
                MGvar$activeplot.labs.T2 <<- "yes"
                MGvar$activeplot.showpoints.T2 <<- "no"
                MGvar$activeplot.pointcol.T2 <<- "black"
                MGvar$activeplot.type.T2 <<- "1"
                MGvar$activeplot.xaxt.T2 <<- "n"
                MGvar$activeplot.yaxt.T2 <<- "n"
                MGvar$activeplot.axescol.T2 <<- "black"
                MGvar$activeplot.showreg.T2 <<- "no"
                MGvar$activeplot.regcol.T2 <<- "red"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$activeplot.title.show.T3 <<- "yes"
                MGvar$activeplot.distmeas.T3 <<- "yes"
                MGvar$activeplot.xlab.T3 <<- ""
                MGvar$activeplot.ylab.T3 <<- ""
                MGvar$activeplot.bg.T3 <<- "white"
                MGvar$activeplot.cex.T3 <<- 0.6
                MGvar$activeplot.labs.T3 <<- "yes"
                MGvar$activeplot.showpoints.T3 <<- "no"
                MGvar$activeplot.pointcol.T3 <<- "black"
                MGvar$activeplot.type.T3 <<- "1"
                MGvar$activeplot.xaxt.T3 <<- "n"
                MGvar$activeplot.yaxt.T3 <<- "n"
                MGvar$activeplot.axescol.T3 <<- "black"
                MGvar$activeplot.showreg.T3 <<- "no"
                MGvar$activeplot.regcol.T3 <<- "red"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$activeplot.title.show.T4 <<- "yes"
                MGvar$activeplot.distmeas.T4 <<- "yes"
                MGvar$activeplot.xlab.T4 <<- ""
                MGvar$activeplot.ylab.T4 <<- ""
                MGvar$activeplot.bg.T4 <<- "white"
                MGvar$activeplot.cex.T4 <<- 0.6
                MGvar$activeplot.labs.T4 <<- "yes"
                MGvar$activeplot.showpoints.T4 <<- "no"
                MGvar$activeplot.pointcol.T4 <<- "black"
                MGvar$activeplot.type.T4 <<- "1"
                MGvar$activeplot.xaxt.T4 <<- "n"
                MGvar$activeplot.yaxt.T4 <<- "n"
                MGvar$activeplot.axescol.T4 <<- "black"
                MGvar$activeplot.showreg.T4 <<- "no"
                MGvar$activeplot.regcol.T4 <<- "red"
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$activeplot.title.show.T5 <<- "yes"
                MGvar$activeplot.distmeas.T5 <<- "yes"
                MGvar$activeplot.xlab.T5 <<- ""
                MGvar$activeplot.ylab.T5 <<- ""
                MGvar$activeplot.bg.T5 <<- "white"
                MGvar$activeplot.cex.T5 <<- 0.6
                MGvar$activeplot.labs.T5 <<- "yes"
                MGvar$activeplot.showpoints.T5 <<- "no"
                MGvar$activeplot.pointcol.T5 <<- "black"
                MGvar$activeplot.type.T5 <<- "1"
                MGvar$activeplot.xaxt.T5 <<- "n"
                MGvar$activeplot.yaxt.T5 <<- "n"
                MGvar$activeplot.axescol.T5 <<- "black"
                MGvar$activeplot.showreg.T5 <<- "no"
                MGvar$activeplot.regcol.T5 <<- "red"
            }
            MGvar$activeplot.title.show <<- "yes"
            MGvar$activeplot.distmeas <<- "yes"
            MGvar$activeplot.xlab <<- ""
            MGvar$activeplot.ylab <<- ""
            MGvar$activeplot.bg <<- "white"
            MGvar$activeplot.bg.temp <<- "white"
            MGvar$activeplot.cex <<- 0.6
            MGvar$activeplot.labs <<- "yes"
            MGvar$activeplot.showpoints <<- "no"
            MGvar$activeplot.pointcol <<- "black"
            MGvar$activeplot.pointcol.temp <<- "black"
            MGvar$activeplot.type <<- "1"
            MGvar$activeplot.xaxt <<- "n"
            MGvar$activeplot.yaxt <<- "n"
            MGvar$activeplot.axescol <<- "black"
            MGvar$activeplot.axescol.temp <<- "black"
            MGvar$activeplot.showreg <<- "no"
            MGvar$activeplot.regcol <<- "red"
            PointColInitialise()
            tabplot()
            tclvalue(MGvar$Conf.Main.val) <<- 1
            tclvalue(MGvar$Conf.Dist.val) <<- 1
            tclvalue(MGvar$Conf.Ylab.val) <<- 0
            tclvalue(MGvar$Conf.Xlab.val) <<- 0
            tclvalue(MGvar$Conf.Points.val) <<- 0
            tclvalue(MGvar$Conf.Labels.val) <<- 1
            tclvalue(MGvar$Conf.PT.var) <<- "Empty Circles"
            tclvalue(MGvar$Conf.AxesMeas.val) <<- 0
            tclvalue(MGvar$Conf.RegLine.val) <<- 0
            MGvar$zoomedplot.title.show <<- "yes"
            MGvar$zoomedplot.distmeas <<- "yes"
            MGvar$zoomedplot.xlab <<- ""
            MGvar$zoomedplot.ylab <<- ""
            MGvar$zoomedplot.bg <<- "white"
            MGvar$zoomedplot.bg.temp <<- "white"
            MGvar$zoomedplot.cex <<- 1
            MGvar$zoomedplot.labs <<- "yes"
            MGvar$zoomedplot.showpoints <<- "no"
            MGvar$zoomedplot.pointcol <<- "black"
            MGvar$zoomedplot.pointcol.temp <<- "black"
            MGvar$zoomedplot.type <<- "1"
            MGvar$zoomedplot.xaxt <<- "n"
            MGvar$zoomedplot.yaxt <<- "n"
            MGvar$zoomedplot.axescol <<- "black"
            MGvar$zoomedplot.axescol.temp <<- "black"
            MGvar$zoomedplot.showreg <<- "no"
            MGvar$zoomedplot.regcol <<- "red"
            MGvar$zoomedplot.regcol.temp <<- "red"
            tclvalue(MGvar$Zoom.Main.val) <<- 1
            tclvalue(MGvar$Zoom.Dist.val) <<- 1
            tclvalue(MGvar$Zoom.Ylab.val) <<- 0
            tclvalue(MGvar$Zoom.Xlab.val) <<- 0
            tclvalue(MGvar$Zoom.Points.val) <<- 0
            tclvalue(MGvar$Zoom.Labels.val) <<- 1
            tclvalue(MGvar$Zoom.PT.var) <<- "Empty Circles"
            tclvalue(MGvar$Zoom.AxesMeas.val) <<- 0
            tclvalue(MGvar$Zoom.RegLine.val) <<- 0
            tkdestroy(CPlotOps)
        }
        tkplace(tkbutton(CPlotOps, text = "OK", width = 15, command = function() OnOK.CP()), 
            relx = 0.15, rely = 0.9, `in` = CPlotOps)
        tkplace(tkbutton(CPlotOps, text = "Default", width = 15, 
            command = function() OnDefault()), relx = 0.55, rely = 0.9, 
            `in` = CPlotOps)
        tkfocus(CPlotOps)
        tkwait.window(CPlotOps)
    }
    ShepPlotOptions <- function() {
        SPlotOps = tktoplevel()
        tkwm.resizable(SPlotOps, "0", "0")
        tkwm.deiconify(SPlotOps)
        tkwm.title(SPlotOps, "Shepard Plot Plotting Options")
        tkwm.geometry(SPlotOps, "350x400")
        SPlotcanvas = tkcanvas(SPlotOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(SPlotcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = SPlotOps)
        SPlotOpsNB <- tk2notebook(SPlotOps, tabs = NULL)
        ShepGen <- tkframe(SPlotOps)
        tkadd(SPlotOpsNB, ShepGen, text = "General")
        frameSGen <- tkwidget(ShepGen, "TitleFrame", text = "Shepard Plot General", 
            background = "white")
        tkplace(frameSGen, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheigh = 0.96, `in` = ShepGen)
        tkplace(tklabel(frameSGen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameSGen)
        Shep.Main.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Main.CB, variable = MGvar$Shep.Main.val)
        tkplace(Shep.Main.CB, relx = 0.7, rely = 0.1, `in` = frameSGen)
        tkplace(tklabel(frameSGen, text = "Display X & Y Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameSGen)
        Shep.Lab.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Lab.CB, variable = MGvar$Shep.Lab.val)
        tkplace(Shep.Lab.CB, relx = 0.7, rely = 0.25, `in` = frameSGen)
        tkplace(tklabel(frameSGen, text = "Display Plot Legend", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameSGen)
        Shep.Leg.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Leg.CB, variable = MGvar$Shep.Leg.val)
        tkplace(Shep.Leg.CB, relx = 0.7, rely = 0.4, `in` = frameSGen)
        tkplace(tklabel(frameSGen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameSGen)
        ChangeColShepBG <- function() {
            MGvar$shepplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$shepplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$shepplot.bg.temp) > 0) 
                tkconfigure(ShepColBG, bg = MGvar$shepplot.bg.temp)
        }
        ShepColBG <- tkbutton(SPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$shepplot.bg, command = function() ChangeColShepBG())
        tkplace(ShepColBG, relx = 0.7, rely = 0.55, `in` = frameSGen)
        ChangeSGen <- function() {
            MainS = as.character(tclvalue(MGvar$Shep.Main.val))
            if (MainS == "0") {
                MGvar$shepplot.title.show <<- "no"
                tclvalue(MGvar$Shep.Main.val) <<- 0
            }
            if (MainS == "1") {
                MGvar$shepplot.title.show <<- "yes"
                tclvalue(MGvar$Shep.Main.val) <<- 1
            }
            LabsS = as.character(tclvalue(MGvar$Shep.Lab.val))
            if (LabsS == "1") {
                MGvar$shepplot.labs.show <<- "yes"
                tclvalue(MGvar$Shep.Lab.val) <<- 1
            }
            if (LabsS == "0") {
                MGvar$shepplot.labs.show <<- "no"
                tclvalue(MGvar$Shep.Lab.val) <<- 0
            }
            LegS = as.character(tclvalue(MGvar$Shep.Leg.val))
            if (LegS == "1") {
                MGvar$shepplot.leg.show <<- "yes"
                tclvalue(MGvar$Shep.Leg.val) <<- 1
            }
            if (LegS == "0") {
                MGvar$shepplot.leg.show <<- "no"
                tclvalue(MGvar$Shep.Leg.val) <<- 0
            }
            MGvar$shepplot.bg <<- MGvar$shepplot.bg.temp
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEShep)
            }
        }
        tkplace(tkbutton(SPlotOps, text = "Change", width = 15, 
            command = function() ChangeSGen()), relx = 0.325, 
            rely = 0.85, `in` = frameSGen)
        ShepPoints <- tk2frame(SPlotOpsNB)
        tkadd(SPlotOpsNB, ShepPoints, text = "Points")
        frameSPoints <- tkwidget(ShepPoints, "TitleFrame", text = "Shepard Plot Points", 
            background = "white")
        tkplace(frameSPoints, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ShepPoints)
        tkplace(tklabel(frameSPoints, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameSPoints)
        Shep.Points.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Points.CB, variable = MGvar$Shep.Points.val)
        tkplace(Shep.Points.CB, relx = 0.7, rely = 0.1, `in` = frameSPoints)
        tkplace(tklabel(frameSPoints, text = "Display Point Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameSPoints)
        Shep.Labels.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Labels.CB, variable = MGvar$Shep.Labels.val)
        tkplace(Shep.Labels.CB, relx = 0.7, rely = 0.25, `in` = frameSPoints)
        tkplace(tklabel(frameSPoints, text = "Point Size", background = "white"), 
            relx = 0.1, rely = 0.4, `in` = frameSPoints)
        Shep.PS.spin <- tk2spinbox(SPlotOps, from = 0.1, to = 2, 
            increment = 0.1, width = 12)
        tkconfigure(Shep.PS.spin, textvariable = MGvar$Shep.PS.var)
        tkplace(Shep.PS.spin, relx = 0.6, rely = 0.4, `in` = frameSPoints)
        tkplace(tklabel(frameSPoints, text = "Point Type", background = "white"), 
            relx = 0.1, rely = 0.55, `in` = frameSPoints)
        shapes <- c("Empty Boxes", "Filled Boxes", "Crosses", 
            "Empty Triangles", "Filled Triangles", "Filled Circles", 
            "Empty Circles")
        Shep.PT.ComboBox <- tkwidget(SPlotOps, "ComboBox", editable = FALSE, 
            values = shapes, width = 12, textvariable = MGvar$Shep.PT.var)
        tkplace(Shep.PT.ComboBox, relx = 0.6, rely = 0.55, `in` = frameSPoints)
        tkplace(tklabel(frameSPoints, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameSPoints)
        ChangeColShepPoints <- function() {
            MGvar$shepplot.pointcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$shepplot.pointcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$shepplot.pointcol.temp) > 0) 
                tkconfigure(ShepColPoints, bg = MGvar$shepplot.pointcol.temp)
        }
        ShepColPoints <- tkbutton(SPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$shepplot.pointcol, command = function() ChangeColShepPoints())
        tkplace(ShepColPoints, relx = 0.7, rely = 0.7, `in` = frameSPoints)
        ChangeShepPoints <- function() {
            SDispP = as.character(tclvalue(MGvar$Shep.Points.val))
            if (SDispP == "0") {
                MGvar$shepplot.showpoints <<- "no"
                tclvalue(MGvar$Shep.Points.val) <<- 0
            }
            if (SDispP == "1") {
                MGvar$shepplot.showpoints <<- "yes"
                tclvalue(MGvar$Shep.Points.val) <<- 1
            }
            SDispL = as.character(tclvalue(MGvar$Shep.Labels.val))
            if (SDispL == "0") {
                MGvar$shepplot.showlabels <<- "no"
                tclvalue(MGvar$Shep.Labels.val) <<- 0
            }
            if (SDispL == "1") {
                MGvar$shepplot.showlabels <<- "yes"
                tclvalue(MGvar$Shep.Labels.val) <<- 1
            }
            SSizeP = as.numeric(tclvalue(MGvar$Shep.PS.var))
            MGvar$shepplot.cex <<- SSizeP
            STypeP = as.character(tclvalue(MGvar$Shep.PT.var))
            if (STypeP == "Empty Circles") {
                MGvar$shepplot.type <<- 1
                tclvalue(MGvar$Shep.PT.var) <<- "Empty Circles"
            }
            if (STypeP == "Filled Circles") {
                MGvar$shepplot.type <<- 16
                tclvalue(MGvar$Shep.PT.var) <<- "Filled Circles"
            }
            if (STypeP == "Empty Boxes") {
                MGvar$shepplot.type <<- 22
                tclvalue(MGvar$Shep.PT.var) <<- "Empty Boxes"
            }
            if (STypeP == "Filled Boxes") {
                MGvar$shepplot.type <<- 15
                tclvalue(MGvar$Shep.PT.var) <<- "Filled Boxes"
            }
            if (STypeP == "Crosses") {
                MGvar$shepplot.type <<- 4
                tclvalue(MGvar$Shep.PT.var) <<- "Crosses"
            }
            if (STypeP == "Empty Triangles") {
                MGvar$shepplot.type <<- 24
                tclvalue(MGvar$Shep.PT.var) <<- "Empty Triangles"
            }
            if (STypeP == "Filled Triangles") {
                MGvar$shepplot.type <<- 17
                tclvalue(MGvar$Shep.PT.var) <<- "Filled Triangles"
            }
            MGvar$shepplot.pointcol <<- MGvar$shepplot.pointcol.temp
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEShep)
            }
        }
        tkplace(tkbutton(SPlotOps, text = "Change", width = 15, 
            command = function() ChangeShepPoints()), relx = 0.325, 
            rely = 0.85, `in` = frameSPoints)
        ShepLines <- tk2frame(SPlotOpsNB)
        tkadd(SPlotOpsNB, ShepLines, text = "Lines")
        frameSLines <- tkwidget(ShepLines, "TitleFrame", text = "Shepard Plot Lines", 
            background = "white")
        tkplace(frameSLines, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ShepLines)
        tkplace(tklabel(frameSLines, text = "Display Shepard Line", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameSLines)
        Shep.Curve.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Curve.CB, variable = MGvar$Shep.Curve.val)
        tkplace(Shep.Curve.CB, relx = 0.7, rely = 0.1, `in` = frameSLines)
        tkplace(tklabel(frameSLines, text = "Line Type", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameSLines)
        Shep.LineT.ComboBox <- tkwidget(SPlotOps, "ComboBox", 
            editable = FALSE, values = c("Solid Line", "Dashed", 
                "Dotted"), width = 12, textvariable = MGvar$Shep.LineT.val)
        tkplace(Shep.LineT.ComboBox, relx = 0.6, rely = 0.25, 
            `in` = frameSLines)
        tkplace(tklabel(frameSLines, text = "Line Colour", background = "white"), 
            relx = 0.1, rely = 0.4, `in` = frameSLines)
        ChangeColShepLine <- function() {
            MGvar$shepplot.curvecol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$shepplot.curvecol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$shepplot.curvecol.temp) > 0) 
                tkconfigure(ShepColCurve, bg = MGvar$shepplot.curvecol.temp)
        }
        ShepColCurve <- tkbutton(SPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$shepplot.curvecol, command = function() ChangeColShepLine())
        tkplace(ShepColCurve, relx = 0.7, rely = 0.4, `in` = frameSLines)
        tkplace(tklabel(frameSLines, text = "Bold", background = "white"), 
            relx = 0.1, rely = 0.55, `in` = frameSLines)
        Shep.Bold.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Bold.CB, variable = MGvar$Shep.Bold.val)
        tkplace(Shep.Bold.CB, relx = 0.7, rely = 0.55, `in` = frameSLines)
        ChangeSLines <- function() {
            ShCurve = as.character(tclvalue(MGvar$Shep.Curve.val))
            if (ShCurve == "1") {
                MGvar$shepplot.curve.show <<- "yes"
                tclvalue(MGvar$Shep.Curve.val) <<- 1
            }
            if (ShCurve == "0") {
                MGvar$shepplot.curve.show <<- "no"
                tclvalue(MGvar$Shep.Curve.val) <<- 0
            }
            STypeL = as.character(tclvalue(MGvar$Shep.LineT.val))
            if (STypeL == "Solid Line") {
                MGvar$shepplot.curve.type <<- 1
                tclvalue(MGvar$Shep.LineT.val) <<- "Solid Line"
            }
            if (STypeL == "Dashed") {
                MGvar$shepplot.curve.type <<- 2
                tclvalue(MGvar$Shep.LineT.val) <<- "Dashed"
            }
            if (STypeL == "Dotted") {
                MGvar$shepplot.curve.type <<- 3
                tclvalue(MGvar$Shep.LineT.val) <<- "Dotted"
            }
            ShBold = as.character(tclvalue(MGvar$Shep.Bold.val))
            if (ShBold == "1") {
                MGvar$shepplot.bold <<- TRUE
                tclvalue(MGvar$Shep.Bold.val) <<- 1
            }
            if (ShBold == "0") {
                MGvar$shepplot.bold <<- FALSE
                tclvalue(MGvar$Shep.Bold.val) <<- 0
            }
            MGvar$shepplot.curvecol <<- MGvar$shepplot.curvecol.temp
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEShep)
            }
        }
        tkplace(tkbutton(SPlotOps, text = "Change", width = 15, 
            command = function() ChangeSLines()), relx = 0.325, 
            rely = 0.85, `in` = frameSLines)
        ShepAxes <- tk2frame(SPlotOpsNB)
        tkadd(SPlotOpsNB, ShepAxes, text = "Axes")
        frameSAxes <- tkwidget(ShepAxes, "TitleFrame", text = "Shepard Plot Axes", 
            background = "white")
        tkplace(frameSAxes, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ShepAxes)
        tkplace(tklabel(frameSAxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameSAxes)
        Shep.Axes.CB <- tk2checkbutton(SPlotOps)
        tkconfigure(Shep.Axes.CB, variable = MGvar$Shep.Axes.val)
        tkplace(Shep.Axes.CB, relx = 0.7, rely = 0.1, `in` = frameSAxes)
        tkplace(tklabel(frameSAxes, text = "Axes Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameSAxes)
        ChangeShepAxescol <- function() {
            MGvar$shepplot.Axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$shepplot.Axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$shepplot.Axescol.temp) > 0) 
                tkconfigure(ShepColAxes, bg = MGvar$shepplot.Axescol.temp)
        }
        ShepColAxes <- tkbutton(SPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$shepplot.Axescol, command = function() ChangeShepAxescol())
        tkplace(ShepColAxes, relx = 0.7, rely = 0.25, `in` = frameSAxes)
        ChangeSAxes <- function() {
            ShAxes = as.character(tclvalue(MGvar$Shep.Axes.val))
            if (ShAxes == "1") {
                MGvar$shepplot.Axes.xaxt <<- "s"
                MGvar$shepplot.Axes.yaxt <<- "s"
                tclvalue(MGvar$Shep.Axes.val) <<- 1
            }
            if (ShAxes == "0") {
                MGvar$shepplot.Axes.xaxt <<- "n"
                MGvar$shepplot.Axes.yaxt <<- "n"
                tclvalue(MGvar$Shep.Axes.val) <<- 0
            }
            MGvar$shepplot.Axescol <<- MGvar$shepplot.Axescol.temp
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEShep)
            }
        }
        tkplace(tkbutton(SPlotOps, text = "Change", width = 15, 
            command = function() ChangeSAxes()), relx = 0.325, 
            rely = 0.85, `in` = frameSAxes)
        tkplace(SPlotOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = SPlotOps)
        OnOK.SP <- function() {
            ChangeSGen()
            ChangeShepPoints()
            ChangeSLines()
            ChangeSAxes()
            tkdestroy(SPlotOps)
        }
        OnCancel.SP <- function() {
            tkdestroy(SPlotOps)
        }
        OnDefault <- function() {
            MGvar$shepplot.title.show <<- "yes"
            MGvar$shepplot.labs.show <<- "yes"
            MGvar$shepplot.leg.show <<- "no"
            MGvar$shepplot.bg <<- "white"
            MGvar$shepplot.bg.temp <<- "white"
            MGvar$shepplot.showpoints <<- "yes"
            MGvar$shepplot.showlabels <<- "no"
            MGvar$shepplot.cex <<- 0.6
            MGvar$shepplot.type <<- 1
            MGvar$shepplot.pointcol.temp <<- "black"
            MGvar$shepplot.pointcol <<- "black"
            MGvar$shepplot.curve.show <<- "yes"
            MGvar$shepplot.curvecol.temp <<- "red"
            MGvar$shepplot.curvecol <<- "red"
            MGvar$shepplot.Axes.xaxt <<- "s"
            MGvar$shepplot.Axes.yaxt <<- "s"
            MGvar$shepplot.Axescol.temp <<- "black"
            MGvar$shepplot.Axescol <<- "black"
            MGvar$Shep.Main.val <<- tclVar("1")
            MGvar$Shep.Lab.val <<- tclVar("1")
            MGvar$Shep.Leg.val <<- tclVar("0")
            MGvar$Shep.Points.val <<- tclVar("1")
            MGvar$Shep.PS.var <<- tclVar(MGvar$shepplot.cex)
            MGvar$Shep.PT.var <<- tclVar("Empty Circles")
            MGvar$Shep.Curve.val <<- tclVar("1")
            MGvar$Shep.LineT.val <<- tclVar("Solid Line")
            MGvar$Shep.Axes.val <<- tclVar("1")
            MGvar$Shep.Labels.val <<- tclVar("0")
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEShep)
            }
            tkdestroy(SPlotOps)
        }
        tkplace(tkbutton(SPlotOps, text = "OK", width = 15, command = function() OnOK.SP()), 
            relx = 0.15, rely = 0.9, `in` = SPlotOps)
        tkplace(tkbutton(SPlotOps, text = "Default", width = 15, 
            command = function() OnDefault()), relx = 0.55, rely = 0.9, 
            `in` = SPlotOps)
        tkfocus(SPlotOps)
        tkwait.window(SPlotOps)
    }
    EnlargedShep <- function() {
        ReusePOW()
        MGvar$EnShep.switch <<- "on"
        MGcomp$EShep <<- tktoplevel()
        tkwm.geometry(MGcomp$EShep, "650x575")
        tkwm.title(MGcomp$EShep, "Enlarged Shepard Plot")
        EShepcanvas = tkcanvas(MGcomp$EShep, width = 773, height = 575, 
            bg = col.sec)
        tkplace(EShepcanvas, `in` = MGcomp$EShep)
        EShscale = 1.6
        ESvscale = 1.3
        MGcomp$imgEShep <<- tkrplot(MGcomp$EShep, function() plotShepard(MGvar$distmat, 
            MGvar$MDSmat), hscale = EShscale, vscale = ESvscale)
        tkplace(MGcomp$imgEShep, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            `in` = MGcomp$EShep)
        CopyShepToClip <- function() {
            tkrreplot(MGcomp$imgEShep)
        }
        tkbind(MGcomp$imgEShep, "<Button-1>", POPShep.OnLeftClick)
        tkconfigure(MGcomp$imgEShep, cursor = "hand2")
        tkplace(tkbutton(MGcomp$EShep, text = "Copy to Clipboard", 
            width = 20, command = function() CopyShepToClip()), 
            relx = 0.22, rely = 0.93, `in` = MGcomp$EShep)
        tkplace(tkbutton(MGcomp$EShep, text = "Plot Options", 
            width = 20, command = function() ShepPlotOptions()), 
            relx = 0.55, rely = 0.93, `in` = MGcomp$EShep)
        tkbind(MGcomp$imgEShep, "<Destroy>", function() {
            MGvar$EnShep.switch <<- "off"
        })
        EShepPopupMenu <- tkmenu(MGcomp$imgEShep, tearoff = FALSE)
        tkadd(EShepPopupMenu, "command", label = "Clear Added Labels", 
            command = ClearShepLabels)
        tkadd(EShepPopupMenu, "command", label = "Label Specific Point", 
            command = Shep.LabelSpecificPoint)
        RightClickPOPShep <- function(x, y) {
            rootx <- as.integer(tkwinfo("rootx", MGcomp$imgEShep))
            rooty <- as.integer(tkwinfo("rooty", MGcomp$imgEShep))
            xTxt <- as.integer(x) + rootx
            yTxt <- as.integer(y) + rooty
            tcl("tk_popup", EShepPopupMenu, xTxt, yTxt)
        }
        tkbind(MGcomp$imgEShep, "<Button-3>", RightClickPOPShep)
        tkbind(MGcomp$EShep, "<Configure>", resize.EShep)
    }
    resize.EShep <- function() {
        ESh.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EShep)))
        ESh.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EShep)))
        WidthScale = 650/1.6
        EShscale <- ESh.width/WidthScale
        dimrat = 1.6/1.3
        ESvscale <- EShscale/dimrat
        tkrreplot(MGcomp$imgEShep, hscale = EShscale, vscale = ESvscale)
    }
    Replot.imgEShep <- function() {
        resize.EShep()
    }
    ClearShepLabels <- function() {
        MGvar$Shep.indexLabeled <<- c()
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
        if (MGvar$EnShep.switch == "on") {
            tkrreplot(MGcomp$imgEShep)
        }
        tabplot()
    }
    POPShep.OnLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", MGcomp$imgEShep)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", MGcomp$imgEShep)))
        xMin <- MGvar$Shep.parPlotSize[1] * width
        xMax <- MGvar$Shep.parPlotSize[2] * width
        yMin <- MGvar$Shep.parPlotSize[3] * height
        yMax <- MGvar$Shep.parPlotSize[4] * height
        rangeX <- MGvar$Shep.usrCoords[2] - MGvar$Shep.usrCoords[1]
        rangeY <- MGvar$Shep.usrCoords[4] - MGvar$Shep.usrCoords[3]
        imgXcoords <- (MGvar$Shepx - MGvar$Shep.usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$Shepy - MGvar$Shep.usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$Shep.usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$Shep.usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        Shep.labelClosestPoint(xClick, yClick, imgXcoords, imgYcoords)
    }
    StressPlotOps = function() {
        StPlotOps = tktoplevel()
        tkwm.resizable(StPlotOps, "0", "0")
        tkwm.deiconify(StPlotOps)
        tkwm.title(StPlotOps, "Stress Plot Options")
        tkwm.geometry(StPlotOps, "350x400")
        StPlotcanvas = tkcanvas(StPlotOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(StPlotcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = StPlotOps)
        StPlotOpsNB <- tk2notebook(StPlotOps, tabs = NULL)
        StressGen <- tk2frame(StPlotOpsNB)
        tkadd(StPlotOpsNB, StressGen, text = "General")
        frameStGen <- tkwidget(StressGen, "TitleFrame", text = "Stress Plot General", 
            background = "white")
        tkplace(frameStGen, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = StressGen)
        tkplace(tklabel(frameStGen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameStGen)
        Stress.Main.CB <- tk2checkbutton(StPlotOps)
        tkconfigure(Stress.Main.CB, variable = MGvar$Stress.Main.val)
        tkplace(Stress.Main.CB, relx = 0.7, rely = 0.1, `in` = frameStGen)
        tkplace(tklabel(frameStGen, text = "Display Process Time", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameStGen)
        Stress.Time.CB <- tk2checkbutton(StPlotOps)
        tkconfigure(Stress.Time.CB, variable = MGvar$Stress.Time.val)
        tkplace(Stress.Time.CB, relx = 0.7, rely = 0.25, `in` = frameStGen)
        tkplace(tklabel(frameStGen, text = "Display X & Y Labels", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameStGen)
        Stress.Lab.CB <- tk2checkbutton(StPlotOps)
        tkconfigure(Stress.Lab.CB, variable = MGvar$Stress.Lab.val)
        tkplace(Stress.Lab.CB, relx = 0.7, rely = 0.4, `in` = frameStGen)
        tkplace(tklabel(frameStGen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameStGen)
        ChangeColStressBG <- function() {
            MGvar$stressplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$stressplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$stressplot.bg.temp) > 0) 
                tkconfigure(StressColBG, bg = MGvar$stressplot.bg.temp)
        }
        StressColBG <- tkbutton(StPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$stressplot.bg, command = function() ChangeColStressBG())
        tkplace(StressColBG, relx = 0.7, rely = 0.55, `in` = frameStGen)
        ChangeStGen <- function() {
            MainSt = as.character(tclvalue(MGvar$Stress.Main.val))
            if (MainSt == "1") {
                MGvar$stressplot.title.show <<- "yes"
                tclvalue(MGvar$Stress.Main.val) <<- 1
            }
            if (MainSt == "0") {
                MGvar$stressplot.title.show <<- "no"
                tclvalue(MGvar$Stress.Main.val) <<- 0
            }
            DTime = as.character(tclvalue(MGvar$Stress.Time.val))
            if (DTime == "1") {
                MGvar$stressplot.time.show <<- "yes"
                tclvalue(MGvar$Stress.Time.val) <<- 1
            }
            if (DTime == "0") {
                MGvar$stressplot.time.show <<- "no"
                tclvalue(MGvar$Stress.Time.val) <<- 0
            }
            LabsSt = as.character(tclvalue(MGvar$Stress.Lab.val))
            if (LabsSt == "1") {
                MGvar$stressplot.labs.show <<- "yes"
                tclvalue(MGvar$Stress.Lab.val) <<- 1
            }
            if (LabsSt == "0") {
                MGvar$stressplot.labs.show <<- "no"
                tclvalue(MGvar$Stress.Lab.val) <<- 0
            }
            MGvar$stressplot.bg <<- MGvar$stressplot.bg.temp
            tkrreplot(imgstress)
            if (MGvar$EnStress.switch == "on") {
                tkrreplot(MGcomp$imgEStress)
            }
            tkrreplot(imgstress2)
            if (MGvar$EnStress2.switch == "on") {
                tkrreplot(MGcomp$imgEStress2)
            }
        }
        tkplace(tkbutton(StPlotOps, text = "Change", width = 15, 
            command = function() ChangeStGen()), relx = 0.325, 
            rely = 0.85, `in` = frameStGen)
        StressLines <- tk2frame(StPlotOpsNB)
        tkadd(StPlotOpsNB, StressLines, text = "Lines")
        frameStLines <- tkwidget(StressLines, "TitleFrame", text = "Stress Plot Lines", 
            background = "white")
        tkplace(frameStLines, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = StressLines)
        tkplace(tklabel(frameStLines, text = "Line Type", background = "white"), 
            relx = 0.1, rely = 0.1, `in` = frameStLines)
        Stress.LineT.ComboBox <- tkwidget(StPlotOps, "ComboBox", 
            editable = FALSE, values = c("Solid Line", "Dashed", 
                "Dotted"), width = 12, textvariable = MGvar$Stress.LineT.val)
        tkplace(Stress.LineT.ComboBox, relx = 0.6, rely = 0.1, 
            `in` = frameStLines)
        tkplace(tklabel(frameStLines, text = "Line Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameStLines)
        ChangeColStressLine <- function() {
            MGvar$stressplot.curvecol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$stressplot.curvecol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$stressplot.curvecol.temp) > 0) 
                tkconfigure(StressColCurve, bg = MGvar$stressplot.curvecol.temp)
        }
        StressColCurve <- tkbutton(StPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$stressplot.curvecol, command = function() ChangeColStressLine())
        tkplace(StressColCurve, relx = 0.7, rely = 0.25, `in` = frameStLines)
        ChangeStLines <- function() {
            StTypeL = as.character(tclvalue(MGvar$Stress.LineT.val))
            if (StTypeL == "Solid Line") {
                MGvar$stressplot.curve.type <<- 1
                tclvalue(MGvar$Stress.LineT.val) <<- "Solid Line"
            }
            if (StTypeL == "Dashed") {
                MGvar$stressplot.curve.type <<- 2
                tclvalue(MGvar$Stress.LineT.val) <<- "Dashed"
            }
            if (StTypeL == "Dotted") {
                MGvar$stressplot.curve.type <<- 3
                tclvalue(MGvar$Stress.LineT.val) <<- "Dotted"
            }
            MGvar$stressplot.curvecol <<- MGvar$stressplot.curvecol.temp
            tkrreplot(imgstress)
            if (MGvar$EnStress.switch == "on") {
                tkrreplot(MGcomp$imgEStress)
            }
            tkrreplot(imgstress2)
            if (MGvar$EnStress2.switch == "on") {
                tkrreplot(MGcomp$imgEStress2)
            }
        }
        tkplace(tkbutton(StPlotOps, text = "Change", width = 15, 
            command = function() ChangeStLines()), relx = 0.325, 
            rely = 0.85, `in` = frameStLines)
        StressAxes <- tk2frame(StPlotOpsNB)
        tkadd(StPlotOpsNB, StressAxes, text = "Axes")
        frameStAxes <- tkwidget(StressAxes, "TitleFrame", text = "Stress Plot Axes", 
            background = "white")
        tkplace(frameStAxes, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = StressAxes)
        tkplace(tklabel(frameStAxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameStAxes)
        Stress.Axes.CB <- tk2checkbutton(StPlotOps)
        tkconfigure(Stress.Axes.CB, variable = MGvar$Stress.Axes.val)
        tkplace(Stress.Axes.CB, relx = 0.7, rely = 0.1, `in` = frameStAxes)
        tkplace(tklabel(frameStAxes, text = "Axes Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameStAxes)
        ChangeStressAxescol <- function() {
            MGvar$stressplot.Axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$stressplot.Axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$stressplot.Axescol.temp) > 0) 
                tkconfigure(StressColAxes, bg = MGvar$stressplot.Axescol.temp)
        }
        StressColAxes <- tkbutton(StPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$stressplot.Axescol, command = function() ChangeStressAxescol())
        tkplace(StressColAxes, relx = 0.7, rely = 0.25, `in` = frameStAxes)
        ChangeStAxes <- function() {
            StAxes = as.character(tclvalue(MGvar$Stress.Axes.val))
            if (StAxes == "1") {
                MGvar$stressplot.Axes.xaxt <<- "s"
                MGvar$stressplot.Axes.yaxt <<- "s"
                tclvalue(MGvar$Stress.Axes.val) <<- 1
            }
            if (StAxes == "0") {
                MGvar$stressplot.Axes.xaxt <<- "n"
                MGvar$stressplot.Axes.yaxt <<- "n"
                tclvalue(MGvar$Stress.Axes.val) <<- 0
            }
            MGvar$stressplot.Axescol <<- MGvar$stressplot.Axescol.temp
            tkrreplot(imgstress)
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEStress)
            }
            tkrreplot(imgstress2)
            if (MGvar$EnStress2.switch == "on") {
                tkrreplot(MGcomp$imgEStress2)
            }
        }
        tkplace(tkbutton(StPlotOps, text = "Change", width = 15, 
            command = function() ChangeStAxes()), relx = 0.325, 
            rely = 0.85, `in` = frameStAxes)
        tkplace(StPlotOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = StPlotOps)
        OnOK.St <- function() {
            ChangeStGen()
            ChangeStLines()
            tkdestroy(StPlotOps)
        }
        OnDefault <- function() {
            MGvar$Stress.Main.val <<- tclVar("1")
            MGvar$Stress.Lab.val <<- tclVar("1")
            MGvar$Stress.Time.val <<- tclVar("1")
            MGvar$Stress.LineT.val <<- tclVar("Solid Line")
            MGvar$Stress.Axes.val <<- tclVar("1")
            MGvar$stressplot.title.show <<- "yes"
            MGvar$stressplot.time.show <<- "yes"
            MGvar$stressplot.labs.show <<- "yes"
            MGvar$stressplot.bg <<- "white"
            MGvar$stressplot.bg.temp <<- "white"
            MGvar$stressplot.curve.type <<- 1
            MGvar$stressplot.curvecol <<- "black"
            MGvar$stressplot.curvecol.temp <<- "black"
            MGvar$stressplot.Axes.xaxt <<- "s"
            MGvar$stressplot.Axes.yaxt <<- "s"
            MGvar$stressplot.Axescol <<- "black"
            MGvar$stressplot.Axescol.temp <<- "black"
            tkrreplot(imgstress)
            if (MGvar$EnShep.switch == "on") {
                tkrreplot(MGcomp$imgEStress)
            }
            tkdestroy(StPlotOps)
        }
        tkplace(tkbutton(StPlotOps, text = "OK", width = 15, 
            command = function() OnOK.St()), relx = 0.15, rely = 0.9, 
            `in` = StPlotOps)
        tkplace(tkbutton(StPlotOps, text = "Default", width = 15, 
            command = function() OnDefault()), relx = 0.55, rely = 0.9, 
            `in` = StPlotOps)
        tkfocus(StPlotOps)
        tkwait.window(StPlotOps)
    }
    ScreePlotOptions <- function() {
        ScPlotOps = tktoplevel()
        tkwm.resizable(ScPlotOps, "0", "0")
        tkwm.deiconify(ScPlotOps)
        tkwm.title(ScPlotOps, "Scree Plot Plotting Options")
        tkwm.geometry(ScPlotOps, "350x400")
        ScPlotcanvas = tkcanvas(ScPlotOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(ScPlotcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = ScPlotOps)
        ScPlotOpsNB <- tk2notebook(ScPlotOps, tabs = NULL)
        ScreeGen <- tk2frame(ScPlotOps)
        tkadd(ScPlotOpsNB, ScreeGen, text = "General")
        frameScGen <- tkwidget(ScreeGen, "TitleFrame", text = "Scree General", 
            background = "white")
        tkplace(frameScGen, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheigh = 0.96, `in` = ScreeGen)
        tkplace(tklabel(frameScGen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameScGen)
        Scree.Main.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.Main.CB, variable = MGvar$Scree.Main.val)
        tkplace(Scree.Main.CB, relx = 0.7, rely = 0.1, `in` = frameScGen)
        tkplace(tklabel(frameScGen, text = "Display X & Y Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameScGen)
        Scree.Lab.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.Lab.CB, variable = MGvar$Scree.Lab.val)
        tkplace(Scree.Lab.CB, relx = 0.7, rely = 0.25, `in` = frameScGen)
        tkplace(tklabel(frameScGen, text = "Display Plot Legend", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameScGen)
        Scree.Leg.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.Leg.CB, variable = MGvar$Scree.Leg.val)
        tkplace(Scree.Leg.CB, relx = 0.7, rely = 0.4, `in` = frameScGen)
        tkplace(tklabel(frameScGen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameScGen)
        ChangeColBG <- function() {
            MGvar$screeplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$screeplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$screeplot.bg.temp) > 0) 
                tkconfigure(ScreeColBG, bg = MGvar$screeplot.bg.temp)
        }
        ScreeColBG <- tkbutton(ScPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$screeplot.bg, command = function() ChangeColBG())
        tkplace(ScreeColBG, relx = 0.7, rely = 0.55, `in` = frameScGen)
        ChangeScGen <- function() {
            MainSc = as.character(tclvalue(MGvar$Scree.Main.val))
            if (MainSc == "0") {
                MGvar$screeplot.title.show <<- "no"
                tclvalue(MGvar$Scree.Main.val) <<- 0
            }
            if (MainSc == "1") {
                MGvar$screeplot.title.show <<- "yes"
            }
            LabsSc = as.character(tclvalue(MGvar$Scree.Lab.val))
            if (LabsSc == "0") {
                MGvar$screeplot.labs.show <<- "no"
                tclvalue(MGvar$Scree.Lab.val) <<- 0
            }
            if (LabsSc == "1") {
                MGvar$screeplot.labs.show <<- "yes"
            }
            LegSc = as.character(tclvalue(MGvar$Scree.Leg.val))
            if (LegSc == "0") {
                MGvar$screeplot.leg.show <<- "no"
                tclvalue(MGvar$Scree.Leg.val) <<- 0
            }
            if (LegSc == "1") {
                MGvar$screeplot.leg.show <<- "yes"
            }
            MGvar$screeplot.bg <<- MGvar$screeplot.bg.temp
            tkrreplot(imgscree)
            if (EnScree.switch == "on") {
                tkrreplot(MGcomp$imgES)
            }
        }
        tkplace(tkbutton(ScPlotOps, text = "Change", width = 15, 
            command = function() ChangeScGen()), relx = 0.325, 
            rely = 0.85, `in` = frameScGen)
        ScreePoints <- tk2frame(ScPlotOpsNB)
        tkadd(ScPlotOpsNB, ScreePoints, text = "Points")
        frameScPoints <- tkwidget(ScreePoints, "TitleFrame", 
            text = "Scree Plot Points", background = "white")
        tkplace(frameScPoints, relx = 0.02, relwidth = 0.96, 
            rely = 0.02, relheight = 0.96, `in` = ScreePoints)
        tkplace(tklabel(frameScPoints, text = "Display All Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameScPoints)
        Scree.Points.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.Points.CB, variable = MGvar$Scree.Points.val)
        tkplace(Scree.Points.CB, relx = 0.8, rely = 0.1, `in` = frameScPoints)
        tkplace(tklabel(frameScPoints, text = "Highlight Current Dimension Point", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameScPoints)
        Scree.CDim.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.CDim.CB, variable = MGvar$Scree.CDim.val)
        tkplace(Scree.CDim.CB, relx = 0.8, rely = 0.25, `in` = frameScPoints)
        tkplace(tklabel(frameScPoints, text = "Highlight Optimum Dimension Point", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameScPoints)
        Scree.ODim.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.ODim.CB, variable = MGvar$Scree.ODim.val)
        tkplace(Scree.ODim.CB, relx = 0.8, rely = 0.4, `in` = frameScPoints)
        tkplace(tklabel(frameScPoints, text = "Current Dimension Point Colour", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameScPoints)
        ChangeColCDim <- function() {
            MGvar$screeplot.Ccol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$screeplot.Ccol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$screeplot.Ccol.temp) > 0) 
                tkconfigure(ScreeCcol, bg = MGvar$screeplot.Ccol.temp)
        }
        ScreeCcol <- tkbutton(ScPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$screeplot.Ccol, command = function() ChangeColCDim())
        tkplace(ScreeCcol, relx = 0.8, rely = 0.55, `in` = frameScPoints)
        tkplace(tklabel(frameScPoints, text = "Optimum Dimension Point Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameScPoints)
        ChangeColODim <- function() {
            MGvar$screeplot.Ocol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$screeplot.Ocol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$screeplot.Ocol.temp) > 0) 
                tkconfigure(ScreeOcol, bg = MGvar$screeplot.Ocol.temp)
        }
        ScreeOcol <- tkbutton(ScPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$screeplot.Ocol, command = function() ChangeColODim())
        tkplace(ScreeOcol, relx = 0.8, rely = 0.7, `in` = frameScPoints)
        ChangeScPoints <- function() {
            PointsSc = as.character(tclvalue(MGvar$Scree.Points.val))
            if (PointsSc == "0") {
                MGvar$screeplot.points.show <<- "no"
            }
            if (PointsSc == "1") {
                MGvar$screeplot.points.show <<- "yes"
                tclvalue(MGvar$Scree.Points.val) <<- 1
            }
            CdimSc = as.character(tclvalue(MGvar$Scree.CDim.val))
            if (CdimSc == "1") {
                MGvar$screeplot.Cdim.show <<- "yes"
            }
            if (CdimSc == "0") {
                MGvar$screeplot.Cdim.show <<- "no"
                tclvalue(MGvar$Scree.CDim.val) <<- 0
            }
            OdimSc = as.character(tclvalue(MGvar$Scree.ODim.val))
            if (OdimSc == "1") {
                MGvar$screeplot.Odim.show <<- "yes"
            }
            if (OdimSc == "0") {
                MGvar$screeplot.Odim.show <<- "no"
                tclvalue(MGvar$Scree.ODim.val) <<- 0
            }
            MGvar$screeplot.Ocol <<- MGvar$screeplot.Ocol.temp
            MGvar$screeplot.Ccol <<- MGvar$screeplot.Ccol.temp
            tkrreplot(imgscree)
            if (EnScree.switch == "on") {
                tkrreplot(MGcomp$imgES)
            }
        }
        tkplace(tkbutton(ScPlotOps, text = "Change", width = 15, 
            command = function() ChangeScPoints()), relx = 0.325, 
            rely = 0.85, `in` = frameScPoints)
        ScreeLines <- tk2frame(ScPlotOpsNB)
        tkadd(ScPlotOpsNB, ScreeLines, text = "Lines")
        frameScLines <- tkwidget(ScreeLines, "TitleFrame", text = "Scree Plot Lines", 
            background = "white")
        tkplace(frameScLines, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ScreeLines)
        tkplace(tklabel(frameScLines, text = "Display Scree Curve", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameScLines)
        Scree.Curve.CB <- tk2checkbutton(ScPlotOpsNB)
        tkconfigure(Scree.Curve.CB, variable = MGvar$Scree.Curve.val)
        tkplace(Scree.Curve.CB, relx = 0.7, rely = 0.1, `in` = frameScLines)
        tkplace(tklabel(frameScLines, text = "Scree Curve Type", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameScLines)
        Scree.LineT.ComboBox <- tkwidget(ScPlotOps, "ComboBox", 
            editable = FALSE, values = c("Solid Line", "Dashed", 
                "Dotted"), width = 12, textvariable = MGvar$Scree.LineT.val)
        tkplace(Scree.LineT.ComboBox, relx = 0.6, rely = 0.25, 
            `in` = frameScLines)
        tkplace(tklabel(frameScLines, text = "Scree Curve Colour", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameScLines)
        ChangeColLine <- function() {
            MGvar$screeplot.curvecol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$screeplot.curvecol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$screeplot.curvecol.temp) > 0) 
                tkconfigure(ScreeColCurve, bg = MGvar$screeplot.curvecol.temp)
        }
        ScreeColCurve <- tkbutton(ScPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$screeplot.curvecol, command = function() ChangeColLine())
        tkplace(ScreeColCurve, relx = 0.7, rely = 0.4, `in` = frameScLines)
        tkplace(tklabel(frameScLines, text = "Display Current Dimension Line", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameScLines)
        Scree.CLine.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.CLine.CB, variable = MGvar$Scree.CLine.val)
        tkplace(Scree.CLine.CB, relx = 0.7, rely = 0.55, `in` = frameScLines)
        tkplace(tklabel(frameScLines, text = "Display Optimum Dimension Line", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameScLines)
        Scree.OLine.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.OLine.CB, variable = MGvar$Scree.OLine.val)
        tkplace(Scree.OLine.CB, relx = 0.7, rely = 0.7, `in` = frameScLines)
        ChangeScLines <- function() {
            SCurve = as.character(tclvalue(MGvar$Scree.Curve.val))
            if (SCurve == "1") {
                MGvar$screeplot.curve.show <<- "yes"
            }
            if (SCurve == "0") {
                MGvar$screeplot.curve.show <<- "no"
                tclvalue(MGvar$Scree.Curve.val) <<- 0
            }
            TypeL = as.character(tclvalue(MGvar$Scree.LineT.val))
            if (TypeL == "Solid Line") {
                MGvar$screeplot.curve.type <<- 1
                tclvalue(MGvar$Scree.LineT.val) <<- "Solid Line"
            }
            if (TypeL == "Dashed") {
                MGvar$screeplot.curve.type <<- 2
                tclvalue(MGvar$Scree.LineT.val) <<- "Dashed"
            }
            if (TypeL == "Dotted") {
                MGvar$screeplot.curve.type <<- 3
                tclvalue(MGvar$Scree.LineT.val) <<- "Dotted"
            }
            CLine = as.character(tclvalue(MGvar$Scree.CLine.val))
            if (CLine == "1") {
                MGvar$screeplot.Cline.show <<- "yes"
            }
            if (CLine == "0") {
                MGvar$screeplot.Cline.show <<- "no"
                tclvalue(MGvar$Scree.CLine.val) <<- 0
            }
            OLine = as.character(tclvalue(MGvar$Scree.OLine.val))
            if (OLine == "1") {
                MGvar$screeplot.Oline.show <<- "yes"
            }
            if (OLine == "0") {
                MGvar$screeplot.Oline.show <<- "no"
                tclvalue(MGvar$Scree.OLine.val) <<- 0
            }
            MGvar$screeplot.curvecol <<- MGvar$screeplot.curvecol.temp
            tkrreplot(imgscree)
            if (EnScree.switch == "on") {
                tkrreplot(MGcomp$imgES)
            }
        }
        tkplace(tkbutton(ScPlotOps, text = "Change", width = 15, 
            command = function() ChangeScLines()), relx = 0.325, 
            rely = 0.85, `in` = frameScLines)
        ScreeAxes <- tk2frame(ScPlotOpsNB)
        tkadd(ScPlotOpsNB, ScreeAxes, text = "Axes")
        frameScAxes <- tkwidget(ScreeAxes, "TitleFrame", text = "Scree Plot Axes", 
            background = "white")
        tkplace(frameScAxes, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ScreeAxes)
        tkplace(tklabel(frameScAxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameScAxes)
        Scree.Axes.CB <- tk2checkbutton(ScPlotOps)
        tkconfigure(Scree.Axes.CB, variable = MGvar$Scree.Axes.val)
        tkplace(Scree.Axes.CB, relx = 0.7, rely = 0.1, `in` = frameScAxes)
        tkplace(tklabel(frameScAxes, text = "Axes Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameScAxes)
        ChangeAxescol <- function() {
            MGvar$screeplot.Axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$screeplot.Axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$screeplot.Axescol.temp) > 0) 
                tkconfigure(ScreeColAxes, bg = MGvar$screeplot.Axescol.temp)
        }
        ScreeColAxes <- tkbutton(ScPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$screeplot.Axescol, command = function() ChangeAxescol())
        tkplace(ScreeColAxes, relx = 0.7, rely = 0.25, `in` = frameScAxes)
        ChangeScAxes <- function() {
            SAxes = as.character(tclvalue(MGvar$Scree.Axes.val))
            if (SAxes == "1") {
                MGvar$screeplot.Axes.xaxt <<- "s"
                MGvar$screeplot.Axes.yaxt <<- "s"
            }
            if (SAxes == "0") {
                MGvar$screeplot.Axes.xaxt <<- "n"
                MGvar$screeplot.Axes.yaxt <<- "n"
                tclvalue(MGvar$Scree.Axes.val) <<- 0
            }
            MGvar$screeplot.Axescol <<- MGvar$screeplot.Axescol.temp
            tkrreplot(imgscree)
            if (EnScree.switch == "on") {
                tkrreplot(MGcomp$imgES)
            }
        }
        tkplace(tkbutton(ScPlotOps, text = "Change", width = 15, 
            command = function() ChangeScAxes()), relx = 0.325, 
            rely = 0.85, `in` = frameScAxes)
        tkplace(ScPlotOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = ScPlotOps)
        OnOK.ScP <- function() {
            ChangeScGen()
            ChangeScPoints()
            ChangeScLines()
            ChangeScAxes()
            tkdestroy(ScPlotOps)
        }
        OnDefault.ScP <- function() {
            MGvar$screeplot.title.show <<- "yes"
            MGvar$screeplot.labs.show <<- "yes"
            MGvar$screeplot.leg.show <<- "yes"
            MGvar$screeplot.bg <<- "white"
            MGvar$screeplot.bg.temp <<- "white"
            MGvar$screeplot.points.show <<- "no"
            MGvar$screeplot.Cdim.show <<- "yes"
            MGvar$screeplot.Odim.show <<- "yes"
            MGvar$screeplot.Ccol <<- "red"
            MGvar$screeplot.Ccol.temp <<- "red"
            MGvar$screeplot.Ocol <<- "blue"
            MGvar$screeplot.Ocol.temp <<- "blue"
            MGvar$screeplot.curve.show <<- "yes"
            MGvar$screeplot.curve.type <<- 1
            MGvar$screeplot.curvecol <<- "black"
            MGvar$screeplot.curvecol.temp <<- "black"
            MGvar$screeplot.Cline.show <<- "yes"
            MGvar$screeplot.Oline.show <<- "yes"
            MGvar$screeplot.Axes.xaxt <<- "s"
            MGvar$screeplot.Axes.yaxt <<- "s"
            MGvar$screeplot.Axescol <<- "black"
            MGvar$screeplot.Axescol.temp <- "black"
            tclvalue(MGvar$Scree.Main.val) <<- 1
            tclvalue(MGvar$Scree.Lab.val) <<- 1
            tclvalue(MGvar$Scree.Leg.val) <<- 1
            tclvalue(MGvar$Scree.Points.val) <<- 0
            tclvalue(MGvar$Scree.CDim.val) <<- 1
            tclvalue(MGvar$Scree.ODim.val) <<- 1
            tclvalue(MGvar$Scree.Curve.val) <<- 1
            tclvalue(MGvar$Scree.LineT.val) <<- "Solid Line"
            tclvalue(MGvar$Scree.CLine.val) <<- 1
            tclvalue(MGvar$Scree.OLine.val) <<- 1
            tclvalue(MGvar$Scree.Axes.val) <<- 1
            tkrreplot(imgscree)
            if (EnScree.switch == "on") {
                tkrreplot(MGcomp$imgES)
            }
            tkdestroy(ScPlotOps)
        }
        tkplace(tkbutton(ScPlotOps, text = "OK", width = 15, 
            command = function() OnOK.ScP()), relx = 0.15, rely = 0.9, 
            `in` = ScPlotOps)
        tkplace(tkbutton(ScPlotOps, text = "Default", width = 15, 
            command = function() OnDefault.ScP()), relx = 0.55, 
            rely = 0.9, `in` = ScPlotOps)
        tkfocus(ScPlotOps)
        tkwait.window(ScPlotOps)
    }
    EnlargedScree <- function() {
        ReusePOW()
        EnScree.switch <<- "on"
        MGcomp$EScree <<- tktoplevel()
        tkwm.geometry(MGcomp$EScree, "773x575")
        tkwm.title(MGcomp$EScree, "Enlarged Scree Plot")
        EScreecanvas = tkcanvas(MGcomp$EScree, width = 818, height = 580, 
            bg = col.sec)
        tkplace(EScreecanvas, `in` = MGcomp$EScree)
        ES.height = as.numeric(tclvalue(tkwinfo("height", MGcomp$EScree)))
        ES.width = as.numeric(tclvalue(tkwinfo("width", MGcomp$EScree)))
        WidthScale = ES.width/1.8
        dimrat = 1.8/1.3
        EShscale = ES.width/WidthScale
        ESvscale = EShscale/dimrat
        MGcomp$imgES <<- tkrplot(MGcomp$EScree, function() plotScree(MGvar$scree.stress, 
            MGvar$screepoints.current, MGvar$screepoints.best, 
            MGvar$MDS.dimensions, MGvar$Opt.dim), hscale = EShscale, 
            vscale = ESvscale)
        tkplace(MGcomp$imgES, relx = 0.05, rely = 0.02, relwidth = 0.9, 
            `in` = MGcomp$EScree)
        CopySToClip <- function() {
            tkrreplot(MGcomp$imgES)
        }
        tkplace(tkbutton(MGcomp$EScree, text = "Copy to Clipboard", 
            width = 20, command = function() CopySToClip()), 
            relx = 0.22, rely = 0.93, `in` = MGcomp$EScree)
        tkplace(tkbutton(MGcomp$EScree, text = "Plot Options", 
            width = 20, command = function() ScreePlotOptions()), 
            relx = 0.6, rely = 0.93, `in` = MGcomp$EScree)
        tkbind(MGcomp$imgES, "<Destroy>", function() {
            EnScree.switch <<- "off"
        })
        tkbind(MGcomp$EScree, "<Configure>", resize.EScree)
    }
    resize.EScree <- function() {
        ES.height <- as.numeric(tclvalue(tkwinfo("height", MGcomp$EScree)))
        ES.width <- as.numeric(tclvalue(tkwinfo("width", MGcomp$EScree)))
        WidthScale = 773/1.8
        EShscale <- ES.width/WidthScale
        dimrat = 1.8/1.3
        ESvscale <- EShscale/dimrat
        tkrreplot(MGcomp$imgES, hscale = EShscale, vscale = ESvscale)
    }
    ZoomPlotOps <- function() {
        ZPlotOps = tktoplevel()
        tkwm.resizable(ZPlotOps, "0", "0")
        tkwm.deiconify(ZPlotOps)
        tkwm.title(ZPlotOps, "Zoomed Configuration Plot Options")
        tkwm.geometry(ZPlotOps, "350x400")
        ZPlotCanvas = tkcanvas(ZPlotOps, width = 420, height = 450, 
            bg = col.sec)
        tkplace(ZPlotCanvas, relx = 0, rely = 0, `in` = ZPlotOps)
        ZPlotOpsNB <- tk2notebook(ZPlotOps, tabs = NULL)
        ZoomGen <- tk2frame(ZPlotOpsNB)
        tkadd(ZPlotOpsNB, ZoomGen, text = "General")
        frameZGen <- tkwidget(ZoomGen, "TitleFrame", text = "Zoomed MDS-Configuration General", 
            background = "white")
        tkplace(frameZGen, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheigh = 0.96, `in` = ZoomGen)
        tkplace(tklabel(frameZGen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameZGen)
        Zoom.Main.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Main.CB, variable = MGvar$Zoom.Main.val)
        tkplace(Zoom.Main.CB, relx = 0.7, rely = 0.1, `in` = frameZGen)
        tkplace(tklabel(frameZGen, text = "Display Distance Measure", 
            background = "white"), relx = 0.1, rely = 0.22, `in` = frameZGen)
        Zoom.Dist.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Dist.CB, variable = MGvar$Zoom.Dist.val)
        tkplace(Zoom.Dist.CB, relx = 0.7, rely = 0.22, `in` = frameZGen)
        tkplace(tklabel(frameZGen, text = "Display Legend", background = "white"), 
            relx = 0.1, rely = 0.34, `in` = frameZGen)
        Zoom.Leg.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Leg.CB, variable = MGvar$Zoom.Leg.val)
        tkplace(Zoom.Leg.CB, relx = 0.7, rely = 0.34, `in` = frameZGen)
        tkplace(tklabel(frameZGen, text = "Provide Label for Y-Axis", 
            background = "white"), relx = 0.1, rely = 0.46, `in` = frameZGen)
        Zoom.Ylab.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Ylab.CB, variable = MGvar$Zoom.Ylab.val)
        tkplace(Zoom.Ylab.CB, relx = 0.7, rely = 0.46, `in` = frameZGen)
        tkplace(tklabel(frameZGen, text = "Provide Label for X-Axis", 
            background = "white"), relx = 0.1, rely = 0.58, `in` = frameZGen)
        Zoom.Xlab.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Xlab.CB, variable = MGvar$Zoom.Xlab.val)
        tkplace(Zoom.Xlab.CB, relx = 0.7, rely = 0.58, `in` = frameZGen)
        tkplace(tklabel(frameZGen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameZGen)
        ChangeZColBG <- function() {
            MGvar$zoomedplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$zoomedplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$zoomedplot.bg.temp) > 0) 
                tkconfigure(ZoomColBG, bg = MGvar$zoomedplot.bg.temp)
        }
        ZoomColBG <- tkbutton(ZPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$zoomedplot.bg, command = function() ChangeZColBG())
        tkplace(ZoomColBG, relx = 0.7, rely = 0.7, `in` = frameZGen)
        ChangeZGen <- function() {
            MainC = as.character(tclvalue(MGvar$Zoom.Main.val))
            if (MainC == "0") {
                MGvar$zoomedplot.title.show <<- "no"
            }
            if (MainC == "1") {
                MGvar$zoomedplot.title.show <<- "yes"
            }
            DistC = as.character(tclvalue(MGvar$Zoom.Dist.val))
            if (DistC == "0") {
                MGvar$zoomedplot.distmeas <<- "no"
                tclvalue(MGvar$Zoom.Dist.val) <<- 0
            }
            if (DistC == "1") {
                MGvar$zoomedplot.distmeas <<- "yes"
            }
            LegC = as.character(tclvalue(MGvar$Zoom.Leg.val))
            if (LegC == "1") {
                MGvar$zoomedplot.showleg <<- "yes"
                tclvalue(MGvar$Zoom.Leg.val) <<- 1
            }
            if (LegC == "0") {
                MGvar$zoomedplot.showleg <<- "no"
                tclvalue(MGvar$Zoom.Leg.val) <<- 0
            }
            YlabC = as.character(tclvalue(MGvar$Zoom.Ylab.val))
            if (YlabC == "0") {
                MGvar$zoomedplot.ylab <<- ""
            }
            if (YlabC == "1") {
                tclvalue(MGvar$Zoom.Ylab.val) <<- 1
                CYlabtt = tktoplevel()
                tkwm.resizable(CYlabtt, "0", "0")
                tkwm.deiconify(CYlabtt)
                tkwm.title(CYlabtt, "Configuration Plot Y Label")
                tkwm.geometry(CYlabtt, "310x100")
                CYlabcanvas = tkcanvas(CYlabtt, width = "310", 
                  height = "100", bg = col.sec)
                tkplace(CYlabcanvas, `in` = CYlabtt)
                frameCY <- tkwidget(CYlabtt, "TitleFrame", text = "Y Label Input", 
                  background = "white")
                tkplace(frameCY, relx = 0.01, relwidth = 0.98, 
                  rely = 0.01, relheight = 0.98, `in` = CYlabtt)
                tkplace(tklabel(frameCY, text = "Enter your Label for Y-Axis", 
                  background = "white"), relx = 0.05, rely = 0.25, 
                  `in` = frameCY)
                YInput = tclVar("")
                Ytext = tkentry(CYlabtt, width = 15, textvariable = YInput)
                tkplace(Ytext, relx = 0.6, rely = 0.25, `in` = frameCY)
                On.Enter <- function() {
                  tkdestroy(CYlabtt)
                  YL = as.character(tclvalue(YInput))
                  MGvar$zoomedplot.ylab <<- YL
                }
                tkplace(tkbutton(CYlabtt, width = 15, text = "Enter", 
                  command = function() On.Enter()), relx = 0.32, 
                  rely = 0.65, `in` = frameCY)
                tkbind(CYlabtt, "<Return>", On.Enter)
                tkwait.window(CYlabtt)
            }
            XlabC = as.character(tclvalue(MGvar$Zoom.Xlab.val))
            if (XlabC == "0") {
                MGvar$zoomedplot.xlab <<- ""
            }
            if (XlabC == "1") {
                tclvalue(MGvar$Conf.Xlab.val) <<- 1
                CXlabtt = tktoplevel()
                tkwm.resizable(CXlabtt, "0", "0")
                tkwm.deiconify(CXlabtt)
                tkwm.title(CXlabtt, "Configuration Plot X Label")
                tkwm.geometry(CXlabtt, "310x100")
                CXlabcanvas = tkcanvas(CXlabtt, width = "310", 
                  height = "100", bg = col.sec)
                tkplace(CXlabcanvas, `in` = CXlabtt)
                frameCX <- tkwidget(CXlabtt, "TitleFrame", text = "X Label Input", 
                  background = "white")
                tkplace(frameCX, relx = 0.01, relwidth = 0.98, 
                  rely = 0.01, relheight = 0.98, `in` = CXlabtt)
                tkplace(tklabel(frameCX, text = "Enter your Label for X-Axis", 
                  background = "white"), relx = 0.05, rely = 0.25, 
                  `in` = frameCX)
                XInput = tclVar("")
                Xtext = tkentry(CXlabtt, width = 15, textvariable = XInput)
                tkplace(Xtext, relx = 0.6, rely = 0.25, `in` = frameCX)
                On.EnterX <- function() {
                  tkdestroy(CXlabtt)
                  XL = as.character(tclvalue(XInput))
                  MGvar$zoomedplot.xlab <<- XL
                }
                tkplace(tkbutton(CXlabtt, width = 15, text = "Enter", 
                  command = function() On.EnterX()), relx = 0.32, 
                  rely = 0.65, `in` = frameCX)
                tkbind(CXlabtt, "<Return>", On.EnterX)
                tkwait.window(CXlabtt)
            }
            MGvar$zoomedplot.bg <<- MGvar$zoomedplot.bg.temp
            if (popzoomswitch == "on") {
                tkrreplot(MGcomp$zooming)
            }
            if (MGvar$seczoomswitch == "on") {
                tkrreplot(imgseczoom)
            }
        }
        tkplace(tkbutton(ZPlotOps, text = "Change", width = 15, 
            command = function() ChangeZGen()), relx = 0.325, 
            rely = 0.85, `in` = frameZGen)
        ZoomPoints <- tk2frame(ZPlotOpsNB)
        tkadd(ZPlotOpsNB, ZoomPoints, text = "Points")
        frameZPoints <- tkwidget(ZoomPoints, "TitleFrame", text = "MDS Configuration Points", 
            background = "white")
        tkplace(frameZPoints, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ZoomPoints)
        tkplace(tklabel(frameZPoints, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameZPoints)
        Zoom.Points.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Points.CB, variable = MGvar$Zoom.Points.val)
        tkplace(Zoom.Points.CB, relx = 0.7, rely = 0.1, `in` = frameZPoints)
        tkplace(tklabel(frameZPoints, text = "Display Point Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameZPoints)
        Zoom.Labels.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.Labels.CB, variable = MGvar$Zoom.Labels.val)
        tkplace(Zoom.Labels.CB, relx = 0.7, rely = 0.25, `in` = frameZPoints)
        tkplace(tklabel(frameZPoints, text = "Point Size", background = "white"), 
            relx = 0.1, rely = 0.4, `in` = frameZPoints)
        Zoom.PS.var <- tclVar(MGvar$zoomedplot.cex)
        Zoom.PS.spin <- tk2spinbox(ZPlotOps, from = 0.1, to = 2, 
            increment = 0.1, width = 12)
        tkconfigure(Zoom.PS.spin, textvariable = Zoom.PS.var)
        tkplace(Zoom.PS.spin, relx = 0.6, rely = 0.4, `in` = frameZPoints)
        tkplace(tklabel(frameZPoints, text = "Point Type", background = "white"), 
            relx = 0.1, rely = 0.55, `in` = frameZPoints)
        shapes <- c("Empty Boxes", "Filled Boxes", "Crosses", 
            "Empty Triangles", "Filled Triangles", "Filled Circles", 
            "Empty Circles")
        Zoom.PT.ComboBox <- tkwidget(ZPlotOps, "ComboBox", editable = FALSE, 
            values = shapes, width = 12, textvariable = MGvar$Zoom.PT.var)
        tkplace(Zoom.PT.ComboBox, relx = 0.6, rely = 0.55, `in` = frameZPoints)
        tkplace(tklabel(frameZPoints, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = frameZPoints)
        ChangeZColPoints <- function() {
            MGvar$zoomedplot.pointcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$zoomedplot.pointcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$zoomedplot.pointcol.temp) > 0) 
                tkconfigure(ZoomColPoints, bg = MGvar$zoomedplot.pointcol.temp)
        }
        ZoomColPoints <- tkbutton(ZPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$zoomedplot.pointcol, command = function() ChangeZColPoints())
        tkplace(ZoomColPoints, relx = 0.7, rely = 0.7, `in` = frameZPoints)
        ChangeZPoints <- function() {
            DispP = as.character(tclvalue(MGvar$Zoom.Points.val))
            if (DispP == "0") {
                MGvar$zoomedplot.showpoints <<- "no"
            }
            if (DispP == "1") {
                MGvar$zoomedplot.showpoints <<- "yes"
                tclvalue(MGvar$Zoom.Points.val) <<- 1
            }
            DispL = as.character(tclvalue(MGvar$Zoom.Labels.val))
            if (DispL == "0") {
                MGvar$zoomedplot.labs <<- "no"
                tclvalue(MGvar$Zoom.Labels.val) <<- 0
            }
            if (DispL == "1") {
                MGvar$zoomedplot.labs <<- "yes"
            }
            SizeP = as.numeric(tclvalue(MGvar$Zoom.PS.var))
            MGvar$zoomedplot.cex <<- SizeP
            TypeP = as.character(tclvalue(MGvar$Zoom.PT.var))
            if (TypeP == "Empty Circles") {
                MGvar$zoomedplot.type <<- 1
                tclvalue(MGvar$Zoom.PT.var) <<- "Empty Circles"
            }
            if (TypeP == "Filled Circles") {
                MGvar$zoomedplot.type <<- 16
                tclvalue(MGvar$Zoom.PT.var) <<- "Filled Circles"
            }
            if (TypeP == "Empty Boxes") {
                MGvar$zoomedplot.type <<- 22
                tclvalue(MGvar$Zoom.PT.var) <<- "Empty Boxes"
            }
            if (TypeP == "Filled Boxes") {
                MGvar$zoomedplot.type <<- 15
                tclvalue(MGvar$Zoom.PT.var) <<- "Filled Boxes"
            }
            if (TypeP == "Crosses") {
                MGvar$zoomedplot.type <<- 4
                tclvalue(MGvar$Zoom.PT.var) <<- "Crosses"
            }
            if (TypeP == "Empty Triangles") {
                MGvar$zoomedplot.type <<- 24
                tclvalue(MGvar$Zoom.PT.var) <<- "Empty Triangles"
            }
            if (TypeP == "Filled Triangles") {
                MGvar$zoomedplot.type <<- 17
                tclvalue(MGvar$Zoom.PT.var) <<- "Filled Triangles"
            }
            MGvar$zoomedplot.pointcol <<- MGvar$zoomedplot.pointcol.temp
            if (popzoomswitch == "on") {
                tkrreplot(MGcomp$zooming)
            }
            if (MGvar$seczoomswitch == "on") {
                tkrreplot(imgseczoom)
            }
        }
        tkplace(tkbutton(ZPlotOps, text = "Change", width = 15, 
            command = function() ChangeZPoints()), relx = 0.325, 
            rely = 0.85, `in` = frameZPoints)
        ZoomLines <- tk2frame(ZPlotOpsNB)
        tkadd(ZPlotOpsNB, ZoomLines, text = "Lines")
        frameZLines <- tkwidget(ZoomLines, "TitleFrame", text = "MDS Configuration Lines", 
            background = "white")
        tkplace(frameZLines, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ZoomLines)
        tkplace(tklabel(frameZLines, text = "Display Regression Axes", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameZLines)
        Zoom.RegLine.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.RegLine.CB, variable = MGvar$Zoom.RegLine.val)
        tkplace(Zoom.RegLine.CB, relx = 0.7, rely = 0.1, `in` = frameZLines)
        tkplace(tklabel(frameZLines, text = "Regression Axes Colour", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameZLines)
        ZoomChangeColReg <- function() {
            MGvar$zoomedplot.regcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$zoomedplot.regcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$zoomedplot.regcol.temp) > 0) 
                tkconfigure(ZoomColReg, bg = MGvar$zoomedplot.regcol.temp)
        }
        ZoomColReg <- tkbutton(ZPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$zoomedplot.regcol, command = function() ZoomChangeColReg())
        tkplace(ZoomColReg, relx = 0.7, rely = 0.25, `in` = frameZLines)
        ZoomChangeLines <- function() {
            ZShowR = as.character(tclvalue(MGvar$Zoom.RegLine.val))
            if (ZShowR == "1") {
                MGvar$zoomedplot.showreg <<- "yes"
                tclvalue(MGvar$Zoom.RegLine.val) <<- 1
            }
            if (ZShowR == "0") {
                MGvar$zoomedplot.showreg <<- "no"
                tclvalue(MGvar$Zoom.RegLine.val) <<- 0
            }
            MGvar$zoomedplot.regcol <<- MGvar$zoomedplot.regcol.temp
            if (popzoomswitch == "on") {
                tkrreplot(MGcomp$zooming)
            }
            if (MGvar$seczoomswitch == "on") {
                tkrreplot(imgseczoom)
            }
        }
        tkplace(tkbutton(ZPlotOps, text = "Change", width = 15, 
            command = function() ZoomChangeLines()), relx = 0.325, 
            rely = 0.85, `in` = frameZLines)
        ZoomAxes <- tk2frame(ZPlotOpsNB)
        tkadd(ZPlotOpsNB, ZoomAxes, text = "Axes")
        frameZAxes <- tkwidget(ZoomAxes, "TitleFrame", text = "MDS Configuration Axes", 
            background = "white")
        tkplace(frameZAxes, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ZoomAxes)
        tkplace(tklabel(frameZAxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameZAxes)
        Zoom.AxesMeas.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.AxesMeas.CB, variable = MGvar$Zoom.AxesMeas.val)
        tkplace(Zoom.AxesMeas.CB, relx = 0.7, rely = 0.1, `in` = frameZAxes)
        tkplace(tklabel(frameZAxes, text = "Axes Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = frameZAxes)
        ChangeZColAxes <- function() {
            MGvar$zoomedplot.axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$zoomedplot.axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$zoomedplot.axescol.temp) > 0) 
                tkconfigure(ZoomColAxes, bg = MGvar$zoomedplot.axescol.temp)
        }
        ZoomColAxes <- tkbutton(ZPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$zoomedplot.axescol, command = function() ChangeZColAxes())
        tkplace(ZoomColAxes, relx = 0.7, rely = 0.25, `in` = frameZAxes)
        tkplace(tklabel(frameZAxes, text = "Display Amount Zoomed", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameZAxes)
        Zoom.ShowZoom.CB <- tk2checkbutton(ZPlotOps)
        tkconfigure(Zoom.ShowZoom.CB, variable = MGvar$Zoom.ShowZoom.val)
        tkplace(Zoom.ShowZoom.CB, relx = 0.7, rely = 0.4, `in` = frameZAxes)
        ChangeZAxes <- function() {
            ShowA = as.character(tclvalue(MGvar$Zoom.AxesMeas.val))
            if (ShowA == "0") {
                MGvar$zoomedplot.yaxt <<- "n"
                MGvar$zoomedplot.yaxt <<- "n"
            }
            if (ShowA == "1") {
                MGvar$zoomedplot.yaxt <<- "s"
                MGvar$zoomedplot.yaxt <<- "s"
                tclvalue(MGvar$Zoom.AxesMeas.val) <<- 1
            }
            MGvar$zoomedplot.axescol <<- MGvar$zoomedplot.axescol.temp
            SZoom = as.character(tclvalue(MGvar$Zoom.ShowZoom.val))
            if (SZoom == "1") {
                MGvar$zoomedplot.showzoom <<- "yes"
                tclvalue(MGvar$Zoom.ShowZoom.val) <<- 1
            }
            if (SZoom == "0") {
                MGvar$zoomedplot.showzoom <<- "no"
                tclvalue(MGvar$Zoom.ShowZoom.val) <<- 0
            }
            if (popzoomswitch == "on") {
                tkrreplot(MGcomp$zooming)
            }
            if (MGvar$seczoomswitch == "on") {
                tkrreplot(imgseczoom)
            }
        }
        tkplace(tkbutton(ZPlotOps, text = "Change", width = 15, 
            command = function() ChangeZAxes()), relx = 0.325, 
            rely = 0.85, `in` = frameZAxes)
        tkplace(ZPlotOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = ZPlotOps)
        OnOK.CP <- function() {
            ChangeZGen()
            ChangeZPoints()
            ZoomChangeLines()
            ChangeZAxes()
            tkdestroy(ZPlotOps)
        }
        OnCancel.CP <- function() {
            tkdestroy(ZPlotOps)
        }
        OnDefault <- function() {
            MGvar$zoomedplot.title.show <<- "yes"
            MGvar$zoomedplot.distmeas <<- "yes"
            MGvar$zoomedplot.xlab <<- ""
            MGvar$zoomedplot.ylab <<- ""
            MGvar$zoomedplot.bg <<- "white"
            MGvar$zoomedplot.bg.temp <<- "white"
            MGvar$zoomedplot.cex <<- 1
            MGvar$zoomedplot.labs <<- "yes"
            MGvar$zoomedplot.showpoints <<- "no"
            MGvar$zoomedplot.pointcol <<- "black"
            MGvar$zoomedplot.pointcol.temp <<- "black"
            MGvar$zoomedplot.type <<- "1"
            MGvar$zoomedplot.yaxt <<- "n"
            MGvar$zoomedplot.yaxt <<- "n"
            MGvar$zoomedplot.axescol <<- "black"
            MGvar$zoomedplot.axescol.temp <<- "black"
            MGvar$zoomedplot.showreg <<- "no"
            tclvalue(MGvar$Zoom.Main.val) <<- 1
            tclvalue(MGvar$Zoom.Dist.val) <<- 1
            tclvalue(MGvar$Zoom.Ylab.val) <<- 0
            tclvalue(MGvar$Zoom.Xlab.val) <<- 0
            tclvalue(MGvar$Zoom.Points.val) <<- 0
            tclvalue(MGvar$Zoom.Labels.val) <<- 1
            tclvalue(MGvar$Zoom.PT.var) <<- "Empty Circles"
            tclvalue(MGvar$Zoom.AxesMeas.val) <<- 0
            tclvalue(MGvar$Zoom.RegLine.val) <<- 0
            if (popzoomswitch == "on") {
                tkrreplot(MGcomp$zooming)
            }
            MGvar$zoomedplot.cex <<- 0.7
            if (MGvar$seczoomswitch == "on") {
                tkrreplot(imgseczoom)
            }
            tkdestroy(ZPlotOps)
        }
        tkplace(tkbutton(ZPlotOps, text = "OK", width = 15, command = function() OnOK.CP()), 
            relx = 0.15, rely = 0.9, `in` = ZPlotOps)
        tkplace(tkbutton(ZPlotOps, text = "Default", width = 15, 
            command = function() OnDefault()), relx = 0.55, rely = 0.9, 
            `in` = ZPlotOps)
        tkfocus(ZPlotOps)
        tkwait.window(ZPlotOps)
    }
    ProcPlotOps <- function() {
        PPlotOps = tktoplevel()
        tkwm.resizable(PPlotOps, "0", "0")
        tkwm.deiconify(PPlotOps)
        tkwm.title(PPlotOps, "Procrustes Plotting Options")
        tkwm.geometry(PPlotOps, "350x400")
        PPlotcanvas = tkcanvas(PPlotOps, width = "420", height = "450", 
            bg = col.sec)
        tkplace(PPlotcanvas, relx = 0, rely = 0, `in` = PPlotOps)
        PPlotOpsNB <- tk2notebook(PPlotOps, tabs = NULL)
        ProcGen <- tk2frame(PPlotOpsNB)
        tkadd(PPlotOpsNB, ProcGen, text = "General")
        framePGen <- tkwidget(ProcGen, "TitleFrame", text = "Procrustes Plot General", 
            background = "white")
        tkplace(framePGen, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ProcGen)
        tkplace(tklabel(framePGen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framePGen)
        Proc.Main.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Main.CB, variable = MGvar$Proc.Main.val)
        tkplace(Proc.Main.CB, relx = 0.7, rely = 0.1, `in` = framePGen)
        tkplace(tklabel(framePGen, text = "Display Plot Legend", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framePGen)
        Proc.Leg.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Leg.CB, variable = MGvar$Proc.Leg.val)
        tkplace(Proc.Leg.CB, relx = 0.7, rely = 0.25, `in` = framePGen)
        tkplace(tklabel(framePGen, text = "Provide Label for Y-Axis", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = framePGen)
        Proc.Ylab.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Ylab.CB, variable = MGvar$Proc.Ylab.val)
        tkplace(Proc.Ylab.CB, relx = 0.7, rely = 0.4, `in` = framePGen)
        tkplace(tklabel(framePGen, text = "Provide Label for X-Axis", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = framePGen)
        Proc.Xlab.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Xlab.CB, variable = MGvar$Proc.Xlab.val)
        tkplace(Proc.Xlab.CB, relx = 0.7, rely = 0.55, `in` = framePGen)
        tkplace(tklabel(framePGen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = framePGen)
        ChangeProcColBG <- function() {
            MGvar$procplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$procplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$procplot.bg.temp) > 0) 
                tkconfigure(ProcColBG, bg = MGvar$procplot.bg.temp)
        }
        ProcColBG <- tkbutton(PPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$procplot.bg, command = function() ChangeProcColBG())
        tkplace(ProcColBG, relx = 0.7, rely = 0.7, `in` = framePGen)
        ChangeProcGen <- function() {
            MainP = as.character(tclvalue(MGvar$Proc.Main.val))
            if (MainP == "0") {
                MGvar$procplot.title.show <<- "no"
                tclvalue(MGvar$Proc.Main.val) <<- 0
            }
            if (MainP == "1") {
                MGvar$procplot.title.show <<- "yes"
                tclvalue(MGvar$Proc.Main.val) <<- 1
            }
            LegP = as.character(tclvalue(MGvar$Proc.Leg.val))
            if (LegP == "0") {
                MGvar$procplot.leg.show <<- "no"
                tclvalue(MGvar$Proc.Leg.val) <<- 0
            }
            if (LegP == "1") {
                MGvar$procplot.leg.show <<- "yes"
                tclvalue(MGvar$Proc.Leg.val) <<- 1
            }
            YlabP = as.character(tclvalue(MGvar$Proc.Ylab.val))
            if (YlabP == "0") {
                MGvar$procplot.ylab <<- ""
                tclvalue(MGvar$Proc.Ylab.val) <<- 0
            }
            if (YlabP == "1") {
                tclvalue(MGvar$Proc.Ylab.val) <<- 1
                PYlabtt = tktoplevel()
                tkwm.resizable(PYlabtt, "0", "0")
                tkwm.deiconify(PYlabtt)
                tkwm.title(PYlabtt, "Procrustes Plot Y Label")
                tkwm.geometry(PYlabtt, "310x100")
                PYlabcanvas = tkcanvas(PYlabtt, width = "310", 
                  height = "100", bg = col.sec)
                tkplace(PYlabcanvas, `in` = PYlabtt)
                framePY <- tkwidget(PYlabtt, "TitleFrame", text = "Y Label Input", 
                  background = "white")
                tkplace(framePY, relx = 0.01, relwidth = 0.98, 
                  rely = 0.01, relheight = 0.98, `in` = PYlabtt)
                tkplace(tklabel(framePY, text = "Enter your Label for Y-Axis", 
                  background = "white"), relx = 0.05, rely = 0.25, 
                  `in` = framePY)
                YInput = tclVar("")
                Ytext = tkentry(PYlabtt, width = 15, textvariable = YInput)
                tkplace(Ytext, relx = 0.6, rely = 0.25, `in` = framePY)
                On.Enter <- function() {
                  tkdestroy(PYlabtt)
                  YL = as.character(tclvalue(YInput))
                  MGvar$procplot.ylab <<- YL
                }
                tkplace(tkbutton(PYlabtt, width = 15, text = "Enter", 
                  command = function() On.Enter()), relx = 0.32, 
                  rely = 0.65, `in` = framePY)
                tkbind(PYlabtt, "<Return>", On.Enter)
                tkwait.window(PYlabtt)
            }
            XlabP = as.character(tclvalue(MGvar$Proc.Xlab.val))
            if (XlabP == "0") {
                MGvar$procplot.xlab <<- ""
                tclvalue(MGvar$Proc.Xlab.val) <<- 0
            }
            if (XlabP == "1") {
                tclvalue(MGvar$Proc.Xlab.val) <<- 1
                PXlabtt = tktoplevel()
                tkwm.resizable(PXlabtt, "0", "0")
                tkwm.deiconify(PXlabtt)
                tkwm.title(PXlabtt, "Procrustes Plot X Label")
                tkwm.geometry(PXlabtt, "310x100")
                PXlabcanvas = tkcanvas(PXlabtt, width = "310", 
                  height = "100", bg = col.sec)
                tkplace(PXlabcanvas, `in` = PXlabtt)
                framePX <- tkwidget(PXlabtt, "TitleFrame", text = "X Label Input", 
                  background = "white")
                tkplace(framePX, relx = 0.01, relwidth = 0.98, 
                  rely = 0.01, relheight = 0.98, `in` = PXlabtt)
                tkplace(tklabel(framePX, text = "Enter your Label for X-Axis", 
                  background = "white"), relx = 0.05, rely = 0.25, 
                  `in` = framePX)
                XInput = tclVar("")
                Xtext = tkentry(PXlabtt, width = 15, textvariable = XInput)
                tkplace(Xtext, relx = 0.6, rely = 0.25, `in` = framePX)
                On.Enter <- function() {
                  tkdestroy(PXlabtt)
                  XL = as.character(tclvalue(XInput))
                  MGvar$procplot.xlab <<- XL
                }
                tkplace(tkbutton(PXlabtt, width = 15, text = "Enter", 
                  command = function() On.Enter()), relx = 0.32, 
                  rely = 0.65, `in` = framePX)
                tkbind(PXlabtt, "<Return>", On.Enter)
                tkwait.window(PXlabtt)
            }
            MGvar$procplot.bg <<- MGvar$procplot.bg.temp
            tkrreplot(procimg)
            if (MGvar$EnProcPlot.switch == "on") {
                tkrreplot(MGcomp$POprocimg)
            }
        }
        tkplace(tkbutton(PPlotOps, text = "Change", width = 15, 
            command = function() ChangeProcGen()), relx = 0.325, 
            rely = 0.85, `in` = framePGen)
        ProcPoints1 <- tk2frame(PPlotOps)
        tkadd(PPlotOpsNB, ProcPoints1, text = "Points 1")
        framePPoints1 <- tkwidget(ProcPoints1, "TitleFrame", 
            text = "Procrustes Plot Points", background = "white")
        tkplace(framePPoints1, relx = 0.02, relwidth = 0.96, 
            rely = 0.02, relheight = 0.96, `in` = ProcPoints1)
        tkplace(tklabel(framePPoints1, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framePPoints1)
        Proc.Points1.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Points1.CB, variable = MGvar$Proc.Points1.val)
        tkplace(Proc.Points1.CB, relx = 0.7, rely = 0.1, `in` = framePPoints1)
        tkplace(tklabel(framePPoints1, text = "Display Point Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framePPoints1)
        Proc.Labels1.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Labels1.CB, variable = MGvar$Proc.Labels1.val)
        tkplace(Proc.Labels1.CB, relx = 0.7, rely = 0.25, `in` = framePPoints1)
        tkplace(tklabel(framePPoints1, text = "Point Size", background = "white"), 
            relx = 0.1, rely = 0.4, `in` = framePPoints1)
        Proc.PS1.val <- tclVar(MGvar$procplot.cex1)
        Proc.PS1.spin <- tk2spinbox(PPlotOps, from = 0.1, to = 2, 
            increment = 0.1, width = 12)
        tkconfigure(Proc.PS1.spin, textvariable = Proc.PS1.val)
        tkplace(Proc.PS1.spin, relx = 0.6, rely = 0.4, `in` = framePPoints1)
        tkplace(tklabel(framePPoints1, text = "Point Type", background = "white"), 
            relx = 0.1, rely = 0.55, `in` = framePPoints1)
        shapes <- c("Empty Boxes", "Filled Boxes", "Crosses", 
            "Empty Triangles", "Filled Triangles", "Filled Circles", 
            "Empty Circles")
        Proc.PT1.ComboBox <- tkwidget(PPlotOps, "ComboBox", editable = FALSE, 
            values = shapes, width = 12, textvariable = MGvar$Proc.PT1.val)
        tkplace(Proc.PT1.ComboBox, relx = 0.6, rely = 0.55, `in` = framePPoints1)
        tkplace(tklabel(framePPoints1, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = framePPoints1)
        ChangeProcCol1Points <- function() {
            MGvar$procplot.point1col.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$procplot.point1col.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$procplot.point1col.temp) > 0) 
                tkconfigure(ProcColPoints1, bg = MGvar$procplot.point1col.temp)
        }
        ProcColPoints1 <- tkbutton(PPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$procplot.point1col, command = function() ChangeProcCol1Points())
        tkplace(ProcColPoints1, relx = 0.7, rely = 0.7, `in` = framePPoints1)
        ChangePoints1 <- function() {
            DispP1 = as.character(tclvalue(MGvar$Proc.Points1.val))
            if (DispP1 == "0") {
                MGvar$procplot.showpoints1 <<- "no"
                tclvalue(MGvar$Proc.Points1.val) <<- 0
            }
            if (DispP1 == "1") {
                MGvar$procplot.showpoints1 <<- "yes"
                tclvalue(MGvar$Proc.Points1.val) <<- 1
            }
            DispL1 = as.character(tclvalue(MGvar$Proc.Labels1.val))
            if (DispL1 == "0") {
                MGvar$procplot.labs1 <<- "no"
                tclvalue(MGvar$Proc.Labels1.val) <<- 0
            }
            if (DispL1 == "1") {
                MGvar$procplot.labs1 <<- "yes"
                tclvalue(MGvar$Proc.Labels1.val) <<- 1
            }
            SizeP1 = as.numeric(tclvalue(Proc.PS1.val))
            MGvar$procplot.cex1 <<- SizeP1
            TypeP1 = as.character(tclvalue(MGvar$Proc.PT1.val))
            if (TypeP1 == "Empty Circles") {
                MGvar$procplot.type1 <<- 1
                tclvalue(MGvar$Proc.PT1.val) <<- "Empty Circles"
            }
            if (TypeP1 == "Filled Circles") {
                MGvar$procplot.type1 <<- 16
                tclvalue(MGvar$Proc.PT1.val) <<- "Filled Circles"
            }
            if (TypeP1 == "Empty Boxes") {
                MGvar$procplot.type1 <<- 22
                tclvalue(MGvar$Proc.PT1.val) <<- "Empty Boxes"
            }
            if (TypeP1 == "Filled Boxes") {
                MGvar$procplot.type1 <<- 15
                tclvalue(MGvar$Proc.PT1.val) <<- "Filled Boces"
            }
            if (TypeP1 == "Crosses") {
                MGvar$procplot.type1 <<- 4
                tclvalue(MGvar$Proc.PT1.val) <<- "Crosses"
            }
            if (TypeP1 == "Empty Triangles") {
                MGvar$procplot.type1 <<- 24
                tclvalue(MGvar$Proc.PT1.val) <<- "Empty Triangles"
            }
            if (TypeP1 == "Filled Triangles") {
                MGvar$procplot.type1 <<- 17
                tclvalue(MGvar$Proc.PT1.val) <<- "Filled Triangles"
            }
            MGvar$procplot.point1col <<- MGvar$procplot.point1col.temp
            tkrreplot(procimg)
            if (MGvar$EnProcPlot.switch == "on") {
                tkrreplot(MGcomp$POprocimg)
            }
        }
        tkplace(tkbutton(PPlotOps, text = "Change", width = 15, 
            command = function() ChangePoints1()), relx = 0.325, 
            rely = 0.85, `in` = framePPoints1)
        ProcPoints2 <- tk2frame(PPlotOps)
        tkadd(PPlotOpsNB, ProcPoints2, text = "Points 2")
        framePPoints2 <- tkwidget(ProcPoints2, "TitleFrame", 
            text = "Procrustes Plot Points", background = "white")
        tkplace(framePPoints2, relx = 0.02, relwidth = 0.96, 
            rely = 0.02, relheight = 0.96, `in` = ProcPoints2)
        tkplace(tklabel(framePPoints2, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framePPoints2)
        Proc.Points2.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Points2.CB, variable = MGvar$Proc.Points2.val)
        tkplace(Proc.Points2.CB, relx = 0.7, rely = 0.1, `in` = framePPoints2)
        tkplace(tklabel(framePPoints2, text = "Display Point Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framePPoints2)
        Proc.Labels2.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.Labels2.CB, variable = MGvar$Proc.Labels2.val)
        tkplace(Proc.Labels2.CB, relx = 0.7, rely = 0.25, `in` = framePPoints2)
        tkplace(tklabel(framePPoints2, text = "Point Size", background = "white"), 
            relx = 0.1, rely = 0.4, `in` = framePPoints2)
        Proc.PS2.val <- tclVar(MGvar$procplot.cex2)
        Proc.PS2.spin <- tk2spinbox(PPlotOps, from = 0.1, to = 2, 
            increment = 0.1, width = 12)
        tkconfigure(Proc.PS2.spin, textvariable = Proc.PS2.val)
        tkplace(Proc.PS2.spin, relx = 0.6, rely = 0.4, `in` = framePPoints2)
        tkplace(tklabel(framePPoints2, text = "Point Type", background = "white"), 
            relx = 0.1, rely = 0.55, `in` = framePPoints2)
        shapes <- c("Empty Boxes", "Filled Boxes", "Crosses", 
            "Empty Triangles", "Filled Triangles", "Filled Circles", 
            "Empty Circles")
        Proc.PT2.ComboBox <- tkwidget(PPlotOps, "ComboBox", editable = FALSE, 
            values = shapes, width = 12, textvariable = MGvar$Proc.PT2.val)
        tkplace(Proc.PT2.ComboBox, relx = 0.6, rely = 0.55, `in` = framePPoints2)
        tkplace(tklabel(framePPoints2, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = framePPoints2)
        ChangeProcCol2Points <- function() {
            MGvar$procplot.point2col.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$procplot.point2col.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$procplot.point2col.temp) > 0) 
                tkconfigure(ProcColPoints2, bg = MGvar$procplot.point2col.temp)
        }
        ProcColPoints2 <- tkbutton(PPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$procplot.point2col, command = function() ChangeProcCol2Points())
        tkplace(ProcColPoints2, relx = 0.7, rely = 0.7, `in` = framePPoints2)
        ChangePoints2 <- function() {
            DispP2 = as.character(tclvalue(MGvar$Proc.Points2.val))
            if (DispP2 == "0") {
                MGvar$procplot.showpoints2 <<- "no"
                tclvalue(MGvar$Proc.Points2.val) <<- 0
            }
            if (DispP2 == "1") {
                MGvar$procplot.showpoints2 <<- "yes"
                tclvalue(MGvar$Proc.Points2.val) <<- 1
            }
            DispL2 = as.character(tclvalue(MGvar$Proc.Labels2.val))
            if (DispL2 == "0") {
                MGvar$procplot.labs2 <<- "no"
                tclvalue(MGvar$Proc.Labels2.val) <<- 0
            }
            if (DispL2 == "1") {
                MGvar$procplot.labs2 <<- "yes"
                tclvalue(MGvar$Proc.Labels2.val) <<- 1
            }
            SizeP2 = as.numeric(tclvalue(Proc.PS2.val))
            MGvar$procplot.cex2 <<- SizeP2
            TypeP2 = as.character(tclvalue(MGvar$Proc.PT2.val))
            if (TypeP2 == "Empty Circles") {
                MGvar$procplot.type2 <<- 1
                tclvalue(MGvar$Proc.PT2.val) <<- "Empty Circles"
            }
            if (TypeP2 == "Filled Circles") {
                MGvar$procplot.type2 <<- 16
                tclvalue(MGvar$Proc.PT2.val) <<- "Filled Circles"
            }
            if (TypeP2 == "Empty Boxes") {
                MGvar$procplot.type2 <<- 22
                tclvalue(MGvar$Proc.PT2.val) <<- "Empty Boxes"
            }
            if (TypeP2 == "Filled Boxes") {
                MGvar$procplot.type2 <<- 15
                tclvalue(MGvar$Proc.PT2.val) <<- "Filled Boces"
            }
            if (TypeP2 == "Crosses") {
                MGvar$procplot.type2 <<- 4
                tclvalue(MGvar$Proc.PT2.val) <<- "Crosses"
            }
            if (TypeP2 == "Empty Triangles") {
                MGvar$procplot.type2 <<- 24
                tclvalue(MGvar$Proc.PT2.val) <<- "Empty Triangles"
            }
            if (TypeP2 == "Filled Triangles") {
                MGvar$procplot.type2 <<- 17
                tclvalue(MGvar$Proc.PT2.val) <<- "Filled Triangles"
            }
            MGvar$procplot.point2col <<- MGvar$procplot.point2col.temp
            tkrreplot(procimg)
            if (MGvar$EnProcPlot.switch == "on") {
                tkrreplot(MGcomp$POprocimg)
            }
        }
        tkplace(tkbutton(PPlotOps, text = "Change", width = 15, 
            command = function() ChangePoints2()), relx = 0.325, 
            rely = 0.85, `in` = framePPoints2)
        ProcLines <- tk2frame(PPlotOpsNB)
        tkadd(PPlotOpsNB, ProcLines, text = "Lines")
        framePLines <- tkwidget(ProcLines, "TitleFrame", text = "Procrustes Plot Lines", 
            background = "white")
        tkplace(framePLines, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ProcLines)
        tkplace(tklabel(framePLines, text = "Display Points 1 Regression Axes", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framePLines)
        Proc.RegLine1.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.RegLine1.CB, variable = MGvar$Proc.RegLine1.val)
        tkplace(Proc.RegLine1.CB, relx = 0.7, rely = 0.1, `in` = framePLines)
        tkplace(tklabel(framePLines, text = "Regression 1 Axes Colour", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framePLines)
        ChangeColPReg1 <- function() {
            MGvar$procplot.regcol1.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$procplot.regcol1.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$procplot.regcol1.temp) > 0) 
                tkconfigure(ProcColReg1, bg = MGvar$procplot.regcol1.temp)
        }
        ProcColReg1 <- tkbutton(PPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$procplot.regcol1, command = function() ChangeColPReg1())
        tkplace(ProcColReg1, relx = 0.7, rely = 0.25, `in` = framePLines)
        tkplace(tklabel(framePLines, text = "Display Points 2 Regression Axes", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = framePLines)
        Proc.RegLine2.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.RegLine2.CB, variable = MGvar$Proc.RegLine2.val)
        tkplace(Proc.RegLine2.CB, relx = 0.7, rely = 0.4, `in` = framePLines)
        tkplace(tklabel(framePLines, text = "Regression 2 Axes Colour", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = framePLines)
        ChangeColPReg2 <- function() {
            MGvar$procplot.regcol2.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$procplot.regcol2.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$procplot.regcol2.temp) > 0) 
                tkconfigure(ProcColReg2, bg = MGvar$procplot.regcol2.temp)
        }
        ProcColReg2 <- tkbutton(PPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$procplot.regcol2, command = function() ChangeColPReg2())
        tkplace(ProcColReg2, relx = 0.7, rely = 0.55, `in` = framePLines)
        ChangeProcLines <- function() {
            PShowR1 = as.character(tclvalue(MGvar$Proc.RegLine1.val))
            if (PShowR1 == "1") {
                MGvar$procplot.showreg1 <<- "yes"
                tclvalue(MGvar$Proc.RegLine1.val) <<- 1
            }
            if (PShowR1 == "0") {
                MGvar$procplot.showreg1 <<- "no"
                tclvalue(MGvar$Proc.RegLine1.val) <<- 0
            }
            PShowR2 = as.character(tclvalue(MGvar$Proc.RegLine2.val))
            if (PShowR2 == "1") {
                MGvar$procplot.showreg2 <<- "yes"
                tclvalue(MGvar$Proc.RegLine2.val) <<- 1
            }
            if (PShowR2 == "0") {
                MGvar$procplot.showreg2 <<- "no"
                tclvalue(MGvar$Proc.RegLine2.val) <<- 0
            }
            MGvar$procplot.regcol1 <<- MGvar$procplot.regcol1.temp
            MGvar$procplot.regcol2 <<- MGvar$procplot.regcol2.temp
            tkrreplot(procimg)
            if (MGvar$EnProcPlot.switch == "on") {
                tkrreplot(MGcomp$POprocimg)
            }
        }
        tkplace(tkbutton(PPlotOps, text = "Change", width = 15, 
            command = function() ChangeProcLines()), relx = 0.325, 
            rely = 0.85, `in` = framePLines)
        ProcAxes <- tk2frame(PPlotOpsNB)
        tkadd(PPlotOpsNB, ProcAxes, text = "Axes")
        framePAxes <- tkwidget(ProcAxes, "TitleFrame", text = "Procrustes Plot Axes", 
            background = "white")
        tkplace(framePAxes, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = ProcAxes)
        tkplace(tklabel(framePAxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framePAxes)
        Proc.AxesMeas.CB <- tk2checkbutton(PPlotOps)
        tkconfigure(Proc.AxesMeas.CB, variable = MGvar$Proc.AxesMeas.val)
        tkplace(Proc.AxesMeas.CB, relx = 0.7, rely = 0.1, `in` = framePAxes)
        tkplace(tklabel(framePAxes, text = "Axes Colour", background = "white"), 
            relx = 0.1, rely = 0.25, `in` = framePAxes)
        ChangeColPAxes <- function() {
            MGvar$procplot.axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$procplot.axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$procplot.axescol.temp) > 0) 
                tkconfigure(ProcColAxes, bg = MGvar$procplot.axescol.temp)
        }
        ProcColAxes <- tkbutton(PPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$procplot.axescol, command = function() ChangeColPAxes())
        tkplace(ProcColAxes, relx = 0.7, rely = 0.25, `in` = framePAxes)
        ChangePAxes <- function() {
            PShowA = as.character(tclvalue(MGvar$Proc.AxesMeas.val))
            if (PShowA == "0") {
                MGvar$procplot.yaxt <<- "n"
                MGvar$procplot.xaxt <<- "n"
                tclvalue(MGvar$Proc.AxesMeas.val) <<- 0
            }
            if (PShowA == "1") {
                MGvar$procplot.yaxt <<- "s"
                MGvar$procplot.xaxt <<- "s"
                tclvalue(MGvar$Proc.AxesMeas.val) <<- 1
            }
            MGvar$procplot.axescol <<- MGvar$procplot.axescol.temp
            tkrreplot(procimg)
            if (MGvar$EnProcPlot.switch == "on") {
                tkrreplot(MGcomp$POprocimg)
            }
        }
        tkplace(tkbutton(PPlotOps, text = "Change", width = 15, 
            command = function() ChangePAxes()), relx = 0.325, 
            rely = 0.85, `in` = framePAxes)
        tkplace(PPlotOpsNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = PPlotOps)
        OnOK.PP <- function() {
            ChangeProcGen()
            ChangePoints1()
            ChangePoints2()
            ChangeProcLines()
            ChangePAxes()
            tkdestroy(PPlotOps)
        }
        OnCancel.PP <- function() {
            tkdestroy(PPlotOps)
        }
        OnDefault.PP <- function() {
            MGvar$procplot.title.show <<- "yes"
            MGvar$procplot.leg.show <<- "no"
            MGvar$procplot.ylab <<- ""
            MGvar$procplot.xlab <<- ""
            MGvar$procplot.bg.temp <<- "white"
            MGvar$procplot.bg <<- "white"
            MGvar$procplot.showpoints1 <<- "no"
            MGvar$procplot.labs1 <<- "yes"
            MGvar$procplot.cex1 <<- 0.7
            MGvar$procplot.type1 <<- 1
            MGvar$procplot.point1col.temp <<- "red"
            MGvar$procplot.point1col <<- "red"
            MGvar$procplot.showpoints2 <<- "no"
            MGvar$procplot.labs2 <<- "yes"
            MGvar$procplot.cex2 <<- 0.7
            MGvar$procplot.type2 <<- 1
            MGvar$procplot.point2col.temp <<- "blue"
            MGvar$procplot.point2col <<- "blue"
            MGvar$procplot.yaxt <<- "n"
            MGvar$procplot.xaxt <<- "n"
            MGvar$procplot.axescol.temp <<- "black"
            MGvar$procplot.axescol <<- "black"
            MGvar$procplot.showreg1 <<- "no"
            MGvar$procplot.showreg2 <<- "no"
            MGvar$procplot.regcol1 <<- "red"
            MGvar$procplot.regcol2 <<- "blue"
            MGvar$Proc.Main.val <<- tclVar("1")
            MGvar$Proc.Leg.val <<- tclVar("0")
            MGvar$Proc.Ylab.val <<- tclVar("0")
            MGvar$Proc.Xlab.val <<- tclVar("0")
            MGvar$Proc.Points1.val <<- tclVar("0")
            MGvar$Proc.Labels1.val <<- tclVar("1")
            MGvar$Proc.PT1.val <<- tclVar("Empty Circles")
            MGvar$Proc.Points2.val <<- tclVar("0")
            MGvar$Proc.Labels2.val <<- tclVar("1")
            MGvar$Proc.PT2.val <<- tclVar("Empty Circles")
            MGvar$Proc.RegLine1.val <<- tclVar("0")
            MGvar$Proc.RegLine2.val <<- tclVar("0")
            MGvar$Proc.AxesMeas.val <<- tclVar("0")
            tkrreplot(procimg)
            tkdestroy(PPlotOps)
        }
        tkplace(tkbutton(PPlotOps, text = "OK", width = 15, command = function() OnOK.PP()), 
            relx = 0.15, rely = 0.9, `in` = PPlotOps)
        tkplace(tkbutton(PPlotOps, text = "Default", width = 15, 
            command = function() OnDefault.PP()), relx = 0.55, 
            rely = 0.9, `in` = PPlotOps)
        tkfocus(PPlotOps)
        tkwait.window(PPlotOps)
    }
    threeDPlotOptions <- function() {
        TDPlotOps <- tktoplevel()
        tkwm.resizable(TDPlotOps, "0", "0")
        tkwm.deiconify(TDPlotOps)
        tkwm.title(TDPlotOps, "3D Plotting Options")
        tkwm.geometry(TDPlotOps, "420x500")
        TDPlotcanvas = tkcanvas(TDPlotOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(TDPlotcanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = TDPlotOps)
        TDmainNB <- tk2notebook(TDPlotOps, tabs = NULL)
        StatOps <- tk2frame(TDmainNB)
        tkadd(TDmainNB, StatOps, text = "Static 3D")
        StatNB <- tk2notebook(TDPlotOps, tabs = NULL)
        Stat3DGen <- tk2frame(StatNB)
        tkadd(StatNB, Stat3DGen, text = "General")
        framestatgen <- tkwidget(Stat3DGen, "TitleFrame", text = "Static 3D Configuration General", 
            background = "white")
        tkplace(framestatgen, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = Stat3DGen)
        tkplace(tklabel(framestatgen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framestatgen)
        Stat.Main.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Main.CB, variable = MGvar$Stat.Main.val)
        tkplace(Stat.Main.CB, relx = 0.75, rely = 0.1, `in` = framestatgen)
        tkplace(tklabel(framestatgen, text = "Display Distance Measure", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framestatgen)
        Stat.Dist.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Dist.CB, variable = MGvar$Stat.Dist.val)
        tkplace(Stat.Dist.CB, relx = 0.75, rely = 0.25, `in` = framestatgen)
        tkplace(tklabel(framestatgen, text = "Display Legend", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = framestatgen)
        Stat.Leg.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Leg.CB, variable = MGvar$Stat.Leg.val)
        tkplace(Stat.Leg.CB, relx = 0.75, rely = 0.4, `in` = framestatgen)
        tkplace(tklabel(framestatgen, text = "Display Axes Labels", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = framestatgen)
        tkplace(tklabel(framestatgen, text = "Background Colour", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = framestatgen)
        Stat.Labs.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Labs.CB, variable = MGvar$Stat.Labs.val)
        tkplace(Stat.Labs.CB, relx = 0.75, rely = 0.55, `in` = framestatgen)
        ChangeStatColBG <- function() {
            MGvar$sthreeDplot.bg.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$sthreeDplot.bg.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$sthreeDplot.bg.temp) > 0) 
                tkconfigure(StatColBG, bg = MGvar$sthreeDplot.bg.temp)
        }
        StatColBG <- tkbutton(TDPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$sthreeDplot.bg, command = function() ChangeStatColBG())
        tkplace(StatColBG, relx = 0.75, rely = 0.7, `in` = framestatgen)
        ChangeStatGen <- function() {
            StMain = as.character(tclvalue(MGvar$Stat.Main.val))
            if (StMain == "0") {
                MGvar$sthreeDplot.title.show <<- "no"
                tclvalue(MGvar$Stat.Main.val) <<- 0
            }
            if (StMain == "1") {
                MGvar$sthreeDplot.title.show <<- "yes"
                tclvalue(MGvar$Stat.Main.val) <<- 1
            }
            StDist = as.character(tclvalue(MGvar$Stat.Dist.val))
            if (StDist == "0") {
                MGvar$sthreeDplot.distmeas <<- "no"
                tclvalue(MGvar$Stat.Dist.val) <<- 0
            }
            if (StDist == "1") {
                MGvar$sthreeDplot.distmeas <<- "yes"
                tclvalue(MGvar$Stat.Dist.val) <<- 1
            }
            StLeg = as.character(tclvalue(MGvar$Stat.Leg.val))
            if (StLeg == "0") {
                MGvar$sthreeDplot.leg.show <<- "no"
                tclvalue(MGvar$Stat.Leg.val) <<- 0
            }
            if (StLeg == "1") {
                MGvar$sthreeDplot.leg.show <<- "yes"
                tclvalue(MGvar$Stat.Leg.val) <<- 1
            }
            StLabs = as.character(tclvalue(MGvar$Stat.Labs.val))
            if (StLabs == "1") {
                MGvar$sthreeDplot.xlab <<- paste("Dim", MGvar$PlottingDimX.3D)
                MGvar$sthreeDplot.ylab <<- paste("Dim", MGvar$PlottingDimY.3D)
                MGvar$sthreeDplot.zlab <<- paste("Dim", MGvar$PlottingDimZ.3D)
                tclvalue(MGvar$Stat.Labs.val) <<- 1
            }
            if (StLabs == "0") {
                MGvar$sthreeDplot.xlab <<- ""
                MGvar$sthreeDplot.ylab <<- ""
                MGvar$sthreeDplot.zlab <<- ""
                tclvalue(MGvar$Stat.Labs.val) <<- 0
            }
            MGvar$sthreeDplot.bg <<- MGvar$sthreeDplot.bg.temp
            tkrreplot(img3Dstat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeStatGen()), relx = 0.325, 
            rely = 0.85, `in` = framestatgen)
        Stat3DPoints <- tk2frame(StatNB)
        tkadd(StatNB, Stat3DPoints, text = "Points")
        framestatpoints <- tkwidget(Stat3DPoints, "TitleFrame", 
            text = "Static 3D Configuration Points", background = "white")
        tkplace(framestatpoints, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = Stat3DPoints)
        tkplace(tklabel(framestatpoints, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framestatpoints)
        Stat.Points.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Points.CB, variable = MGvar$Stat.Points.val)
        tkplace(Stat.Points.CB, relx = 0.75, rely = 0.1, `in` = framestatpoints)
        tkplace(tklabel(framestatpoints, text = "Display Point Labels", 
            background = "white"), relx = 0.1, rely = 0.22, `in` = framestatpoints)
        Stat.Labels.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Labels.CB, variable = MGvar$Stat.Labels.val)
        tkplace(Stat.Labels.CB, relx = 0.75, rely = 0.22, `in` = framestatpoints)
        tkplace(tklabel(framestatpoints, text = "Display Point Size", 
            background = "white"), relx = 0.1, rely = 0.34, `in` = framestatpoints)
        Stat.PS.var <- tclVar(MGvar$sthreeDplot.cex)
        Stat.PS.spin <- tk2spinbox(TDPlotOps, from = 0.1, to = 2, 
            increment = 0.1, width = 15)
        tkconfigure(Stat.PS.spin, textvariable = Stat.PS.var)
        tkplace(Stat.PS.spin, relx = 0.65, rely = 0.34, `in` = framestatpoints)
        tkplace(tklabel(framestatpoints, text = "Point Type", 
            background = "white"), relx = 0.1, rely = 0.46, `in` = framestatpoints)
        shapes <- c("Empty Boxes", "Filled Boxes", "Crosses", 
            "Empty Triangles", "Filled Triangles", "Filled Circles", 
            "Empty Circles")
        Stat.PT.ComboBox <- tkwidget(TDPlotOps, "ComboBox", editable = FALSE, 
            values = shapes, width = 12, textvariable = MGvar$Stat.PT.var)
        tkplace(Stat.PT.ComboBox, relx = 0.65, rely = 0.46, `in` = framestatpoints)
        tkplace(tklabel(framestatpoints, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.58, `in` = framestatpoints)
        ChangeStatColPoints <- function() {
            MGvar$sthreeDplot.pointcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$sthreeDplot.pointcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$sthreeDplot.pointcol.temp) > 0) 
                tkconfigure(StatColPoints, bg = MGvar$sthreeDplot.pointcol.temp)
        }
        StatColPoints <- tkbutton(TDPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$sthreeDplot.pointcol, command = function() ChangeStatColPoints())
        tkplace(StatColPoints, relx = 0.75, rely = 0.58, `in` = framestatpoints)
        tkplace(tklabel(framestatpoints, text = "Highlight 3D element", 
            background = "white"), relx = 0.1, rely = 0.7, `in` = framestatpoints)
        Stat.HL.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.HL.CB, variable = MGvar$Stat.HL.var)
        tkplace(Stat.HL.CB, relx = 0.75, rely = 0.7, `in` = framestatpoints)
        ChangeStatPoints <- function() {
            StDispP = as.character(tclvalue(MGvar$Stat.Points.val))
            if (StDispP == "0") {
                MGvar$sthreeDplot.showpoints <<- "no"
                tclvalue(MGvar$Stat.Points.val) <<- 0
            }
            if (StDispP == "1") {
                MGvar$sthreeDplot.showpoints <<- "yes"
                tclvalue(MGvar$Stat.Points.val) <<- 1
            }
            StDispL = as.character(tclvalue(MGvar$Stat.Labels.val))
            if (StDispL == "0") {
                MGvar$sthreeDplot.showlabels <<- "no"
                tclvalue(MGvar$Stat.Labels.val) <<- 0
            }
            if (StDispL == "1") {
                MGvar$sthreeDplot.showlabels <<- "yes"
                tclvalue(MGvar$Stat.Labels.val) <<- 1
            }
            StSizeP = as.numeric(tclvalue(Stat.PS.var))
            MGvar$sthreeDplot.cex <<- StSizeP
            StTypeP = as.character(tclvalue(MGvar$Stat.PT.var))
            if (StTypeP == "Empty Circles") {
                MGvar$sthreeDplot.type <<- 1
                tclvalue(MGvar$Stat.PT.var) <<- "Empty Circles"
            }
            if (StTypeP == "Filled Circles") {
                MGvar$sthreeDplot.type <<- 16
                tclvalue(MGvar$Stat.PT.var) <<- "Filled Circles"
            }
            if (StTypeP == "Empty Boxes") {
                MGvar$sthreeDplot.type <<- 22
                tclvalue(MGvar$Stat.PT.var) <<- "Empty Boxes"
            }
            if (StTypeP == "Filled Boxes") {
                MGvar$sthreeDplot.type <<- 15
                tclvalue(MGvar$Stat.PT.var) <<- "Filled Boxes"
            }
            if (StTypeP == "Crosses") {
                MGvar$sthreeDplot.type <<- 4
                tclvalue(MGvar$Stat.PT.var) <<- "Crosses"
            }
            if (StTypeP == "Empty Triangles") {
                MGvar$sthreeDplot.type <<- 24
                tclvalue(MGvar$Stat.PT.var) <<- "Empty Triangles"
            }
            if (StTypeP == "Filled Triangles") {
                MGvar$sthreeDplot.type <<- 17
                tclvalue(MGvar$Stat.PT.var) <<- "Filled Triangles"
            }
            Highlight = as.character(tclvalue(MGvar$Stat.HL.var))
            if (Highlight == "0") {
                MGvar$sthreeDplot.HL <<- "no"
                tclvalue(MGvar$Stat.HL.var) <<- 0
            }
            if (Highlight == "1") {
                MGvar$sthreeDplot.HL <<- "yes"
                tclvalue(MGvar$Stat.HL.var) <<- 1
            }
            MGvar$sthreeDplot.pointcol <<- MGvar$sthreeDplot.pointcol.temp
            tkrreplot(img3Dstat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeStatPoints()), relx = 0.325, 
            rely = 0.85, `in` = framestatpoints)
        Stat3DLines <- tk2frame(StatNB)
        tkadd(StatNB, Stat3DLines, text = "Lines")
        framestatlines <- tkwidget(Stat3DLines, "TitleFrame", 
            text = "Static 3D Configuration Lines", background = "white")
        tkplace(framestatlines, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = Stat3DLines)
        tkplace(tklabel(framestatlines, text = "Display Regression on X Plane", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framestatlines)
        Stat.RegLine.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.RegLine.CB, variable = MGvar$Stat.RegLine.val)
        tkplace(Stat.RegLine.CB, relx = 0.75, rely = 0.1, `in` = framestatlines)
        tkplace(tklabel(framestatlines, text = "Regression X Plane Colour", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framestatlines)
        ChangeStatColReg <- function() {
            MGvar$sthreeDplot.regXcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$sthreeDplot.regXcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$sthreeDplot.regXcol.temp) > 0) 
                tkconfigure(StatColReg, bg = MGvar$sthreeDplot.regXcol.temp)
        }
        StatColReg <- tkbutton(TDPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$sthreeDplot.regXcol, command = function() ChangeStatColReg())
        tkplace(StatColReg, relx = 0.75, rely = 0.25, `in` = framestatlines)
        ChangeStatLines <- function() {
            StShowR = as.character(tclvalue(MGvar$Stat.RegLine.val))
            if (StShowR == "1") {
                MGvar$sthreeDplot.showregX <<- "yes"
                tclvalue(MGvar$Stat.RegLine.val) <<- 1
            }
            if (StShowR == "0") {
                MGvar$sthreeDplot.showregX <<- "no"
                tclvalue(MGvar$Stat.RegLine.val) <<- 0
            }
            MGvar$sthreeDplot.regXcol <<- MGvar$sthreeDplot.regXcol.temp
            tkrreplot(img3Dstat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeStatLines()), relx = 0.325, 
            rely = 0.85, `in` = framestatlines)
        Stat3DAxes <- tk2frame(StatNB)
        tkadd(StatNB, Stat3DAxes, text = "Axes")
        framestataxes <- tkwidget(Stat3DAxes, "TitleFrame", text = "Static 3D Configuration Axes", 
            background = "white")
        tkplace(framestataxes, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = Stat3DAxes)
        tkplace(tklabel(framestataxes, text = "Display Axes", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = framestataxes)
        Stat.AxesMeas.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.AxesMeas.CB, variable = MGvar$Stat.AxesMeas.val)
        tkplace(Stat.AxesMeas.CB, relx = 0.75, rely = 0.1, `in` = framestataxes)
        tkplace(tklabel(framestataxes, text = "Axes Colour", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = framestataxes)
        ChangeStatColAxes <- function() {
            MGvar$sthreeDplot.axescol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$sthreeDplot.axescol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$sthreeDplot.axescol.temp) > 0) 
                tkconfigure(StatColAxes, bg = MGvar$sthreeDplot.axescol.temp)
        }
        StatColAxes <- tkbutton(TDPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$sthreeDplot.axescol, command = function() ChangeStatColAxes())
        tkplace(StatColAxes, relx = 0.75, rely = 0.25, `in` = framestataxes)
        tkplace(tklabel(framestataxes, text = "Display Grid", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = framestataxes)
        Stat.Grid.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(Stat.Grid.CB, variable = MGvar$Stat.Grid.val)
        tkplace(Stat.Grid.CB, relx = 0.75, rely = 0.4, `in` = framestataxes)
        ChangeStatAxes <- function() {
            StShowA = as.character(tclvalue(MGvar$Stat.AxesMeas.val))
            if (StShowA == "0") {
                MGvar$sthreeDplot.showaxes <<- "no"
                tclvalue(MGvar$Stat.AxesMeas.val) <<- 0
            }
            if (StShowA == "1") {
                MGvar$sthreeDplot.showaxes <<- "yes"
                tclvalue(MGvar$Stat.AxesMeas.val) <<- 1
            }
            MGvar$sthreeDplot.axescol <<- MGvar$sthreeDplot.axescol.temp
            StGrid = as.character(tclvalue(MGvar$Stat.Grid.val))
            if (StGrid == "1") {
                MGvar$sthreeDplot.showgrid <<- "yes"
                tclvalue(MGvar$Stat.Grid.val) <<- 1
            }
            if (StGrid == "0") {
                MGvar$sthreeDplot.showgrid <<- "no"
                tclvalue(MGvar$Stat.Grid.val) <<- 0
            }
            tkrreplot(img3Dstat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeStatAxes()), relx = 0.325, 
            rely = 0.85, `in` = framestataxes)
        RGLNB <- tk2notebook(TDPlotOps, tabs = NULL)
        RGL3DGen <- tk2frame(RGLNB)
        tkadd(RGLNB, RGL3DGen, text = "General")
        frameRGLgen <- tkwidget(RGL3DGen, "TitleFrame", text = "RGL 3D Configuration General", 
            background = "white")
        tkplace(frameRGLgen, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = RGL3DGen)
        tkplace(tklabel(frameRGLgen, text = "Display Main Title", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameRGLgen)
        RGL.Main.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(RGL.Main.CB, variable = MGvar$RGL.Main.val)
        tkplace(RGL.Main.CB, relx = 0.75, rely = 0.1, `in` = frameRGLgen)
        tkplace(tklabel(frameRGLgen, text = "Display Axes Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameRGLgen)
        RGL.axlabs.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(RGL.axlabs.CB, variable = MGvar$RGL.axlabs.val)
        tkplace(RGL.axlabs.CB, relx = 0.75, rely = 0.25, `in` = frameRGLgen)
        ChangeRGLGen <- function() {
            RGLMain = as.character(tclvalue(MGvar$RGL.Main.val))
            if (RGLMain == "0") {
                MGvar$rglplot.title.show <<- "no"
                tclvalue(MGvar$RGL.Main.val) <<- 0
            }
            if (RGLMain == "1") {
                MGvar$rglplot.title.show <<- "yes"
                tclvalue(MGvar$RGL.Main.val) <<- 1
            }
            RGLLabs = as.character(tclvalue(MGvar$RGL.axlabs.val))
            if (RGLLabs == "1") {
                MGvar$rglplot.xlab <<- paste("Dim", MGvar$PlottingDimX.3D)
                MGvar$rglplot.ylab <<- paste("Dim", MGvar$PlottingDimY.3D)
                MGvar$rglplot.zlab <<- paste("Dim", MGvar$PlottingDimZ.3D)
                tclvalue(MGvar$RGL.axlabs.val) <<- 1
            }
            if (RGLLabs == "0") {
                MGvar$rglplot.xlab <<- ""
                MGvar$rglplot.ylab <<- ""
                MGvar$rglplot.zlab <<- ""
                tclvalue(MGvar$RGL.axlabs.val) <<- 0
            }
            plotting3D(MGvar$MDSmat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeRGLGen()), relx = 0.325, 
            rely = 0.85, `in` = frameRGLgen)
        RGL3DPoints <- tk2frame(RGLNB)
        tkadd(RGLNB, RGL3DPoints, text = "Points")
        frameRGLpoints <- tkwidget(RGL3DPoints, "TitleFrame", 
            text = "RGL 3D Configuration Points", background = "white")
        tkplace(frameRGLpoints, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = RGL3DPoints)
        tkplace(tklabel(frameRGLpoints, text = "Display Points", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameRGLpoints)
        RGL.points.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(RGL.points.CB, variable = MGvar$RGL.points.val)
        tkplace(RGL.points.CB, relx = 0.75, rely = 0.1, `in` = frameRGLpoints)
        tkplace(tklabel(frameRGLpoints, text = "Display Points Labels", 
            background = "white"), relx = 0.1, rely = 0.25, `in` = frameRGLpoints)
        RGL.ptlabels.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(RGL.ptlabels.CB, variable = MGvar$RGL.ptlabels.val)
        tkplace(RGL.ptlabels.CB, relx = 0.75, rely = 0.25, `in` = frameRGLpoints)
        tkplace(tklabel(frameRGLpoints, text = "Point Size", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameRGLpoints)
        RGL.PS.var <- tclVar(MGvar$rglplot.ptsize)
        RGL.PS.spin <- tk2spinbox(TDPlotOps, from = 0, to = 4, 
            increment = 0.5, width = 15)
        tkconfigure(RGL.PS.spin, textvariable = RGL.PS.var)
        tkplace(RGL.PS.spin, relx = 0.65, rely = 0.4, `in` = frameRGLpoints)
        tkplace(tklabel(frameRGLpoints, text = "Point Colour", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameRGLpoints)
        ChangeRGLptCol <- function() {
            MGvar$rglplot.ptcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = MGvar$rglplot.ptcol.temp, 
                  title = "Choose a colour"))))
            if (nchar(MGvar$rglplot.ptcol.temp) > 0) 
                tkconfigure(RGLptCOl, bg = MGvar$rglplot.ptcol.temp)
        }
        RGLptCOl <- tkbutton(TDPlotOps, text = "", width = 2, 
            height = 1, bg = MGvar$rglplot.ptcol, command = function() ChangeRGLptCol())
        tkplace(RGLptCOl, relx = 0.75, rely = 0.55, `in` = frameRGLpoints)
        ChangeRGLPoints <- function() {
            RGLpoints = as.character(tclvalue(MGvar$RGL.points.val))
            if (RGLpoints == "0") {
                MGvar$rglplot.showpoints <<- "no"
                tclvalue(MGvar$RGL.points.val) <<- 0
            }
            if (RGLpoints == "1") {
                MGvar$rglplot.showpoints <<- "yes"
                tclvalue(MGvar$RGL.points.val) <<- 1
            }
            RGLlabels = as.character(tclvalue(MGvar$RGL.ptlabels.val))
            if (RGLlabels == "0") {
                MGvar$rglplot.showlabels <<- "no"
                tclvalue(MGvar$RGL.ptlabels.val) <<- 0
            }
            if (RGLlabels == "1") {
                MGvar$rglplot.showlabels <<- "yes"
                tclvalue(MGvar$RGL.ptlabels.val) <<- 1
            }
            MGvar$rglplot.ptsize <<- as.numeric(tclvalue(RGL.PS.var))
            MGvar$rglplot.ptcol <<- MGvar$rglplot.ptcol.temp
            plotting3D(MGvar$MDSmat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeRGLPoints()), relx = 0.325, 
            rely = 0.85, `in` = frameRGLpoints)
        RGL3DAxes <- tk2frame(RGLNB)
        tkadd(RGLNB, RGL3DAxes, text = "Axes")
        frameRGLaxes <- tkwidget(RGL3DAxes, "TitleFrame", text = "RGL 3D Configuration Axes", 
            background = "white")
        tkplace(frameRGLaxes, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = RGL3DAxes)
        tkplace(tklabel(frameRGLaxes, text = "Display Axes Measures", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameRGLaxes)
        RGL.AxesMeas.CB <- tk2checkbutton(TDPlotOps)
        tkconfigure(RGL.AxesMeas.CB, variable = MGvar$RGL.AxesMeas.val)
        tkplace(RGL.AxesMeas.CB, relx = 0.75, rely = 0.1, `in` = frameRGLaxes)
        ChangeRGLAxes <- function() {
            RGLaxmeas = as.character(tclvalue(MGvar$RGL.AxesMeas.val))
            if (RGLaxmeas == "0") {
                MGvar$rglplot.axmeas <<- "no"
                tclvalue(MGvar$RGL.AxesMeas.val) <<- 0
            }
            if (RGLaxmeas == "1") {
                MGvar$rglplot.axmeas <<- "yes"
                tclvalue(MGvar$RGL.AxesMeas.val) <<- 1
            }
            plotting3D(MGvar$MDSmat)
        }
        tkplace(tkbutton(TDPlotOps, text = "Change", width = 15, 
            command = function() ChangeRGLAxes()), relx = 0.325, 
            rely = 0.85, `in` = frameRGLaxes)
        RGLOps <- tk2frame(TDmainNB)
        tkadd(TDmainNB, RGLOps, text = "RGL Plot")
        tkplace(StatNB, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = StatOps)
        tkplace(RGLNB, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = RGLOps)
        tkplace(TDmainNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.88, `in` = TDPlotOps)
        OnOK.3D <- function() {
            ChangeStatGen()
            ChangeStatPoints()
            ChangeStatLines()
            ChangeStatAxes()
            ChangeRGLGen()
            ChangeRGLPoints()
            ChangeRGLAxes()
            tkdestroy(TDPlotOps)
        }
        OnDefault.3D <- function() {
            tkdestroy(TDPlotOps)
        }
        tkplace(tkbutton(TDPlotOps, text = "OK", width = 15, 
            command = function() OnOK.3D()), relx = 0.15, rely = 0.92, 
            `in` = TDPlotOps)
        tkplace(tkbutton(TDPlotOps, text = "Default", width = 15, 
            command = function() OnDefault.3D()), relx = 0.55, 
            rely = 0.92, `in` = TDPlotOps)
        tkfocus(TDPlotOps)
        tkwait.window(TDPlotOps)
    }
    DataOptions <- function(Title = MGvar$datatitle) {
        DataOps = tktoplevel()
        tkwm.resizable(DataOps, "0", "0")
        tkwm.deiconify(DataOps)
        tkwm.title(DataOps, "Data Options")
        tkwm.geometry(DataOps, "420x450")
        Datacanvas = tkcanvas(DataOps, width = "1128", height = "756", 
            bg = col.sec)
        tkplace(Datacanvas, relx = 0, rely = 0, relwidth = 1, 
            relheight = 1, `in` = DataOps)
        DataOpsNB <- tk2notebook(DataOps, tabs = NULL)
        MGvar$activedata <- tk2frame(DataOpsNB)
        tkadd(DataOpsNB, MGvar$activedata, text = "Active Data")
        frameAD.T <- tkwidget(MGvar$activedata, "TitleFrame", 
            text = "Transpose", background = "white")
        tkplace(frameAD.T, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.38, `in` = MGvar$activedata)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frameAD.T, text = "All procedures in this package require that the active data have\nobjects as rows and variables as columns. If your data is not in\nthis format then please transpose.", 
            font = fontsmall), relx = 0.08, rely = 0.15, `in` = frameAD.T)
        tkplace(tklabel(frameAD.T, text = "Transpose Active Data?", 
            background = "white"), relx = 0.15, rely = 0.56, 
            `in` = frameAD.T)
        cbtrans <- tkcheckbutton(DataOps)
        cbTValue <- tclVar("0")
        tkconfigure(cbtrans, variable = cbTValue)
        tkplace(cbtrans, relx = 0.7, rely = 0.55, `in` = frameAD.T)
        ChangeT <- function() {
            TOp = as.character(tclvalue(cbTValue))
            if (TOp == "1") {
                MGvar$activedata <<- as.matrix(MGvar$activedata)
                if (nrow(MGvar$activedata) == 1 && ncol(MGvar$activedata) == 
                  1) {
                  tkmessageBox(message = "No data to transpose!", 
                    icon = "error")
                }
                else {
                  MGvar$activedata <<- t(MGvar$activedata)
                  MGvar$maxdims <<- nrow(MGvar$activedata) - 
                    1
                  MGvar$tShepx <<- as.vector(0)
                  tkmessageBox(message = paste("Your data has been transposed"))
                }
            }
        }
        tkplace(tkbutton(DataOps, text = "Change", width = 15, 
            command = function() ChangeT()), relx = 0.325, rely = 0.75, 
            `in` = frameAD.T)
        frameAD.Name <- tkwidget(MGvar$activedata, "TitleFrame", 
            text = "Change Name", background = "white")
        tkplace(frameAD.Name, relx = 0.02, relwidth = 0.96, rely = 0.41, 
            relheight = 0.3, `in` = MGvar$activedata)
        tclvalue(MGvar$datnam) <- paste("Current name of active data is:               ", 
            MGvar$datatitle)
        tkplace(tklabel(frameAD.Name, text = tclvalue(MGvar$datnam), 
            textvariable = MGvar$datnam, background = "white"), 
            relx = 0.1, rely = 0.2, `in` = frameAD.Name)
        tkplace(tklabel(frameAD.Name, text = "Enter new name of Dataset", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameAD.Name)
        ChangingName = tclVar("")
        NewNameBox = tkentry(DataOps, width = 15, textvariable = ChangingName)
        tkplace(NewNameBox, relx = 0.6, rely = 0.4, `in` = frameAD.Name)
        ChangeN <- function() {
            MGvar$activedata <<- as.matrix(MGvar$activedata)
            if (nrow(MGvar$activedata) == 1 && ncol(MGvar$activedata) == 
                1) {
                tkmessageBox(message = "No data to name!", icon = "error")
            }
            else {
                NName = as.character(tclvalue(ChangingName))
                MGvar$datatitle <<- NName
                tclvalue(MGvar$datnam) <<- paste("Current name of active data is:               ", 
                  MGvar$datatitle)
                tclvalue(labelText) <<- paste("Active Dataset is", 
                  MGvar$datatitle)
            }
        }
        tkplace(tkbutton(DataOps, text = "Change Name", width = 15, 
            command = function() ChangeN()), relx = 0.325, rely = 0.7, 
            `in` = frameAD.Name)
        frameAD.Scale <- tkwidget(MGvar$activedata, "TitleFrame", 
            text = "Scale Data", background = "white")
        tkplace(frameAD.Scale, relx = 0.02, relwidth = 0.96, 
            rely = 0.72, relheight = 0.27, `in` = MGvar$activedata)
        tkplace(tklabel(frameAD.Scale, text = "Scaling Data will scale the columns of your data between 0 and 1.", 
            font = fontsmall), relx = 0.07, rely = 0.15, `in` = frameAD.Scale)
        tkplace(tklabel(frameAD.Scale, text = "Scale your active data?", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = frameAD.Scale)
        ScDat.val <- tclVar(0)
        ScDat.CB <- tk2checkbutton(DataOps)
        tkconfigure(ScDat.CB, variable = ScDat.val)
        tkplace(ScDat.CB, relx = 0.75, rely = 0.4, `in` = frameAD.Scale)
        ChangeS <- function() {
            Scl <- as.character(tclvalue(ScDat.val))
            if (Scl == "1") {
                MGvar$activedata <<- stand.1.range(MGvar$activedata)
            }
        }
        tkplace(tkbutton(DataOps, text = "Change", width = 15, 
            command = function() ChangeS()), relx = 0.325, rely = 0.65, 
            `in` = frameAD.Scale)
        tkplace(DataOpsNB, relx = 0.05, relwidth = 0.9, rely = 0.01, 
            relheight = 0.85, `in` = DataOps)
        OnOK.Dat <- function() {
            tkdestroy(DataOps)
        }
        OnCancel.Dat <- function() {
            tkdestroy(DataOps)
        }
        tkplace(tkbutton(DataOps, text = "OK", width = 15, command = function() OnOK.Dat()), 
            relx = 0.15, rely = 0.9, `in` = DataOps)
        tkplace(tkbutton(DataOps, text = "Cancel", width = 15, 
            command = function() OnCancel.Dat()), relx = 0.55, 
            rely = 0.9, `in` = DataOps)
        tkfocus(DataOps)
        tkbind(DataOps, "<Return>", OnOK.Dat)
        tkwait.window(DataOps)
    }
    Appearance.Change <- function() {
        Appearance = tktoplevel()
        tkwm.resizable(Appearance, "0", "0")
        tkwm.deiconify(Appearance)
        tkwm.title(Appearance, "Appearance Options")
        tkwm.geometry(Appearance, "320x150")
        AppearanceCanvas = tkcanvas(Appearance, bg = col.sec)
        tkplace(AppearanceCanvas, relx = 0, rely = 0, relheight = 1, 
            relwidth = 1, `in` = Appearance)
        frameApp <- tkwidget(Appearance, "TitleFrame", text = "Appearance Changes", 
            background = "white")
        tkplace(frameApp, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = Appearance)
        tkplace(tklabel(frameApp, text = "Change Primary Background Colour", 
            background = "white"), relx = 0.05, rely = 0.15, 
            `in` = frameApp)
        ChangeColPrim <- function() {
            col.prim.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = col.prim, title = "Choose a colour"))))
            if (nchar(col.prim.temp) > 0) 
                tkconfigure(ColPrim, bg = col.prim.temp)
        }
        ColPrim <- tkbutton(Appearance, text = "", width = 2, 
            height = 1, bg = col.prim, command = function() ChangeColPrim())
        tkplace(ColPrim, relx = 0.85, rely = 0.15, `in` = frameApp)
        tkplace(tklabel(frameApp, text = "Change Secondary Background Colour", 
            background = "white"), relx = 0.05, rely = 0.45, 
            `in` = frameApp)
        ChangeColSec <- function() {
            col.sec.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = col.sec, title = "Choose a colour"))))
            if (nchar(col.sec.temp) > 0) 
                tkconfigure(ColSec, bg = col.sec.temp)
        }
        ColSec <- tkbutton(Appearance, text = "", width = 2, 
            height = 1, bg = col.sec, command = function() ChangeColSec())
        tkplace(ColSec, relx = 0.85, rely = 0.45, `in` = frameApp)
        OnOK.App <- function() {
            surechange = tkmessageBox(message = "Are you sure you want to make these changes?", 
                type = "yesno", default = "no")
            if (as.character(surechange) == "yes") {
                tkdestroy(Appearance)
                col.prim <<- col.prim.temp
                tkplace(tklabel(mytt, text = "MDS-GUI under Construction", 
                  foreground = "black", font = fontHeading, background = col.prim), 
                  relx = 0.52, rely = 0.01, `in` = mytt)
                tkplace(tklabel(mytt, text = tclvalue(MGvar$UserName), 
                  textvariable = MGvar$UserName, font = UNHeading.Font, 
                  background = col.prim), relx = 0.01, rely = 0.005, 
                  `in` = mytt)
                col.sec <<- col.sec.temp
                tkconfigure(backcanvas, bg = col.prim)
                tkconfigure(AppearanceCanvas, bg = col.sec)
            }
            else {
                tkdestroy(Appearance)
            }
        }
        Default <- function() {
            col.prim <<- "#dae0f1"
            col.sec <<- "#dae0f1"
            tkconfigure(backcanvas, bg = col.prim)
            tkconfigure(AppearanceCanvas, bg = col.sec)
            tkplace(tklabel(mytt, text = "MDS-GUI under Construction", 
                foreground = "black", font = fontHeading, background = col.prim), 
                relx = 0.52, rely = 0.01, `in` = mytt)
            tkplace(tklabel(mytt, text = tclvalue(MGvar$UserName), 
                textvariable = MGvar$UserName, font = UNHeading.Font, 
                background = col.prim), relx = 0.01, rely = 0.005, 
                `in` = mytt)
            tkdestroy(Appearance)
        }
        tkplace(tkbutton(Appearance, text = "Change", width = 10, 
            command = function() OnOK.App()), relx = 0.2, rely = 0.75, 
            `in` = frameApp)
        tkplace(tkbutton(Appearance, text = "Default", width = 10, 
            command = function() Default()), relx = 0.55, rely = 0.75, 
            `in` = frameApp)
    }
    plotting2DZoom <- function(data, showtitle = MGvar$zoomedplot.title.show, 
        showmeas = MGvar$zoomedplot.distmeas, xlabel = MGvar$zoomedplot.xlab, 
        ylabel = MGvar$zoomedplot.ylab, bgcol = MGvar$zoomedplot.bg, 
        pointcex = MGvar$zoomedplot.cex, showlabs = MGvar$zoomedplot.labs, 
        showpoints = MGvar$zoomedplot.showpoints, pointcol = MGvar$zoomedplot.pointcol, 
        pointshape = MGvar$zoomedplot.type, ymeas = MGvar$zoomedplot.yaxt, 
        xmeas = MGvar$zoomedplot.yaxt, axcol = MGvar$zoomedplot.axescol, 
        showzoomX = MGvar$zoomedplot.showzoom, showreg = MGvar$zoomedplot.showreg, 
        regcol = MGvar$zoomedplot.regcol, showleg = MGvar$zoomedplot.showleg) {
        if (showtitle == "yes") {
            graphtitle = MGvar$activeplot.title
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showmeas == "yes") {
            distanceMeasure = MGvar$dMeas
        }
        if (showmeas == "no") {
            distanceMeasure = ""
        }
        if (showpoints == "yes") {
            pointtype = "p"
        }
        if (showpoints == "no") {
            pointtype = "n"
        }
        params <- par(bg = bgcol)
        par(mar = c(3, 3, 3, 3))
        par(cex.axis = 0.8)
        Xmin <- min(data[, 1])
        Xmax <- max(data[, 1])
        Ymin <- min(data[, 2])
        Ymax <- max(data[, 2])
        Xrange <- Xmax - Xmin
        Yrange <- Ymax - Ymin
        if (showleg == "yes") {
            Ymax = Ymax + 0.03 * Yrange
            Xmax = Xmax + 0.015 * Xrange
            Xmin = Xmin - 0.015 * Xrange
        }
        plot(data, type = pointtype, ylim = c(Ymin, Ymax), xlim = c(Xmin, 
            Xmax), asp = 1, sub = distanceMeasure, xaxt = xmeas, 
            yaxt = ymeas, main = graphtitle, bg = bgcol, cex = pointcex, 
            col = pointcol, fg = axcol, pch = pointshape, ylab = "")
        par(params)
        if (showreg == "yes") {
            regline = lm(data[, 2] ~ data[, 1])
            abline(regline, col = regcol)
        }
        if (showleg == "yes") {
            legend("topleft", "MDS Configuration Points", pch = pointshape, 
                bty = "n", cex = 0.8, pt.bg = "white", col = pointcol)
            if (showreg == "yes") {
                leglab = paste("Regression Line")
                legend("topright", 1.09 * Ymax, leglab, bty = "n", 
                  cex = 0.8, pt.bg = "white", lty = 1, col = regcol)
            }
        }
        mtext(text = xlabel, side = 1, line = 1.9, adj = 0.5, 
            cex = 1, col = "black")
        mtext(text = ylabel, side = 2, line = 1.9, adj = 0.5, 
            cex = 1, col = "black")
        if (nrow(data) == 1 && ncol(data) == 1) {
        }
        else {
            if (showlabs == "yes") {
                text(data, rownames(data), cex = pointcex, col = pointcol)
            }
            if (showlabs == "no") {
            }
        }
        if (showmeas == "yes") {
            mtext(text = MGvar$dMeas, side = 3, line = 2.4, adj = 0.01, 
                cex = 0.7, col = "black")
        }
        MGvar$parPlotSizeZ <<- par("plt")
        MGvar$usrCoordsZ <<- par("usr")
        MGvar$ZlabelsVec <<- rownames(data)
        if (length(MGvar$indexZLabeled) > 0) {
            for (i in (1:length(MGvar$indexZLabeled))) {
                ZClosest <- MGvar$indexZLabeled[i]
                text(data[ZClosest, 1], data[ZClosest, 2], labels = MGvar$ZlabelsVec[ZClosest], 
                  pos = 3)
            }
        }
        ZoomAmount <- tclVar("")
        if (showzoomX == "yes") {
            tclvalue(ZoomAmount) <- paste("Zoom: ", tclvalue(MGvar$Zoom.Rat.var), 
                "X")
            mtext(text = tclvalue(ZoomAmount), side = 3, line = 2.4, 
                adj = 1, cex = 0.7, col = "black")
        }
    }
    labelClosestPointZoom <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        indexZClosest <- which.min(squared.Distance)
        MGvar$indexZLabeled <<- c(MGvar$indexZLabeled, indexZClosest)
        tkrreplot(MGcomp$zooming)
    }
    OnZoomLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        width <- as.numeric(tclvalue(tkwinfo("reqwidth", MGcomp$zooming)))
        height <- as.numeric(tclvalue(tkwinfo("reqheight", MGcomp$zooming)))
        xMin <- MGvar$parPlotSizeZ[1] * width
        xMax <- MGvar$parPlotSizeZ[2] * width
        yMin <- MGvar$parPlotSizeZ[3] * height
        yMax <- MGvar$parPlotSizeZ[4] * height
        rangeX <- MGvar$usrCoordsZ[2] - MGvar$usrCoordsZ[1]
        rangeY <- MGvar$usrCoordsZ[4] - MGvar$usrCoordsZ[3]
        imgXcoords <- (MGvar$newCoords[, 1] - MGvar$usrCoordsZ[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$newCoords[, 2] - MGvar$usrCoordsZ[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoordsZ[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoordsZ[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        msg <- paste("Label the point closest to these approximate plot coordinates: \n", 
            "x =", format(xPlotCoord, digits = 2), ",y =", format(yPlotCoord, 
                digits = 2), "?")
        mbval <- tkmessageBox(title = "Label Point Closest to These Approximate Plot Coordinates", 
            message = msg, type = "yesno", icon = "question")
        if (tclvalue(mbval) == "yes") {
            labelClosestPointZoom(xClick, yClick, imgXcoords, 
                imgYcoords)
        }
    }
    GetClosestPointZoom <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        indexZMClosest <- which.min(squared.Distance)
        MGvar$activeX <<- MGvar$MDSmat[indexZMClosest, 1]
        MGvar$activeY <<- MGvar$MDSmat[indexZMClosest, 2]
    }
    GetCoordsLeftClick <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$MDSmat[, 1] - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$MDSmat[, 2] - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        msg <- paste("Zoom with focal point at: \n", "x =", format(xPlotCoord, 
            digits = 2), ",y =", format(yPlotCoord, digits = 2), 
            "?")
        mbval <- tkmessageBox(title = "Zoom cursor selection", 
            message = msg, type = "yesno", icon = "question")
        if (tclvalue(mbval) == "yes") {
            GetClosestPointZoom(xClick, yClick, imgXcoords, imgYcoords)
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                tkbind(img, "<Button-1>", OnPlotLeftClick)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                tkbind(img2, "<Button-1>", OnPlotLeftClick)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                tkbind(img3, "<Button-1>", OnPlotLeftClick)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                tkbind(img4, "<Button-1>", OnPlotLeftClick)
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                tkbind(img5, "<Button-1>", OnPlotLeftClick)
            }
            tkdestroy(MGcomp$orp)
        }
    }
    seczoomplot <- function(data, showtitle = MGvar$zoomedplot.title.show, 
        showmeas = MGvar$zoomedplot.distmeas, bgcol = MGvar$zoomedplot.bg, 
        pointcex = MGvar$zoomedplot.cex, showlabs = MGvar$zoomedplot.labs, 
        showpoints = MGvar$zoomedplot.showpoints, pointcol = MGvar$zoomedplot.pointcol, 
        pointshape = MGvar$zoomedplot.type, axcol = MGvar$zoomedplot.axescol, 
        showzoomX = MGvar$zoomedplot.showzoom, Tab = tclvalue(MGvar$ActivePlottingTab), 
        showreg = MGvar$zoomedplot.showreg, regcol = MGvar$zoomedplot.regcol, 
        showleg = MGvar$zoomedplot.showleg) {
        if (showtitle == "yes") {
            graphtitle = MGvar$activeplot.title
        }
        if (showtitle == "no") {
            graphtitle = ""
        }
        if (showmeas == "yes") {
            distanceMeasure = MGvar$dMeas
        }
        if (showmeas == "no") {
            distanceMeasure = ""
        }
        if (showpoints == "yes") {
            pointtype = "p"
        }
        if (showpoints == "no") {
            pointtype = "n"
        }
        params <- par(bg = bgcol)
        par(mar = c(1.5, 3, 1.5, 3))
        par(cex.axis = 0.8)
        if (nrow(data) == 1 && ncol(data) == 1) {
            plot(data, type = "n", xaxt = "n", yaxt = "n", ylab = "n", 
                xlab = "n")
            par(params)
        }
        else {
            Xmin <- min(data[, 1])
            Xmax <- max(data[, 1])
            Ymin <- min(data[, 2])
            Ymax <- max(data[, 2])
            Xrange <- Xmax - Xmin
            Yrange <- Ymax - Ymin
            if (showleg == "yes") {
                Ymax = Ymax + 0.03 * Yrange
                Xmax = Xmax + 0.015 * Xrange
                Xmin = Xmin - 0.015 * Xrange
            }
            plot(data, type = pointtype, ylim = c(Ymin, Ymax), 
                xlim = c(Xmin, Xmax), asp = 1, bg = bgcol, cex = pointcex, 
                xaxt = "n", yaxt = "n", ylab = "n", xlab = "n", 
                col = pointcol, fg = axcol, pch = pointshape)
            if (nrow(data) == 1 && ncol(data) == 1) {
            }
            else {
                if (showlabs == "yes") {
                  text(data, rownames(data), cex = pointcex, 
                    col = pointcol)
                }
                if (showlabs == "no") {
                }
            }
            if (showleg == "yes") {
                legend("topleft", "MDS Configuration Points", 
                  pch = pointshape, bty = "n", cex = 0.8, pt.bg = "white", 
                  col = pointcol)
                if (showreg == "yes") {
                  leglab = paste("Regression Line")
                  legend("topright", 1.09 * Ymax, leglab, bty = "n", 
                    cex = 0.8, pt.bg = "white", lty = 1, col = regcol)
                }
            }
            mtext(text = graphtitle, side = 2, line = 1.5, adj = 0.5, 
                cex = 1, col = "black")
            mtext(text = distanceMeasure, side = 2, line = 2.4, 
                adj = 0, cex = 0.7, col = "black")
            if (showreg == "yes") {
                regline = lm(data[, 2] ~ data[, 1])
                abline(regline, col = regcol)
            }
            ZoomAmount <- tclVar("")
            if (showzoomX == "yes") {
                tclvalue(ZoomAmount) <- paste("Zoom: ", tclvalue(MGvar$Zoom.Rat.var), 
                  "X")
                if (MGvar$seczoomswitch == "on") {
                  mtext(text = tclvalue(ZoomAmount), side = 2, 
                    line = 2.4, adj = 1, cex = 0.7, col = "black")
                }
            }
            ZoomedTab <- tclVar("")
            tclvalue(ZoomedTab) <- paste("Zoomed image of ", 
                Tab)
            if (MGvar$seczoomswitch == "on") {
                mtext(text = tclvalue(ZoomedTab), side = 4, line = 1.5, 
                  adj = 0.5, cex = 0.9, col = "black")
            }
            MGvar$parPlotSizeZ <<- par("plt")
            MGvar$usrCoordsZ <<- par("usr")
            MGvar$ZlabelsVec <<- rownames(data)
            if (length(MGvar$indexZLabeled) > 0) {
                for (i in (1:length(MGvar$indexZLabeled))) {
                  ZClosest <- MGvar$indexZLabeled[i]
                  text(data[ZClosest, 1], data[ZClosest, 2], 
                    labels = MGvar$ZlabelsVec[ZClosest], pos = 3)
                }
            }
        }
    }
    ZoomPlot <- function() {
        Zoomtt = tktoplevel()
        tkwm.resizable(Zoomtt, "0", "0")
        tkwm.deiconify(Zoomtt)
        tkwm.title(Zoomtt, "Advanced Zooming on Active Plot")
        tkwm.geometry(Zoomtt, "350x400")
        Zoomcanvas = tkcanvas(Zoomtt, width = 420, height = 450, 
            bg = col.sec)
        tkplace(Zoomcanvas, relx = 0, rely = 0, `in` = Zoomtt)
        frameZ1 <- tkwidget(Zoomtt, "TitleFrame", text = "Zooming Amount", 
            background = "white")
        tkplace(frameZ1, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.2, `in` = Zoomtt)
        tkplace(tklabel(frameZ1, text = "Zoom Ratio", background = "white"), 
            relx = 0.04, rely = 0.4, `in` = frameZ1)
        MGvar$Zoom.Rat.var <<- tclVar("1")
        Zoom.Rat.spin <- tk2spinbox(Zoomtt, from = 1, to = 20, 
            increment = 0.5, width = 10)
        tkconfigure(Zoom.Rat.spin, textvariable = MGvar$Zoom.Rat.var)
        tkplace(Zoom.Rat.spin, relx = 0.7, rely = 0.4, `in` = frameZ1)
        ZrbValue <- tclVar("Centre")
        frameZ2 <- tkwidget(Zoomtt, "TitleFrame", text = "Focal Point", 
            background = "white")
        tkplace(frameZ2, relx = 0.05, rely = 0.23, relwidth = 0.9, 
            relheight = 0.3, `in` = Zoomtt)
        tkplace(tklabel(frameZ2, text = "Focal Point:", background = "white"), 
            relx = 0.04, rely = 0.2, `in` = frameZ2)
        tkplace(tklabel(frameZ2, text = "Centre", background = "white"), 
            relx = 0.35, rely = 0.2, `in` = frameZ2)
        rb.centre <- tk2radiobutton(Zoomtt)
        tkconfigure(rb.centre, variable = ZrbValue, value = "Centre")
        tkplace(rb.centre, relx = 0.9, rely = 0.2, `in` = frameZ2)
        tkplace(tklabel(frameZ2, text = "Cursor Location", background = "white"), 
            relx = 0.35, rely = 0.45, `in` = frameZ2)
        rb.loc <- tk2radiobutton(Zoomtt)
        tkconfigure(rb.loc, variable = ZrbValue, value = "Location")
        tkplace(rb.loc, relx = 0.9, rely = 0.45, `in` = frameZ2)
        tkplace(tklabel(frameZ2, text = "Point Name", background = "white"), 
            relx = 0.35, rely = 0.7, `in` = frameZ2)
        objnames = rownames(MGvar$activedata)
        PointEntry = tclVar(objnames[1])
        Pointtext.ComboBox <- tkwidget(Zoomtt, "ComboBox", editable = FALSE, 
            values = objnames, width = 10, textvariable = PointEntry)
        tkplace(Pointtext.ComboBox, relx = 0.6, rely = 0.7, `in` = frameZ2)
        rb.point <- tk2radiobutton(Zoomtt)
        tkconfigure(rb.point, variable = ZrbValue, value = "Point")
        tkplace(rb.point, relx = 0.9, rely = 0.7, `in` = frameZ2)
        frameZ3 <- tkwidget(Zoomtt, "TitleFrame", text = "Focal Point", 
            background = "white")
        tkplace(frameZ3, relx = 0.05, rely = 0.55, relwidth = 0.9, 
            relheight = 0.3, `in` = Zoomtt)
        tkplace(tklabel(frameZ3, text = "Output Location:", background = "white"), 
            relx = 0.04, rely = 0.2, `in` = frameZ3)
        tkplace(tklabel(frameZ3, text = "Replace Active Plot", 
            background = "white"), relx = 0.35, rely = 0.2, `in` = frameZ3)
        Zoom.AP.var <- tclVar("0")
        Zoom.AP.CB <- tk2checkbutton(Zoomtt)
        tkconfigure(Zoom.AP.CB, variable = Zoom.AP.var)
        tkplace(Zoom.AP.CB, relx = 0.85, rely = 0.2, `in` = frameZ3)
        tkplace(tklabel(frameZ3, text = "Secondary Tabbed-Book", 
            background = "white"), relx = 0.35, rely = 0.45, 
            `in` = frameZ3)
        Zoom.SP.var <- tclVar("1")
        Zoom.SP.CB <- tk2checkbutton(Zoomtt)
        tkconfigure(Zoom.SP.CB, variable = Zoom.SP.var)
        tkplace(Zoom.SP.CB, relx = 0.85, rely = 0.45, `in` = frameZ3)
        tkplace(tklabel(frameZ3, text = "Popped Out Window", 
            background = "white"), relx = 0.35, rely = 0.7, `in` = frameZ3)
        Zoom.POW.var <- tclVar("0")
        Zoom.POW.CB <- tk2checkbutton(Zoomtt)
        tkconfigure(Zoom.POW.CB, variable = Zoom.POW.var)
        tkplace(Zoom.POW.CB, relx = 0.85, rely = 0.7, `in` = frameZ3)
        Zoombut <- function() {
            surechange = tkmessageBox(message = "Are you sure you want to make these changes?", 
                type = "yesno", default = "no")
            if (as.character(surechange) == "yes") {
                ZoomRat = as.numeric(tclvalue(MGvar$Zoom.Rat.var))
                focpoint = as.character(tclvalue(ZrbValue))
                Pname = as.character(tclvalue(PointEntry))
                RepPlot = as.character(tclvalue(Zoom.AP.var))
                SecPlot = as.character(tclvalue(Zoom.SP.var))
                PopPlot = as.character(tclvalue(Zoom.POW.var))
                MGvar$indexZLabeled <<- c()
                MGvar$labeledZPoints <<- list()
                zoomfunc(MGvar$MDSmat, focpoint, Pname, ZoomRat, 
                  RepPlot, SecPlot, PopPlot)
            }
            else {
                tkdestroy(Zoomtt)
            }
        }
        zoomfunc <- function(Coords, focalpoint, pointname, ratio, 
            rep, sec, pop) {
            vertrange = max(Coords[, 2]) - min(Coords[, 2])
            horrange = max(Coords[, 1]) - min(Coords[, 1])
            zoom = ratio
            newrangearea = 1/zoom
            newvertrange = vertrange * newrangearea
            newhorrange = horrange * newrangearea
            datlength = nrow(Coords)
            if (focalpoint == "Centre") {
                MGvar$activeX <<- max(Coords[, 1]) - horrange/2
                MGvar$activeY <<- max(Coords[, 2]) - vertrange/2
            }
            if (focalpoint == "Location") {
                tkdestroy(Zoomtt)
                MGcomp$orp <<- tktoplevel()
                tkwm.resizable(MGcomp$orp, "0", "0")
                tkwm.title(MGcomp$orp, "Cursor Zoom")
                tkwm.geometry(MGcomp$orp, "300x90")
                orpcanvas = tkcanvas(MGcomp$orp, width = 300, 
                  height = 120, bg = col.sec)
                tkplace(orpcanvas, relx = 0, rely = 0, `in` = MGcomp$orp)
                frameOrp <- tkwidget(MGcomp$orp, "TitleFrame", 
                  text = "Active Cursor Zoom", background = "white")
                tkplace(frameOrp, relx = 0.05, rely = 0.05, relheight = 0.9, 
                  relwidth = 0.9, `in` = MGcomp$orp)
                tkplace(tklabel(frameOrp, text = "The cursor zoom function is active, please \nleft click a point on the active plot to perform\n the zoom.", 
                  background = "white"), relx = 0.05, rely = 0.2, 
                  `in` = frameOrp)
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  tkbind(img, "<Button-1>", GetCoordsLeftClick)
                  tkconfigure(img, cursor = "crosshair")
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  tkbind(img2, "<Button-1>", GetCoordsLeftClick)
                  tkconfigure(img2, cursor = "crosshair")
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  tkbind(img3, "<Button-1>", GetCoordsLeftClick)
                  tkconfigure(img3, cursor = "crosshair")
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  tkbind(img4, "<Button-1>", GetCoordsLeftClick)
                  tkconfigure(img4, cursor = "crosshair")
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  tkbind(img5, "<Button-1>", GetCoordsLeftClick)
                  tkconfigure(img5, cursor = "crosshair")
                }
                tkwait.window(MGcomp$orp)
                ActiveArrowCursor()
            }
            if (focalpoint == "Point") {
                for (i in 1:datlength) if (rownames(Coords)[i] == 
                  pointname) {
                  activepoint = i
                }
                MGvar$activeX <<- Coords[activepoint, 1]
                MGvar$activeY <<- Coords[activepoint, 2]
            }
            newX = as.vector(0)
            newY = as.vector(0)
            names = as.vector(0)
            count = 1
            for (i in 1:datlength) {
                if (abs(MGvar$activeX - Coords[i, 1]) < 0.5 * 
                  newhorrange && abs(MGvar$activeY - Coords[i, 
                  2]) < 0.5 * newvertrange) {
                  newX[count] = Coords[i, 1]
                  newY[count] = Coords[i, 2]
                  names[count] = rownames(Coords)[i]
                  count = count + 1
                }
            }
            MGvar$newCoords <<- as.matrix(0)
            MGvar$newCoords <<- matrix(nrow = length(newX), ncol = 2)
            MGvar$newCoords[, 1] <<- newX
            MGvar$newCoords[, 2] <<- newY
            rownames(MGvar$newCoords) <<- names
            count = 1
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                preindexLabeled = MGvar$indexLabeled.T1
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                preindexLabeled = MGvar$indexLabeled.T2
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                preindexLabeled = MGvar$indexLabeled.T3
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                preindexLabeled = MGvar$indexLabeled.T4
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                preindexLabeled = MGvar$indexLabeled.T5
            }
            if (length(preindexLabeled) > 0) {
                for (i in 1:length(rownames(MGvar$newCoords))) {
                  for (j in 1:length(preindexLabeled)) {
                    if (rownames(MGvar$newCoords)[i] == names(preindexLabeled)[j]) {
                      MGvar$indexZLabeled[count] <<- i
                      names(MGvar$indexZLabeled)[count] <<- names(preindexLabeled)[j]
                      count = count + 1
                    }
                  }
                }
            }
            if (rep == "1") {
                MGvar$MDSmat <<- MGvar$newCoords
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$activeplot.cex.T1 <<- 0.8
                  MGvar$Tab1.zoomedswitch <<- "on"
                  MGvar$newCoords.T1 <<- MGvar$newCoords
                  MGvar$indexLabeled.T1 <<- MGvar$indexZLabeled
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$activeplot.cex.T2 <<- 0.8
                  MGvar$Tab2.zoomedswitch <<- "on"
                  MGvar$newCoords.T2 <<- MGvar$newCoords
                  MGvar$indexLabeled.T2 <<- MGvar$indexZLabeled
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$activeplot.cex.T3 <<- 0.8
                  MGvar$Tab3.zoomedswitch <<- "on"
                  MGvar$newCoords.T3 <<- MGvar$newCoords
                  MGvar$indexLabeled.T3 <<- MGvar$indexZLabeled
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$activeplot.cex.T4 <<- 0.8
                  MGvar$Tab4.zoomedswitch <<- "on"
                  MGvar$newCoords.T4 <<- MGvar$newCoords
                  MGvar$indexLabeled.T4 <<- MGvar$indexZLabeled
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$activeplot.cex.T5 <<- 0.8
                  MGvar$Tab5.zoomedswitch <<- "on"
                  MGvar$newCoords.T5 <<- MGvar$newCoords
                  MGvar$indexLabeled.T5 <<- MGvar$indexZLabeled
                }
                tabplot()
            }
            if (sec == "1") {
                MGvar$seczoomswitch <<- "on"
                tk2notetab.select(mySecondaryNB, "Zoomed Plot")
                MGvar$zoomedplot.cex <<- 0.7
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$newCoords.T1 <<- MGvar$newCoords
                  MGvar$indexLabeled.T1 <<- MGvar$indexZLabeled
                  tkrreplot(imgseczoom, function() seczoomplot(MGvar$newCoords, 
                    Tab = "Tab1"))
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$newCoords.T2 <<- MGvar$newCoords
                  MGvar$indexLabeled.T2 <<- MGvar$indexZLabeled
                  tkrreplot(imgseczoom, function() seczoomplot(MGvar$newCoords, 
                    Tab = "Tab2"))
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$newCoords.T3 <<- MGvar$newCoords
                  MGvar$indexLabeled.T3 <<- MGvar$indexZLabeled
                  tkrreplot(imgseczoom, function() seczoomplot(MGvar$newCoords, 
                    Tab = "Tab3"))
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$newCoords.T4 <<- MGvar$newCoords
                  MGvar$indexLabeled.T4 <<- MGvar$indexZLabeled
                  tkrreplot(imgseczoom, function() seczoomplot(MGvar$newCoords, 
                    Tab = "Tab4"))
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$newCoords.T5 <<- MGvar$newCoords
                  MGvar$indexLabeled.T5 <<- MGvar$indexZLabeled
                  tkrreplot(imgseczoom, function() seczoomplot(MGvar$newCoords, 
                    Tab = "Tab5"))
                }
            }
            if (pop == "1") {
                ReusePOW()
                popzoomswitch <<- "on"
                zoomhscale <- 1.5
                zoomvscale <- 1.5
                MGcomp$zoomplottt <<- tktoplevel()
                tkwm.geometry(MGcomp$zoomplottt, "650x650")
                zoomplotcanvas = tkcanvas(MGcomp$zoomplottt, 
                  width = 650, height = 650, bg = col.sec)
                tkplace(zoomplotcanvas, relx = 0, rely = 0, `in` = MGcomp$zoomplottt)
                tkwm.title(MGcomp$zoomplottt, "Zoomed Coordinates")
                MGcomp$zooming <<- tkrplot(MGcomp$zoomplottt, 
                  function() plotting2DZoom(MGvar$newCoords), 
                  hscale = zoomhscale, vscale = zoomvscale)
                MGcomp$ZoomPlotMenu <<- tkmenu(MGcomp$zooming, 
                  tearoff = FALSE)
                tkadd(MGcomp$ZoomPlotMenu, "command", label = "Clear Added Point Labels", 
                  command = ClearZoomPoints)
                tkadd(MGcomp$ZoomPlotMenu, "command", label = "Copy Zoomed Plot to Clipboard", 
                  command = function() {
                    tkrreplot(MGcomp$zooming)
                  })
                tkbind(MGcomp$zooming, "<Button-3>", RightClickZoom)
                CopyToClip <- function() {
                  tkrreplot(MGcomp$zooming)
                }
                copy.but <- tkbutton(MGcomp$zoomplottt, text = "Copy to Clipboard", 
                  command = CopyToClip)
                tkplace(MGcomp$zooming, relx = 0.05, rely = 0.02, 
                  `in` = MGcomp$zoomplottt)
                tkplace(tkbutton(MGcomp$zoomplottt, text = "Copy to Clipboard", 
                  width = 20, command = function() CopyToClip()), 
                  relx = 0.2, rely = 0.95, `in` = MGcomp$zoomplottt)
                tkplace(tkbutton(MGcomp$zoomplottt, text = "Plot Options", 
                  width = 20, command = function() ZoomPlotOps()), 
                  relx = 0.6, rely = 0.95, `in` = MGcomp$zoomplottt)
                tkconfigure(MGcomp$zooming, cursor = "hand2")
                tkbind(MGcomp$zooming, "<Button-1>", OnZoomLeftClick)
                tkbind(MGcomp$zooming, "<Destroy>", function() {
                  popzoomswitch <<- "off"
                })
            }
            tkdestroy(Zoomtt)
        }
        OnCancel.Zoom <- function() {
            tkdestroy(Zoomtt)
        }
        tkplace(tkbutton(Zoomtt, text = "Zoom", width = 15, command = function() Zoombut()), 
            relx = 0.15, rely = 0.9, `in` = Zoomtt)
        tkplace(tkbutton(Zoomtt, text = "Cancel", width = 15, 
            command = function() OnCancel.Zoom()), relx = 0.55, 
            rely = 0.9, `in` = Zoomtt)
    }
    Clear <- function() {
        ClearOption = tkmessageBox(message = "Are you sure you want to clear all data and plots?", 
            type = "yesno", default = "yes")
        if (as.character(ClearOption) == "yes") {
            tkconfigure(mytt, cursor = "watch")
            MGvar$distmat <<- as.matrix(0)
            MGvar$activedata <<- as.matrix(0)
            MGvar$MDSmat <<- as.matrix(0)
            MGvar$distdat <<- as.matrix(0)
            MGvar$distmat.T1 <<- as.matrix(0)
            MGvar$distmat.T2 <<- as.matrix(0)
            MGvar$distmat.T3 <<- as.matrix(0)
            MGvar$distmat.T4 <<- as.matrix(0)
            MGvar$distmat.T5 <<- as.matrix(0)
            MGvar$MDSmat <<- as.matrix(0)
            MGvar$MDSmat.T1 <<- as.matrix(0)
            MGvar$MDSmat.T2 <<- as.matrix(0)
            MGvar$MDSmat.T3 <<- as.matrix(0)
            MGvar$MDSmat.T4 <<- as.matrix(0)
            MGvar$MDSmat.T5 <<- as.matrix(0)
            MGvar$MDS.dimensions <<- 2
            MGvar$Opt.dim <<- 1
            MGvar$Opt.dim.T1 <<- 1
            MGvar$Opt.dim.T2 <<- 1
            MGvar$Opt.dim.T3 <<- 1
            MGvar$Opt.dim.T4 <<- 1
            MGvar$Opt.dim.T5 <<- 1
            MGvar$scree.stress <<- as.vector(0)
            MGvar$scree.stress.T1 <<- as.vector(0)
            MGvar$scree.stress.T2 <<- as.vector(0)
            MGvar$scree.stress.T3 <<- as.vector(0)
            MGvar$scree.stress.T4 <<- as.vector(0)
            MGvar$scree.stress.T5 <<- as.vector(0)
            MGvar$screepoints.current <<- as.vector(0)
            MGvar$screepoints.current.T1 <<- as.vector(0)
            MGvar$screepoints.current.T2 <<- as.vector(0)
            MGvar$screepoints.current.T3 <<- as.vector(0)
            MGvar$screepoints.current.T4 <<- as.vector(0)
            MGvar$screepoints.current.T5 <<- as.vector(0)
            MGvar$screepoints.best <<- as.vector(0)
            MGvar$screepoints.best.T1 <<- as.vector(0)
            MGvar$screepoints.best.T2 <<- as.vector(0)
            MGvar$screepoints.best.T3 <<- as.vector(0)
            MGvar$screepoints.best.T4 <<- as.vector(0)
            MGvar$screepoints.best.T5 <<- as.vector(0)
            MGvar$distmatinternal <<- as.matrix(0)
            MGvar$activeplot.title <<- " "
            MGvar$activeplot.title.show <<- "yes"
            MGvar$activeplot.title.show.T1 <<- "yes"
            MGvar$activeplot.title.show.T2 <<- "yes"
            MGvar$activeplot.title.show.T3 <<- "yes"
            MGvar$activeplot.title.show.T4 <<- "yes"
            MGvar$activeplot.title.show.T5 <<- "yes"
            MGvar$activeplot.distmeas <<- "yes"
            MGvar$activeplot.distmeas.T1 <<- "yes"
            MGvar$activeplot.distmeas.T2 <<- "yes"
            MGvar$activeplot.distmeas.T3 <<- "yes"
            MGvar$activeplot.distmeas.T4 <<- "yes"
            MGvar$activeplot.distmeas.T5 <<- "yes"
            MGvar$activeplot.xlab <<- ""
            MGvar$activeplot.xlab.T1 <<- ""
            MGvar$activeplot.xlab.T2 <<- ""
            MGvar$activeplot.xlab.T3 <<- ""
            MGvar$activeplot.xlab.T4 <<- ""
            MGvar$activeplot.xlab.T5 <<- ""
            MGvar$activeplot.ylab <<- ""
            MGvar$activeplot.ylab.T1 <<- ""
            MGvar$activeplot.ylab.T2 <<- ""
            MGvar$activeplot.ylab.T3 <<- ""
            MGvar$activeplot.ylab.T4 <<- ""
            MGvar$activeplot.ylab.T5 <<- ""
            MGvar$activeplot.bg <<- "white"
            MGvar$activeplot.bg.temp <<- "white"
            MGvar$activeplot.bg.T1 <<- "white"
            MGvar$activeplot.bg.T2 <<- "white"
            MGvar$activeplot.bg.T3 <<- "white"
            MGvar$activeplot.bg.T4 <<- "white"
            MGvar$activeplot.bg.T5 <<- "white"
            MGvar$activeplot.cex <<- 0.6
            MGvar$activeplot.cex.T1 <<- 0.6
            MGvar$activeplot.cex.T2 <<- 0.6
            MGvar$activeplot.cex.T3 <<- 0.6
            MGvar$activeplot.cex.T4 <<- 0.6
            MGvar$activeplot.cex.T5 <<- 0.6
            MGvar$activeplot.labs <<- "yes"
            MGvar$activeplot.labs.T1 <<- "yes"
            MGvar$activeplot.labs.T2 <<- "yes"
            MGvar$activeplot.labs.T3 <<- "yes"
            MGvar$activeplot.labs.T4 <<- "yes"
            MGvar$activeplot.labs.T5 <<- "yes"
            MGvar$activeplot.showpoints <<- "no"
            MGvar$activeplot.showpoints.T1 <<- "no"
            MGvar$activeplot.showpoints.T2 <<- "no"
            MGvar$activeplot.showpoints.T3 <<- "no"
            MGvar$activeplot.showpoints.T4 <<- "no"
            MGvar$activeplot.showpoints.T5 <<- "no"
            MGvar$activeplot.pointcol <<- "black"
            MGvar$activeplot.pointcol.temp <<- "black"
            MGvar$activeplot.pointcol.T1 <<- "black"
            MGvar$activeplot.pointcol.T2 <<- "black"
            MGvar$activeplot.pointcol.T3 <<- "black"
            MGvar$activeplot.pointcol.T4 <<- "black"
            MGvar$activeplot.pointcol.T5 <<- "black"
            MGvar$activeplot.type <<- 1
            MGvar$activeplot.type.T1 <<- 1
            MGvar$activeplot.type.T2 <<- 1
            MGvar$activeplot.type.T3 <<- 1
            MGvar$activeplot.type.T4 <<- 1
            MGvar$activeplot.type.T5 <<- 1
            MGvar$activeplot.xaxt <<- "n"
            MGvar$activeplot.xaxt.T1 <<- "n"
            MGvar$activeplot.xaxt.T2 <<- "n"
            MGvar$activeplot.xaxt.T3 <<- "n"
            MGvar$activeplot.xaxt.T4 <<- "n"
            MGvar$activeplot.xaxt.T5 <<- "n"
            MGvar$activeplot.yaxt <<- "n"
            MGvar$activeplot.yaxt.T1 <<- "n"
            MGvar$activeplot.yaxt.T2 <<- "n"
            MGvar$activeplot.yaxt.T3 <<- "n"
            MGvar$activeplot.yaxt.T4 <<- "n"
            MGvar$activeplot.yaxt.T5 <<- "n"
            MGvar$activeplot.axescol <<- "black"
            MGvar$activeplot.axescol.temp <<- "black"
            MGvar$activeplot.axescol.T1 <<- "black"
            MGvar$activeplot.axescol.T2 <<- "black"
            MGvar$activeplot.axescol.T3 <<- "black"
            MGvar$activeplot.axescol.T4 <<- "black"
            MGvar$activeplot.axescol.T5 <<- "black"
            MGvar$zoomedplot.title.show <<- "yes"
            MGvar$zoomedplot.distmeas <<- "yes"
            MGvar$zoomedplot.xlab <<- ""
            MGvar$zoomedplot.ylab <<- ""
            MGvar$zoomedplot.bg <<- "white"
            MGvar$zoomedplot.bg.temp <<- "white"
            MGvar$zoomedplot.cex <<- 1
            MGvar$zoomedplot.labs <<- "yes"
            MGvar$zoomedplot.showpoints <<- "no"
            MGvar$zoomedplot.pointcol <<- "black"
            MGvar$zoomedplot.pointcol.temp <<- "black"
            MGvar$zoomedplot.type <<- 1
            MGvar$zoomedplot.yaxt <<- "n"
            MGvar$zoomedplot.yaxt <<- "n"
            MGvar$zoomedplot.axescol <<- "black"
            MGvar$zoomedplot.axescol.temp <<- "black"
            MGvar$zoomedplot.showzoom <<- "yes"
            MGvar$shepplot.title.show <<- "yes"
            MGvar$shepplot.labs.show <<- "yes"
            MGvar$shepplot.leg.show <<- "no"
            MGvar$shepplot.bg <<- "white"
            MGvar$shepplot.bg.temp <<- "white"
            MGvar$shepplot.showpoints <<- "yes"
            MGvar$shepplot.cex <<- 0.6
            MGvar$shepplot.type <<- 1
            MGvar$shepplot.pointcol.temp <<- "black"
            MGvar$shepplot.pointcol <<- "black"
            MGvar$shepplot.curve.show <<- "yes"
            MGvar$shepplot.curve.type <<- 1
            MGvar$shepplot.curvecol.temp <<- "red"
            MGvar$shepplot.curvecol <<- "red"
            MGvar$shepplot.Axes.xaxt <<- "s"
            MGvar$shepplot.Axes.yaxt <<- "s"
            MGvar$shepplot.Axescol.temp <<- "black"
            MGvar$shepplot.Axescol <<- "black"
            MGvar$screeplot.title.show <<- "yes"
            MGvar$screeplot.labs.show <<- "yes"
            MGvar$screeplot.leg.show <<- "yes"
            MGvar$screeplot.bg <<- "white"
            MGvar$screeplot.bg.temp <<- "white"
            MGvar$screeplot.points.show <<- "no"
            MGvar$screeplot.Cdim.show <<- "yes"
            MGvar$screeplot.Odim.show <<- "yes"
            MGvar$screeplot.Ccol <<- "red"
            MGvar$screeplot.Ccol.temp <<- "red"
            MGvar$screeplot.Ocol <<- "blue"
            MGvar$screeplot.Ocol.temp <<- "blue"
            MGvar$screeplot.curve.show <<- "yes"
            MGvar$screeplot.curve.type <<- 1
            MGvar$screeplot.curvecol <<- "black"
            MGvar$screeplot.curvecol.temp <<- "black"
            MGvar$screeplot.Cline.show <<- "yes"
            MGvar$screeplot.Oline.show <<- "yes"
            MGvar$screeplot.Axes.xaxt <<- "s"
            MGvar$screeplot.Axes.yaxt <<- "s"
            MGvar$screeplot.Axescol <<- "black"
            MGvar$screeplot.Axescol.temp <<- "black"
            if (popzoomswitch == "on") {
                tkdestroy(MGcomp$zoomplottt)
            }
            MGvar$datatitle <<- " "
            tclvalue(labelText) <<- paste("No Active Dataset")
            tclvalue(MGvar$datnam) <<- paste("Current name of active data is:               ", 
                MGvar$datatitle)
            MGvar$MDSStress <<- ""
            MGvar$MDSStress.T1 <<- "-"
            MGvar$MDSStress.T2 <<- "-"
            MGvar$MDSStress.T3 <<- "-"
            MGvar$MDSStress.T4 <<- "-"
            MGvar$MDSStress.T5 <<- "-"
            MGvar$dMeas <<- " "
            MGvar$dMeas.T1 <<- "-"
            MGvar$dMeas.T2 <<- "-"
            MGvar$dMeas.T3 <<- "-"
            MGvar$dMeas.T4 <<- "-"
            MGvar$dMeas.T5 <<- "-"
            MGvar$MDStype = ""
            MGvar$MDStype.T1 <<- "-"
            MGvar$MDStype.T2 <<- "-"
            MGvar$MDStype.T3 <<- "-"
            MGvar$MDStype.T4 <<- "-"
            MGvar$MDStype.T5 <<- "-"
            tableupdate()
            MGvar$Zoom.Main.val <<- tclVar("1")
            MGvar$Zoom.Dist.val <<- tclVar("1")
            MGvar$Zoom.Ylab.val <<- tclVar("0")
            MGvar$Zoom.Xlab.val <<- tclVar("0")
            MGvar$Zoom.Points.val <<- tclVar("0")
            MGvar$Zoom.Labels.val <<- tclVar("1")
            MGvar$Zoom.PT.var <<- tclVar("Empty Circles")
            MGvar$Zoom.AxesMeas.val <<- tclVar("0")
            MGvar$Zoom.ShowZoom.val <<- tclVar("1")
            MGvar$Scree.Main.val <<- tclVar("1")
            MGvar$Scree.Lab.val <<- tclVar("1")
            MGvar$Scree.Leg.val <<- tclVar("1")
            MGvar$Scree.Points.val <<- tclVar("0")
            MGvar$Scree.CDim.val <<- tclVar("1")
            MGvar$Scree.ODim.val <<- tclVar("1")
            MGvar$Scree.Curve.val <<- tclVar("1")
            MGvar$Scree.LineT.val <<- tclVar("Solid Line")
            MGvar$Scree.CLine.val <<- tclVar("1")
            MGvar$Scree.OLine.val <<- tclVar("1")
            MGvar$Scree.Axes.val <<- tclVar("1")
            MGvar$Shep.Main.val <<- tclVar("1")
            MGvar$Shep.Lab.val <<- tclVar("1")
            MGvar$Shep.Leg.val <<- tclVar("0")
            MGvar$Shep.Points.val <<- tclVar("1")
            MGvar$Shep.PS.var <<- tclVar(MGvar$shepplot.cex)
            MGvar$Shep.PT.var <<- tclVar("Empty Circles")
            MGvar$Shep.Curve.val <<- tclVar("1")
            MGvar$Shep.LineT.val <<- tclVar("Solid Line")
            MGvar$Shep.Axes.val <<- tclVar("1")
            MGvar$Conf.Main.val <<- tclVar("1")
            MGvar$Conf.Dist.val <<- tclVar("1")
            MGvar$Conf.Ylab.val <<- tclVar("0")
            MGvar$Conf.Xlab.val <<- tclVar("0")
            MGvar$Conf.Points.val <<- tclVar("0")
            MGvar$Conf.Labels.val <<- tclVar("1")
            MGvar$Conf.PT.var <<- tclVar("Empty Circles")
            MGvar$Conf.AxesMeas.val <<- tclVar("0")
            if (MGvar$EnShep.switch == "on") {
                tkdestroy(MGcomp$EShep)
            }
            MGvar$EnShep.switch <<- "off"
            if (EnScree.switch == "on") {
                tkdestroy(MGcomp$EScree)
            }
            EnScree.switch <<- "off"
            MGvar$seczoomswitch <<- "off"
            popzoomswitch <<- "off"
            MGvar$weights.vec <<- as.matrix(0)
            MGvar$sigma <<- as.matrix(0)
            MGvar$threeDconf <<- as.matrix(0)
            MGvar$temp.dims <<- 2
            MGvar$activeplot.title <<- ""
            MGvar$activeplot.title.T1 <<- ""
            MGvar$activeplot.title.T2 <<- ""
            MGvar$activeplot.title.T3 <<- ""
            MGvar$activeplot.title.T4 <<- ""
            MGvar$activeplot.title.T5 <<- ""
            MGvar$activeshepplot.title <<- " "
            MGvar$activescreeplot.title <<- " "
            MGvar$ActivePlottingTab <<- tclVar("Tab1")
            MGvar$rb.PT.Value <<- tclVar("PTab1")
            MGvar$MDS.dimensions <<- 2
            MGvar$Tab1.zoomedswitch <<- "off"
            MGvar$Tab2.zoomedswitch <<- "off"
            MGvar$Tab3.zoomedswitch <<- "off"
            MGvar$Tab4.zoomedswitch <<- "off"
            MGvar$Tab5.zoomedswitch <<- "off"
            MGvar$newCoords <<- as.matrix(0)
            MGvar$newCoords.T1 <<- as.matrix(0)
            MGvar$newCoords.T2 <<- as.matrix(0)
            MGvar$newCoords.T3 <<- as.matrix(0)
            MGvar$newCoords.T4 <<- as.matrix(0)
            MGvar$newCoords.T5 <<- as.matrix(0)
            MGvar$indexLabeled <<- c()
            MGvar$indexLabeled.T1 <<- c()
            MGvar$indexLabeled.T2 <<- c()
            MGvar$indexLabeled.T3 <<- c()
            MGvar$indexLabeled.T4 <<- c()
            MGvar$indexLabeled.T5 <<- c()
            MGvar$labeledPoints <<- list()
            if (MGvar$EnActivePlot.switch.T1 == "on") {
                tkdestroy(MGcomp$EActive)
            }
            if (MGvar$EnActivePlot.switch.T2 == "on") {
                tkdestroy(MGcomp$EActive)
            }
            if (MGvar$EnActivePlot.switch.T3 == "on") {
                tkdestroy(MGcomp$EActive)
            }
            if (MGvar$EnActivePlot.switch.T4 == "on") {
                tkdestroy(MGcomp$EActive)
            }
            if (MGvar$EnActivePlot.switch.T5 == "on") {
                tkdestroy(MGcomp$EActive)
            }
            MGvar$EnActivePlot.switch.T1 <<- "off"
            MGvar$EnActivePlot.switch.T2 <<- "off"
            MGvar$EnActivePlot.switch.T3 <<- "off"
            MGvar$EnActivePlot.switch.T4 <<- "off"
            MGvar$EnActivePlot.switch.T5 <<- "off"
            tkconfigure(mytt, cursor = "arrow")
            tkrreplot(img)
            tkrreplot(img2)
            tkrreplot(img3)
            tkrreplot(img4)
            tkrreplot(img5)
            tabplot()
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            tkrreplot(imgscree)
            tkrreplot(imgseczoom)
        }
        else {
        }
    }
    PerformProcrustes <- function() {
        PProc <- tktoplevel()
        tkwm.resizable(PProc, "0", "0")
        tkwm.deiconify(PProc)
        tkwm.title(PProc, "Procrustes Analysis")
        tkwm.geometry(PProc, "350x400")
        PProccanvas <- tkcanvas(PProc, width = 420, height = 450, 
            bg = col.sec)
        tkplace(PProccanvas, relx = 0, rely = 0, `in` = PProc)
        PProcNB <- tk2notebook(PProc, tabs = NULL)
        PProc.1 <- tk2frame(PProcNB)
        tkadd(PProcNB, PProc.1, text = "Procrustes")
        framePProc.1 <- tkwidget(PProc.1, "TitleFrame", text = "Perform Procrustes", 
            background = "white")
        tkplace(framePProc.1, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = PProc.1)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(framePProc.1, text = "Procrustes Analysis is a process by which one set of \ncoordinates is manipulated to match another as best as \npossible. This function will plot the manipulated set of\ncoordinates with the target coordinates to show their \nlevel of similiarity.", 
            font = fontsmall), relx = 0.08, rely = 0.1, `in` = framePProc.1)
        tkplace(tklabel(framePProc.1, text = "First Set of Coordinates", 
            background = "white"), relx = 0.1, rely = 0.43, `in` = framePProc.1)
        tkplace(tklabel(framePProc.1, text = "(Coordinates to be transformed)", 
            font = fontsmall), relx = 0.1, rely = 0.48, `in` = framePProc.1)
        PProc.C1.spin <- tk2spinbox(PProc, values = c("Plot5", 
            "Plot4", "Plot3", "Plot2", "Plot1"), width = 8)
        PProc.C1.val <- tclVar("Plot1")
        tkconfigure(PProc.C1.spin, textvariable = PProc.C1.val)
        tkplace(PProc.C1.spin, relx = 0.7, rely = 0.45, `in` = framePProc.1)
        tkplace(tklabel(framePProc.1, text = "Second Set of Coordinates", 
            background = "white"), relx = 0.1, rely = 0.63, `in` = framePProc.1)
        tkplace(tklabel(framePProc.1, text = "(Coordinates to be matched)", 
            font = fontsmall), relx = 0.1, rely = 0.68, `in` = framePProc.1)
        PProc.C2.spin <- tk2spinbox(PProc, values = c("Plot5", 
            "Plot4", "Plot3", "Plot2", "Plot1"), width = 8)
        PProc.C2.val <- tclVar("Plot2")
        tkconfigure(PProc.C2.spin, textvariable = PProc.C2.val)
        tkplace(PProc.C2.spin, relx = 0.7, rely = 0.65, `in` = framePProc.1)
        OnProcPlot <- function() {
            C1 <- tclvalue(PProc.C1.val)
            if (C1 == "Plot1") {
                ProcConf1 <<- MGvar$MDSmat.T1
                tclvalue(PProc.C1.val) <<- "Plot1"
            }
            if (C1 == "Plot2") {
                ProcConf1 <<- MGvar$MDSmat.T2
                tclvalue(PProc.C1.val) <<- "Plot2"
            }
            if (C1 == "Plot3") {
                ProcConf1 <<- MGvar$MDSmat.T3
                tclvalue(PProc.C1.val) <<- "Plot3"
            }
            if (C1 == "Plot4") {
                ProcConf1 <<- MGvar$MDSmat.T4
                tclvalue(PProc.C1.val) <<- "Plot4"
            }
            if (C1 == "Plot5") {
                ProcConf1 <<- MGvar$MDSmat.T5
                tclvalue(PProc.C1.val) <<- "Plot5"
            }
            C2 <- tclvalue(PProc.C2.val)
            if (C2 == "Plot1") {
                ProcConf2 <<- MGvar$MDSmat.T1
                tclvalue(PProc.C2.val) <<- "Plot1"
            }
            if (C2 == "Plot2") {
                ProcConf2 <<- MGvar$MDSmat.T2
                tclvalue(PProc.C2.val) <<- "Plot2"
            }
            if (C2 == "Plot3") {
                ProcConf2 <<- MGvar$MDSmat.T3
                tclvalue(PProc.C2.val) <<- "Plot3"
            }
            if (C2 == "Plot4") {
                ProcConf2 <<- MGvar$MDSmat.T4
                tclvalue(PProc.C2.val) <<- "Plot4"
            }
            if (C2 == "Plot5") {
                ProcConf2 <<- MGvar$MDSmat.T5
                tclvalue(PProc.C2.val) <<- "Plot5"
            }
            if (MGvar$GenSet.ClearIL == "yes") {
                MGvar$Proc.indexLabeled <<- c()
            }
            PerformP <- TRUE
            if (nrow(ProcConf1) != nrow(ProcConf2)) {
                PerformP <- FALSE
            }
            if (PerformP) {
                for (i in 1:nrow(ProcConf1)) {
                  if (rownames(ProcConf1)[i] != rownames(ProcConf2)[i]) {
                    PerformP <- FALSE
                  }
                }
            }
            if (!PerformP) {
                tkmessageBox(message = "Configuration Points to do not match! Unable to perform Procrustes Analyses.")
            }
            if (PerformP) {
                tkrreplot(procimg, function() ProcrustesMDSGUIplot(conf1 = ProcConf1, 
                  conf2 = ProcConf2, tab1 = C1, tab2 = C2))
                tk2notetab.select(myPlottingNB, "Procrustes")
            }
            tkdestroy(PProc)
        }
        tkplace(tkbutton(PProc, text = "Plot", width = 15, command = function() OnProcPlot()), 
            relx = 0.325, rely = 0.85, `in` = framePProc.1)
        OnOK.Proc <- function() {
            tkdestroy(PProc)
        }
        OnCancel.Proc <- function() {
            tkdestroy(PProc)
        }
        tkplace(PProcNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = PProc)
        tkplace(tkbutton(PProc, text = "OK", width = 15, command = function() OnOK.Proc()), 
            relx = 0.15, rely = 0.9, `in` = PProc)
        tkplace(tkbutton(PProc, text = "Cancel", width = 15, 
            command = function() OnCancel.Proc()), relx = 0.55, 
            rely = 0.9, `in` = PProc)
    }
    GeneralSettings <- function() {
        GenSet <- tktoplevel()
        tkwm.resizable(GenSet, "0", "0")
        tkwm.deiconify(GenSet)
        tkwm.geometry(GenSet, "420x450")
        tkwm.title(GenSet, "General Settings")
        GenSetcanvas <- tkcanvas(GenSet, width = "420", height = "450", 
            bg = col.sec)
        tkplace(GenSetcanvas, `in` = GenSet)
        GenSetNB <- tk2notebook(GenSet, tabs = NULL)
        GSGen <- tk2frame(GenSetNB)
        tkadd(GenSetNB, GSGen, text = "General")
        frameGSComp <- tkwidget(GSGen, "TitleFrame", text = "Computation Options", 
            background = "white")
        tkplace(frameGSComp, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.67, `in` = GSGen)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frameGSComp, text = "Calculations within MDSGUI can be very time consuming \nwhen datasets become large. This computation time can be\ndecreased if various plots are ommited from the output.\nSelect below which plots you wish to be included during\nthe MDS procedure.", 
            font = fontsmall), relx = 0.1, rely = 0.1, `in` = frameGSComp)
        tkplace(tklabel(frameGSComp, text = "MDS Point Configuration", 
            background = "white"), relx = 0.1, rely = 0.45, `in` = frameGSComp)
        GS.MDSConfig.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.MDSConfig.CB, variable = MGvar$GS.MDSConfig.val)
        tkplace(GS.MDSConfig.CB, relx = 0.75, rely = 0.45, `in` = frameGSComp)
        tkplace(tklabel(frameGSComp, text = "Shepard Plot", background = "white"), 
            relx = 0.1, rely = 0.58, `in` = frameGSComp)
        GS.Shep.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Shep.CB, variable = MGvar$GS.Shep.val)
        tkplace(GS.Shep.CB, relx = 0.75, rely = 0.58, `in` = frameGSComp)
        tkplace(tklabel(frameGSComp, text = "Stress Plot", background = "white"), 
            relx = 0.1, rely = 0.71, `in` = frameGSComp)
        GS.Stress.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Stress.CB, variable = MGvar$GS.Stress.val)
        tkplace(GS.Stress.CB, relx = 0.75, rely = 0.71, `in` = frameGSComp)
        tkplace(tklabel(frameGSComp, text = "Scree Plot", background = "white"), 
            relx = 0.1, rely = 0.84, `in` = frameGSComp)
        GS.Scree.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Scree.CB, variable = MGvar$GS.Scree.val)
        tkplace(GS.Scree.CB, relx = 0.75, rely = 0.84, `in` = frameGSComp)
        frameGSWin <- tkwidget(GSGen, "TitleFrame", text = "Window Options", 
            background = "white")
        tkplace(frameGSWin, relx = 0.02, relwidth = 0.96, rely = 0.71, 
            relheight = 0.27, `in` = GSGen)
        tkplace(tklabel(frameGSWin, text = "Multiple popped-out windows can be messy. Choose below\n whether or not to have multiple windows or use only one.", 
            font = fontsmall), relx = 0.08, rely = 0.2, `in` = frameGSWin)
        tkplace(tklabel(frameGSWin, text = "Re-Use Popped Out Window", 
            background = "white"), relx = 0.1, rely = 0.65, `in` = frameGSWin)
        GS.Win.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Win.CB, variable = MGvar$GS.Win.val)
        tkplace(GS.Win.CB, relx = 0.75, rely = 0.65, `in` = frameGSWin)
        GSCon <- tk2frame(GenSetNB)
        tkadd(GenSetNB, GSCon, text = "Convergence")
        frameGSConMDS <- tkwidget(GSCon, "TitleFrame", text = "MDS Convergence Settings", 
            background = "white")
        tkplace(frameGSConMDS, relx = 0.02, relwidth = 0.96, 
            rely = 0.02, relheight = 0.96, `in` = GSCon)
        tkplace(tklabel(frameGSConMDS, text = "MDS Relative Convergence", 
            background = "white"), relx = 0.1, rely = 0.15, `in` = frameGSConMDS)
        GS.MDStol.val <- tclVar(MGvar$MDS.tol)
        GS.MDStol.Input = tkentry(GenSet, width = 15, textvariable = GS.MDStol.val, 
            border = 2)
        tkplace(GS.MDStol.Input, relx = 0.65, rely = 0.15, `in` = frameGSConMDS)
        tkplace(tklabel(frameGSConMDS, text = "Maximum MDS Iterations", 
            background = "white"), relx = 0.1, rely = 0.3, `in` = frameGSConMDS)
        GS.MDSiter.val <- tclVar(MGvar$MDS.iter.max)
        GS.MDSiter.Input = tkentry(GenSet, width = 15, textvariable = GS.MDSiter.val, 
            border = 2)
        tkplace(GS.MDSiter.Input, relx = 0.65, rely = 0.3, `in` = frameGSConMDS)
        GSGraph <- tk2frame(GenSetNB)
        tkadd(GenSetNB, GSGraph, text = "Graphical")
        frameGSGraph <- tkwidget(GSGraph, "TitleFrame", text = "Graphing Settings", 
            background = "white")
        tkplace(frameGSGraph, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheight = 0.47, `in` = GSGraph)
        tkplace(tklabel(frameGSGraph, text = "MDS plot defaults with point labels", 
            background = "white"), relx = 0.1, rely = 0.1, `in` = frameGSGraph)
        GS.Lab.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Lab.CB, variable = MGvar$GS.Lab.val)
        tkplace(GS.Lab.CB, relx = 0.8, rely = 0.1, `in` = frameGSGraph)
        tkplace(tklabel(frameGSGraph, text = "MDS plot defaults with points", 
            background = "white"), relx = 0.1, rely = 0.28, `in` = frameGSGraph)
        GS.Pt.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Pt.CB, variable = MGvar$GS.Pt.val)
        tkplace(GS.Pt.CB, relx = 0.8, rely = 0.28, `in` = frameGSGraph)
        tkplace(tklabel(frameGSGraph, text = "MDS plot defaults with Variable Axes", 
            background = "white"), relx = 0.1, rely = 0.46, `in` = frameGSGraph)
        GS.Reg.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Reg.CB, variable = MGvar$GS.Reg.val)
        tkplace(GS.Reg.CB, relx = 0.8, rely = 0.46, `in` = frameGSGraph)
        tkplace(tklabel(frameGSGraph, text = "Added Shepard points clear with new plot", 
            background = "white"), relx = 0.1, rely = 0.64, `in` = frameGSGraph)
        GS.ShepP.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.ShepP.CB, variable = MGvar$GS.ShepP.val)
        tkplace(GS.ShepP.CB, relx = 0.8, rely = 0.64, `in` = frameGSGraph)
        tkplace(tklabel(frameGSGraph, text = "Individual point colours clear with new plot", 
            background = "white"), relx = 0.1, rely = 0.82, `in` = frameGSGraph)
        GS.Pcol.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.Pcol.CB, variable = MGvar$GS.Pcol.val)
        tkplace(GS.Pcol.CB, relx = 0.8, rely = 0.82, `in` = frameGSGraph)
        frameGSLabels <- tkwidget(GSGraph, "TitleFrame", text = "Point Label Settings", 
            background = "white")
        tkplace(frameGSLabels, relx = 0.02, rely = 0.51, relwidth = 0.96, 
            relheight = 0.47, `in` = GSGraph)
        tkplace(tklabel(frameGSLabels, text = "Added point labels refer to the labels manually added to the\nplot by the user. By Default these are added above the chosen\npoint and are cleared when a new plot is created.", 
            font = fontsmall), relx = 0.1, rely = 0.2, `in` = frameGSLabels)
        tkplace(tklabel(frameGSLabels, text = "Clear labels index for new plots", 
            background = "white"), relx = 0.1, rely = 0.55, `in` = frameGSLabels)
        GS.IL.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.IL.CB, variable = MGvar$GS.IL.val)
        tkplace(GS.IL.CB, relx = 0.8, rely = 0.55, `in` = frameGSLabels)
        tkplace(tklabel(frameGSLabels, text = "Added Point Label Postion", 
            background = "white"), relx = 0.1, rely = 0.75, `in` = frameGSLabels)
        GS.IP.ComboBox <- tkwidget(GenSet, "ComboBox", editable = FALSE, 
            values = c("Left", "Bottom", "Right", "Top"), width = 10, 
            textvariable = MGvar$GS.IP.val)
        tkplace(GS.IP.ComboBox, relx = 0.65, rely = 0.75, `in` = frameGSLabels)
        GSVis <- tk2frame(GenSetNB)
        tkadd(GenSetNB, GSVis, text = "Visualisation")
        frameGSVis <- tkwidget(GSVis, "TitleFrame", text = "Visualisation Settings", 
            background = "white")
        tkplace(frameGSVis, relx = 0.02, rely = 0.02, relwidth = 0.96, 
            relheight = 0.96, `in` = GSVis)
        tkplace(tklabel(frameGSVis, text = "Many of the MDS procedures in MDSGUI are iterative\nprocesses. Due to this, it is possible to view the changing\nresults throughout the iterations. This may, in some cases,\nbe time consuming so these options are adjustable.", 
            font = fontsmall), relx = 0.1, rely = 0.1, `in` = frameGSVis)
        tkplace(tklabel(frameGSVis, text = "Update MDS configuration", 
            background = "white"), relx = 0.1, rely = 0.3, `in` = frameGSVis)
        GS.UpConf.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.UpConf.CB, variable = MGvar$GS.UpConf.val)
        tkplace(GS.UpConf.CB, relx = 0.75, rely = 0.3, `in` = frameGSVis)
        tkplace(tklabel(frameGSVis, text = "Update Shepard Plot", 
            background = "white"), relx = 0.1, rely = 0.45, `in` = frameGSVis)
        GS.UpShep.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.UpShep.CB, variable = MGvar$GS.UpShep.val)
        tkplace(GS.UpShep.CB, relx = 0.75, rely = 0.45, `in` = frameGSVis)
        tkplace(tklabel(frameGSVis, text = "Update Stress Plot", 
            background = "white"), relx = 0.1, rely = 0.6, `in` = frameGSVis)
        GS.UpStress.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.UpStress.CB, variable = MGvar$GS.UpStress.val)
        tkplace(GS.UpStress.CB, relx = 0.75, rely = 0.6, `in` = frameGSVis)
        tkplace(tklabel(frameGSVis, text = "MDS Iterations Progress-Bar", 
            background = "white"), relx = 0.1, rely = 0.75, `in` = frameGSVis)
        GS.UpProg.CB <- tk2checkbutton(GenSet)
        tkconfigure(GS.UpProg.CB, variable = MGvar$GS.UpProg.val)
        tkplace(GS.UpProg.CB, relx = 0.75, rely = 0.75, `in` = frameGSVis)
        tkplace(GenSetNB, relx = 0.05, rely = 0.01, relwidth = 0.9, 
            relheight = 0.85, `in` = GenSet)
        OnOK.GenSet <- function() {
            surechange = tkmessageBox(message = "Are you sure you want to make these changes?", 
                type = "yesno", default = "no")
            if (as.character(surechange) == "yes") {
                MDSConf = as.character(tclvalue(MGvar$GS.MDSConfig.val))
                if (MDSConf == "1") {
                  MGvar$GenSet.CalcMDS <<- "yes"
                  tclvalue(MGvar$GS.MDSConfig.val) <<- 1
                }
                if (MDSConf == "0") {
                  MGvar$GenSet.CalcMDS <<- "no"
                  tclvalue(MGvar$GS.MDSConfig.val) <<- 0
                }
                MDSShep = as.character(tclvalue(MGvar$GS.Shep.val))
                if (MDSShep == "1") {
                  MGvar$GenSet.CalcShep <<- "yes"
                  tclvalue(MGvar$GS.Shep.val) <<- 1
                }
                if (MDSShep == "0") {
                  MGvar$GenSet.CalcShep <<- "no"
                  tclvalue(MGvar$GS.Shep.val) <<- 0
                }
                MGvar$MDSStress = as.character(tclvalue(MGvar$GS.Stress.val))
                if (MGvar$MDSStress == "1") {
                  MGvar$GenSet.CalcStress <<- "yes"
                  tclvalue(MGvar$GS.Stress.val) <<- 1
                }
                if (MGvar$MDSStress == "0") {
                  MGvar$GenSet.CalcStress <<- "no"
                  tclvalue(MGvar$GS.Stress.val) <<- 0
                }
                MDSScree = as.character(tclvalue(MGvar$GS.Scree.val))
                if (MDSScree == "1") {
                  MGvar$GenSet.CalcScree <<- "yes"
                  tclvalue(MGvar$GS.Scree.val) <<- 1
                }
                if (MDSScree == "0") {
                  MGvar$GenSet.CalcScree <<- "no"
                  tclvalue(MGvar$GS.Scree.val) <<- 0
                }
                RUWin <- as.character(tclvalue(MGvar$GS.Win.val))
                if (RUWin == "1") {
                  MGvar$GenSet.RUWin <<- "yes"
                  tclvalue(MGvar$GS.Win.val) <<- 1
                }
                if (RUWin == "0") {
                  MGvar$GenSet.RUWin <<- "no"
                  tclvalue(MGvar$GS.Win.val) <<- 0
                }
                mdstol <- as.numeric(tclvalue(GS.MDStol.val))
                MGvar$MDS.tol <<- mdstol
                mdsiter <- as.numeric(tclvalue(GS.MDSiter.val))
                MGvar$MDS.iter.max <<- mdsiter
                tclvalue(iterrat) <<- paste(MGvar$MDS.iter.max, 
                  "iters")
                DLab = as.character(tclvalue(MGvar$GS.Lab.val))
                if (DLab == "1") {
                  MGvar$activeplot.labs <<- "yes"
                  MGvar$activeplot.labs.T1 <<- "yes"
                  MGvar$activeplot.labs.T2 <<- "yes"
                  MGvar$activeplot.labs.T3 <<- "yes"
                  MGvar$activeplot.labs.T4 <<- "yes"
                  MGvar$activeplot.labs.T5 <<- "yes"
                  MGvar$procplot.labs1 <<- "yes"
                  MGvar$procplot.labs2 <<- "yes"
                  MGvar$zoomedplot.labs <<- "yes"
                  MGvar$GenSet.DefLabs <<- "yes"
                  tclvalue(MGvar$GS.Lab.val) <<- 1
                  tclvalue(MGvar$Conf.Labels.val) <<- 1
                  tclvalue(MGvar$Proc.Labels1.val) <<- 1
                  tclvalue(MGvar$Proc.Labels2.val) <<- 1
                  tclvalue(MGvar$Zoom.Labels.val) <<- 1
                }
                if (DLab == "0") {
                  MGvar$activeplot.labs <<- "no"
                  MGvar$activeplot.labs.T1 <<- "no"
                  MGvar$activeplot.labs.T2 <<- "no"
                  MGvar$activeplot.labs.T3 <<- "no"
                  MGvar$activeplot.labs.T4 <<- "no"
                  MGvar$activeplot.labs.T5 <<- "no"
                  MGvar$procplot.labs1 <<- "no"
                  MGvar$procplot.labs2 <<- "no"
                  MGvar$zoomedplot.labs <<- "no"
                  MGvar$GenSet.DefLabs <<- "no"
                  tclvalue(MGvar$GS.Lab.val) <<- 0
                  tclvalue(MGvar$Conf.Labels.val) <<- 0
                  tclvalue(MGvar$Proc.Labels1.val) <<- 0
                  tclvalue(MGvar$Proc.Labels2.val) <<- 0
                  tclvalue(MGvar$Zoom.Labels.val) <<- 0
                }
                DPt = as.character(tclvalue(MGvar$GS.Pt.val))
                if (DPt == "1") {
                  MGvar$activeplot.showpoints <<- "yes"
                  MGvar$activeplot.showpoints.T1 <<- "yes"
                  MGvar$activeplot.showpoints.T2 <<- "yes"
                  MGvar$activeplot.showpoints.T3 <<- "yes"
                  MGvar$activeplot.showpoints.T4 <<- "yes"
                  MGvar$activeplot.showpoints.T5 <<- "yes"
                  MGvar$procplot.showpoints1 <<- "yes"
                  MGvar$procplot.showpoints2 <<- "yes"
                  MGvar$zoomedplot.showpoints <<- "yes"
                  MGvar$GenSet.DefPts <<- "yes"
                  tclvalue(MGvar$GS.Pt.val) <<- 1
                  tclvalue(MGvar$Conf.Points.val) <<- 1
                  tclvalue(MGvar$Proc.Points1.val) <<- 1
                  tclvalue(MGvar$Proc.Points2.val) <<- 1
                  tclvalue(MGvar$Zoom.Points.val) <<- 1
                }
                if (DPt == "0") {
                  MGvar$activeplot.showpoints <<- "no"
                  MGvar$activeplot.showpoints.T1 <<- "no"
                  MGvar$activeplot.showpoints.T2 <<- "no"
                  MGvar$activeplot.showpoints.T3 <<- "no"
                  MGvar$activeplot.showpoints.T4 <<- "no"
                  MGvar$activeplot.showpoints.T5 <<- "no"
                  MGvar$procplot.showpoints1 <<- "no"
                  MGvar$procplot.showpoints2 <<- "no"
                  MGvar$zoomedplot.showpoints <<- "no"
                  MGvar$GenSet.DefPts <<- "no"
                  tclvalue(MGvar$GS.Pt.val) <<- 0
                  tclvalue(MGvar$Conf.Points.val) <<- 0
                  tclvalue(MGvar$Proc.Points1.val) <<- 0
                  tclvalue(MGvar$Proc.Points2.val) <<- 0
                  tclvalue(MGvar$Zoom.Points.val) <<- 0
                }
                DReg = as.character(tclvalue(MGvar$GS.Reg.val))
                if (DReg == "1") {
                  MGvar$activeplot.showreg <<- "yes"
                  MGvar$activeplot.showreg.T1 <<- "yes"
                  MGvar$activeplot.showreg.T2 <<- "yes"
                  MGvar$activeplot.showreg.T3 <<- "yes"
                  MGvar$activeplot.showreg.T4 <<- "yes"
                  MGvar$activeplot.showreg.T5 <<- "yes"
                  MGvar$procplot.showreg1 <<- "yes"
                  MGvar$procplot.showreg2 <<- "yes"
                  MGvar$zoomedplot.showreg <<- "yes"
                  MGvar$GenSet.DefReg <<- "yes"
                  tclvalue(MGvar$GS.Reg.val) <<- 1
                  tclvalue(MGvar$Conf.RegLine.val) <<- 1
                  tclvalue(MGvar$Proc.RegLine1.val) <<- 1
                  tclvalue(MGvar$Proc.RegLine2.val) <<- 1
                  tclvalue(MGvar$Zoom.RegLine.val) <<- 1
                }
                if (DReg == "0") {
                  MGvar$activeplot.showreg <<- "no"
                  MGvar$activeplot.showreg.T1 <<- "no"
                  MGvar$activeplot.showreg.T2 <<- "no"
                  MGvar$activeplot.showreg.T3 <<- "no"
                  MGvar$activeplot.showreg.T4 <<- "no"
                  MGvar$activeplot.showreg.T5 <<- "no"
                  MGvar$procplot.showreg1 <<- "no"
                  MGvar$procplot.showreg2 <<- "no"
                  MGvar$zoomedplot.showreg <<- "no"
                  MGvar$GenSet.DefReg <<- "no"
                  tclvalue(MGvar$GS.Reg.val) <<- 0
                  tclvalue(MGvar$Conf.RegLine.val) <<- 0
                  tclvalue(MGvar$Proc.RegLine1.val) <<- 0
                  tclvalue(MGvar$Proc.RegLine2.val) <<- 0
                  tclvalue(MGvar$Zoom.RegLine.val) <<- 0
                }
                ShepL = as.character(tclvalue(MGvar$GS.ShepP.val))
                if (ShepL == "0") {
                  MGvar$GenSet.KeepShep <<- "yes"
                  tclvalue(MGvar$GS.ShepP.val) <<- 0
                }
                if (ShepL == "1") {
                  MGvar$GenSet.KeepShep <<- "no"
                  tclvalue(MGvar$GS.ShepP.val) <<- 1
                }
                ClearC = as.character(tclvalue(MGvar$GS.Pcol.val))
                if (ClearC == "1") {
                  MGvar$GenSet.ClearCols <<- "yes"
                  tclvalue(MGvar$GS.Pcol.val) <<- 1
                }
                if (ClearC == "0") {
                  MGvar$GenSet.ClearCols <<- "no"
                  tclvalue(MGvar$GS.Pcol.val) <<- 0
                }
                CIL = as.character(tclvalue(MGvar$GS.IL.val))
                if (CIL == "1") {
                  MGvar$GenSet.ClearIL <<- "yes"
                  tclvalue(MGvar$GS.IL.val) <<- 1
                }
                if (CIL == "0") {
                  MGvar$GenSet.ClearIL <<- "no"
                  tclvalue(MGvar$GS.IL.val) <<- 0
                }
                ILPos = as.character(tclvalue(MGvar$GS.IP.val))
                if (ILPos == "Bottom") {
                  MGvar$GenSet.ILPos <<- 1
                  tclvalue(MGvar$GS.IP.val) <<- "Bottom"
                }
                if (ILPos == "Left") {
                  MGvar$GenSet.ILPos <<- 2
                  tclvalue(MGvar$GS.IP.val) <<- "Left"
                }
                if (ILPos == "Top") {
                  MGvar$GenSet.ILPos <<- 3
                  tclvalue(MGvar$GS.IP.val) <<- "Top"
                }
                if (ILPos == "Right") {
                  MGvar$GenSet.ILPos <<- 4
                  tclvalue(MGvar$GS.IP.val) <<- "Right"
                }
                UC <- as.character(tclvalue(MGvar$GS.UpConf.val))
                if (UC == "1") {
                  MGvar$GenSet.UpConf <<- "yes"
                  tclvalue(MGvar$GS.UpConf.val) <<- 1
                }
                if (UC == "0") {
                  MGvar$GenSet.UpConf <<- "no"
                  tclvalue(MGvar$GS.UpConf.val) <<- 0
                }
                USh <- as.character(tclvalue(MGvar$GS.UpShep.val))
                if (USh == "1") {
                  MGvar$GenSet.UpShep <<- "yes"
                  tclvalue(MGvar$GS.UpShep.val) <<- 1
                }
                if (USh == "0") {
                  MGvar$GenSet.UpShep <<- "no"
                  tclvalue(MGvar$GS.UpShep.val) <<- 0
                }
                USt <- as.character(tclvalue(MGvar$GS.UpStress.val))
                if (USt == "1") {
                  MGvar$GenSet.UpStress <<- "yes"
                  tclvalue(MGvar$GS.UpStress.val) <<- 1
                }
                if (USt == "0") {
                  MGvar$GenSet.UpStress <<- "no"
                  tclvalue(MGvar$GS.UpStress.val) <<- 0
                }
                UPb <- as.character(tclvalue(MGvar$GS.UpProg.val))
                if (UPb == "1") {
                  MGvar$GenSet.UpProg <<- "yes"
                  tclvalue(MGvar$GS.UpProg.val) <<- 1
                }
                if (UPb == "0") {
                  MGvar$GenSet.UpProg <<- "no"
                  tclvalue(MGvar$GS.UpProg.val) <<- 0
                }
                tkdestroy(GenSet)
            }
        }
        OnDefault.GenSet <- function() {
            MGvar$GS.MDSConfig.val <<- tclVar("1")
            MGvar$GS.Shep.val <<- tclVar("1")
            MGvar$GS.Scree.val <<- tclVar("1")
            MGvar$GS.Win.val <<- tclVar("0")
            MGvar$GS.Lab.val <<- tclVar("1")
            MGvar$GS.Pt.val <<- tclVar("0")
            MGvar$GS.Reg.val <<- tclVar("0")
            MGvar$GS.IL.val <<- tclVar("1")
            MGvar$GS.IP.val <<- tclVar("Top")
            MGvar$GS.UpConf.val <<- tclVar("1")
            MGvar$GS.UpStress.val <<- tclVar("1")
            MGvar$GS.UpShep.val <<- tclVar("1")
            MGvar$GS.UpProg.val <<- tclVar("1")
            MGvar$GenSet.CalcMDS <<- "yes"
            MGvar$GenSet.CalcShep <<- "yes"
            MGvar$GenSet.CalcScree <<- "yes"
            MGvar$GenSet.RUWin <<- "no"
            MGvar$GenSet.DefLabs <<- "yes"
            MGvar$GenSet.DefPts <<- "no"
            MGvar$GenSet.DefReg <<- "no"
            MGvar$GenSet.ClearIL <<- "yes"
            MGvar$GenSet.ILPos <<- 3
            MGvar$GenSet.UpConf <<- "yes"
            MGvar$GenSet.UpShep <<- "yes"
            MGvar$GenSet.UpStress <<- "yes"
            MGvar$GenSet.UpProg <<- "yes"
            tkdestroy(GenSet)
        }
        tkplace(tkbutton(GenSet, text = "OK", width = 15, command = function() OnOK.GenSet()), 
            relx = 0.15, rely = 0.9, `in` = GenSet)
        tkplace(tkbutton(GenSet, text = "Default", width = 15, 
            command = function() OnDefault.GenSet()), relx = 0.55, 
            rely = 0.9, `in` = GenSet)
        tkfocus(GenSet)
        tkwait.window(GenSet)
    }
    my.tk2edit = function(x, title = "Matrix Editor", header = NULL, 
        maxHeight = 600, maxWidth = 800, fontsize = 9, ...) {
        if (!is.tk()) 
            stop("Package Tk is required but not loaded")
        if (!inherits(tclRequire("Tktable", warn = FALSE), "tclObj")) 
            stop("Tcl package 'Tktable' must be installed first")
        .Tcl(paste("option add *Table.font {courier", fontsize, 
            "bold}"))
        old <- options(scipen = 7)
        on.exit(options(old))
        numrows = nrow(x)
        numcols = ncol(x)
        makeCharMat <- function(x) {
            mat <- matrix(unlist(x), nrow = nrow(as.matrix(x)))
            dm <- dim(mat)
            hasRownames <- length(rn <- rownames(x)) > 0
            hasColnames <- length(cn <- colnames(x)) > 0
            if (!hasRownames) 
                rn <- paste("[", 1:nrow(x), ",]", sep = "")
            if (!hasColnames) 
                cn <- paste("[,", 1:ncol(x), "]", sep = "")
            mat[] <- apply(unclass(mat), 2, format, justify = "right")
            mat <- rbind(cn, mat)
            mat <- cbind(c("", rn), mat)
            mat
        }
        fillTclArrayFromCharMat <- function(ta, cm) {
            for (j in 2:ncol(cm)) ta[[0, j - 1]] <- as.tclObj(cm[1, 
                j], drop = TRUE)
            for (i in 2:nrow(cm)) for (j in 1:ncol(cm)) ta[[i - 
                1, j - 1]] <- as.tclObj(cm[i, j], drop = TRUE)
        }
        tA <- tclArray()
        cmat <- makeCharMat(x)
        fillTclArrayFromCharMat(tA, cmat)
        tt <- tktoplevel()
        tkwm.title(tt, title)
        colwidths <- apply(cmat, 2, function(x) max(nchar(x)) + 
            1)
        nTableCols <- ncol(cmat)
        if ((moreWidth <- 60 - sum(colwidths)) > 0) {
            addEach <- moreWidth%/%length(colwidths)
            if (addEach < 5) 
                colwidths <- colwidths + addEach + 1
            else nTableCols <- nTableCols + ceiling(moreWidth/10)
        }
        tktable <- tkwidget(tt, "table", variable = tA, rows = nrow(cmat), 
            cols = nTableCols, titlerows = 1, titlecols = 1, 
            selecttitle = 1, anchor = "e", multiline = 0, selectmode = "extended", 
            rowseparator = dQuote("\n"), colseparator = dQuote("\t"), 
            background = "white", maxheight = maxHeight, maxwidth = maxWidth, 
            xscrollcommand = function(...) tkset(xscr, ...), 
            yscrollcommand = function(...) tkset(yscr, ...))
        xscr <- tkscrollbar(tt, orient = "horizontal", command = function(...) tkxview(tktable, 
            ...))
        yscr <- tkscrollbar(tt, command = function(...) tkyview(tktable, 
            ...))
        for (i in 1:ncol(cmat)) tcl(tktable, "width", i - 1, 
            colwidths[i])
        string <- "bind Table <BackSpace> {\n\t\tset ::tk::table::Priv(junk) [%W icursor]\n    \tif {[string compare {} $::tk::table::Priv(junk)] && $::tk::table::Priv(junk)} {\n\t\t%W delete active [expr {$::tk::table::Priv(junk)-1}]\n    \t}}"
        .Tcl(string)
        activeRow <- function() as.numeric(tkindex(tktable, "active", 
            "row"))
        activeCol <- function() as.numeric(tkindex(tktable, "active", 
            "col"))
        undoEdits <- function() {
            ta <- tclArray()
            fillTclArrayFromCharMat(ta, cmat)
            assign("tA", ta, inherits = TRUE)
            tkconfigure(tktable, variable = tA)
        }
        finish <- function() tkdestroy(tt)
        cancel <- function() {
            undoEdits()
            tkdestroy(tt)
        }
        insertRow <- function() {
            row <- activeRow()
            col <- activeCol()
            tkinsert(tktable, "rows", row, 1)
            newCell <- paste(row + 1, col, sep = ",")
            tkactivate(tktable, newCell)
            tksee(tktable, newCell)
            numrows <<- numrows + 1
        }
        insertCol <- function() {
            row <- activeRow()
            col <- activeCol()
            tkinsert(tktable, "cols", col, 1)
            newCell <- paste(row, col + 1, sep = ",")
            tkactivate(tktable, newCell)
            tksee(tktable, newCell)
            numcols <<- numcols + 1
        }
        deleteRow <- function() {
            if ((row <- activeRow()) != 0) {
                tkdelete(tktable, "rows", row, 1)
                numrows <<- numrows - 1
            }
        }
        deleteCol <- function() {
            if ((col <- activeCol()) != 0) {
                tkdelete(tktable, "cols", col, 1)
                numcols <<- numcols - 1
            }
        }
        copyRow <- function() {
            src <- activeRow()
            if (src != 0) {
                insertRow()
                dst <- activeRow()
                for (j in 0:(ncol(tA) - 1)) tA[[dst, j]] <- tA[[src, 
                  j]]
            }
        }
        copyCol <- function() {
            src <- activeCol()
            if (src != 0) {
                insertCol()
                dst <- activeCol()
                for (i in 0:(nrow(tA) - 1)) tA[[i, dst]] <- tA[[i, 
                  src]]
            }
        }
        finishButton <- tkbutton(tt, text = "Finish", command = finish)
        cancelButton <- tkbutton(tt, text = "Cancel", command = cancel)
        undoEditsButton <- tkbutton(tt, text = "Undo Edits", 
            command = undoEdits)
        insertRowButton <- tkbutton(tt, text = "Insert Row", 
            command = insertRow)
        copyRowButton <- tkbutton(tt, text = "Copy Row", command = copyRow)
        deleteRowButton <- tkbutton(tt, text = "Delete Row", 
            command = deleteRow)
        insertColButton <- tkbutton(tt, text = "Insert Col", 
            command = insertCol)
        copyColButton <- tkbutton(tt, text = "Copy Col", command = copyCol)
        deleteColButton <- tkbutton(tt, text = "Delete Col", 
            command = deleteCol)
        if (length(header) > 0) {
            for (label in header) tkgrid(tklabel(tt, text = label), 
                columnspan = 7, sticky = "nw")
        }
        tkgrid(tktable, yscr, columnspan = 8)
        tkgrid.configure(tktable, sticky = "news")
        tkgrid.configure(yscr, sticky = "nsw")
        tkgrid(xscr, sticky = "new", columnspan = 8)
        tkgrid(insertRowButton, copyRowButton, deleteRowButton, 
            sticky = "news")
        tkgrid(insertColButton, copyColButton, deleteColButton, 
            "x", cancelButton, undoEditsButton, finishButton, 
            sticky = "news")
        tkgrid.columnconfigure(tt, 3, weight = 1)
        tkgrid.rowconfigure(tt, length(header), weight = 1)
        tkactivate(tktable, "0,0")
        tktag.configure(tktable, "active", background = "lightyellow2")
        tktag.configure(tktable, "title", state = "normal")
        tkfocus(tt)
        tkwait.window(tt)
        outMat <- matrix("", nrow = numrows + 1, ncol = numcols + 
            1)
        for (i in 1:nrow(outMat)) {
            for (j in 1:ncol(outMat)) {
                val <- tA[[i - 1, j - 1]]
                if (is.null(val)) 
                  val <- ""
                else val <- tclvalue(val)
                outMat[i, j] <- val
            }
        }
        rn <- outMat[, 1][-1]
        cn <- outMat[1, ][-1]
        outMat <- outMat[-1, -1, drop = FALSE]
        badRownames <- c(grep("\\[.*\\]", rn), (1:length(rn))[is.na(rn)])
        if (length(badRownames) != length(rn)) {
            rn[badRownames] <- ""
            rownames(outMat) <- rn
        }
        badColnames <- c(grep("\\[.*\\]", cn), (1:length(cn))[is.na(cn)])
        if (length(badColnames) != length(cn)) {
            cn[badColnames] <- ""
            colnames(outMat) <- cn
        }
        mode(outMat) <- mode(x)
        Sys.sleep(0.1)
        return(outMat)
    }
    Destroytoplevel <- function() {
        QuitOption = tkmessageBox(message = "Are you sure you want to quit?", 
            type = "yesno", default = "yes")
        if (as.character(QuitOption) == "yes") {
            tkdestroy(myPlottingNB)
            tkdestroy(mySecondaryNB)
            tkdestroy(mytt)
        }
        else {
        }
    }
    simplezoomout <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$zoominrat.T1 <<- MGvar$zoominrat.T1 * 1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$zoominrat.T2 <<- MGvar$zoominrat.T2 * 1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$zoominrat.T3 <<- MGvar$zoominrat.T3 * 1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$zoominrat.T4 <<- MGvar$zoominrat.T4 * 1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$zoominrat.T5 <<- MGvar$zoominrat.T5 * 1.1
        }
        tabplot()
    }
    simplezoomin <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$zoominrat.T1 <<- MGvar$zoominrat.T1/1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$zoominrat.T2 <<- MGvar$zoominrat.T2/1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$zoominrat.T3 <<- MGvar$zoominrat.T3/1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$zoominrat.T4 <<- MGvar$zoominrat.T4/1.1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$zoominrat.T5 <<- MGvar$zoominrat.T5/1.1
        }
        tabplot()
    }
    moveupfunc <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$moveup.T1 <<- MGvar$moveup.T1 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$moveup.T2 <<- MGvar$moveup.T2 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$moveup.T3 <<- MGvar$moveup.T3 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$moveup.T4 <<- MGvar$moveup.T4 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$moveup.T5 <<- MGvar$moveup.T5 * 1.05
        }
        tabplot()
    }
    movedownfunc <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$movedown.T1 <<- MGvar$movedown.T1 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$movedown.T2 <<- MGvar$movedown.T2 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$movedown.T3 <<- MGvar$movedown.T3 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$movedown.T4 <<- MGvar$movedown.T4 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$movedown.T5 <<- MGvar$movedown.T5 * 1.05
        }
        tabplot()
    }
    moveleftfunc <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$moveleft.T1 <<- MGvar$moveleft.T1 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$moveleft.T2 <<- MGvar$moveleft.T2 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$moveleft.T3 <<- MGvar$moveleft.T3 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$moveleft.T4 <<- MGvar$moveleft.T4 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$moveleft.T5 <<- MGvar$moveleft.T5 * 1.05
        }
        tabplot()
    }
    moverightfunc <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$moveright.T1 <<- MGvar$moveright.T1 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$moveright.T2 <<- MGvar$moveright.T2 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$moveright.T3 <<- MGvar$moveright.T3 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$moveright.T4 <<- MGvar$moveright.T4 * 1.05
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$moveright.T5 <<- MGvar$moveright.T5 * 1.05
        }
        tabplot()
    }
    origpos <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$zoominrat.T1 <<- 1
            MGvar$moveup.T1 <<- 1
            MGvar$movedown.T1 <<- 1
            MGvar$moveleft.T1 <<- 1
            MGvar$moveright.T1 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$zoominrat.T2 <<- 1
            MGvar$moveup.T2 <<- 1
            MGvar$movedown.T2 <<- 1
            MGvar$moveleft.T2 <<- 1
            MGvar$moveright.T2 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$zoominrat.T3 <<- 1
            MGvar$moveup.T3 <<- 1
            MGvar$movedown.T3 <<- 1
            MGvar$moveleft.T3 <<- 1
            MGvar$moveright.T3 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$zoominrat.T4 <<- 1
            MGvar$moveup.T4 <<- 1
            MGvar$movedown.T4 <<- 1
            MGvar$moveleft.T4 <<- 1
            MGvar$moveright.T4 <<- 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$zoominrat.T5 <<- 1
            MGvar$moveup.T5 <<- 1
            MGvar$movedown.T5 <<- 1
            MGvar$moveleft.T5 <<- 1
            MGvar$moveright.T5 <<- 1
        }
        tabplot()
    }
    proc.origpos <- function() {
        MGvar$proc.zoominrat <<- 1
        MGvar$proc.moveup <<- 1
        MGvar$proc.movedown <<- 1
        MGvar$proc.moveleft <<- 1
        MGvar$proc.moveright <<- 1
        tkrreplot(procimg)
    }
    RegLines <- function() {
        Regtt <- tktoplevel()
        tkwm.resizable(Regtt, "0", "0")
        tkwm.deiconify(Regtt)
        tkwm.title(Regtt, "Regression Axes")
        tkwm.geometry(Regtt, "420x250")
        Regcanvas = tkcanvas(Regtt, width = 420, height = 250, 
            bg = col.sec)
        tkplace(Regcanvas, `in` = Regtt)
        Regframe <- tkwidget(Regtt, "TitleFrame", text = "Add Regression Axes", 
            background = "white")
        tkplace(Regframe, relx = 0.05, relwidth = 0.9, rely = 0.01, 
            relheight = 0.82, `in` = Regtt)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(Regframe, text = "Add a Regression Line to your plot relating to one of the \nvariables from your active dataset.", 
            font = fontsmall), relx = 0.1, rely = 0.1, `in` = Regframe)
        tkplace(tklabel(Regframe, text = "Dependent Variable", 
            background = "white"), relx = 0.1, rely = 0.4, `in` = Regframe)
        regvalues = colnames(MGvar$activedata)
        Reg.var.val <- tclVar(regvalues[1])
        Reg.var.ComboBox <- tkwidget(Regtt, "ComboBox", editable = FALSE, 
            values = regvalues, width = 15, textvariable = Reg.var.val)
        tkplace(Reg.var.ComboBox, relx = 0.58, rely = 0.4, `in` = Regframe)
        tkplace(tklabel(Regframe, text = "Replace Previous Regression Line?", 
            background = "white"), relx = 0.1, rely = 0.65, `in` = Regframe)
        Reg.RegRep.val <- tclVar("1")
        Reg.RegRep.CB <- tk2checkbutton(Regtt)
        tkconfigure(Reg.RegRep.CB, variable = Reg.RegRep.val)
        tkplace(Reg.RegRep.CB, relx = 0.8, rely = 0.65, `in` = Regframe)
        tkplace(tklabel(Regframe, text = "Label Regression Line on Plot?", 
            background = "white"), relx = 0.1, rely = 0.8, `in` = Regframe)
        Reg.RegLab.val <- tclVar("1")
        Reg.RegLab.CB <- tk2checkbutton(Regtt)
        tkconfigure(Reg.RegLab.CB, variable = Reg.RegLab.val)
        tkplace(Reg.RegLab.CB, relx = 0.8, rely = 0.8, `in` = Regframe)
        OnOK.Reg <- function() {
            for (i in 1:length(regvalues)) {
                if (as.character(tclvalue(Reg.var.val)) == regvalues[i]) {
                  RegLineCalculation(i)
                }
            }
            tkdestroy(Regtt)
        }
        OnCancel.Reg <- function() {
            tkdestroy(Regtt)
        }
        tkplace(tkbutton(Regtt, width = 15, text = "Add Reg Line", 
            command = function() OnOK.Reg()), relx = 0.15, rely = 0.86, 
            `in` = Regtt)
        tkplace(tkbutton(Regtt, width = 15, text = "Cancel", 
            command = function() OnCancel.Reg()), relx = 0.58, 
            rely = 0.86, `in` = Regtt)
    }
    RegLineCalculation <- function() {
        plot(MGvar$MDSmat)
        X = MGvar$MDSmat
        for (i in 1:ncol(MGvar$activedata)) {
            y = MGvar$activedata[, i]
            B = solve(t(X) %*% X) %*% t(X) %*% y
            abline(0, (B[2]/B[1]))
        }
    }
    PointtoMove <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$MDSmat[, 1] - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$MDSmat[, 2] - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        GetClosestPointMovePoint(xClick, yClick, imgXcoords, 
            imgYcoords)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", OnPlotLeftClick)
        }
        tkdestroy(MGcomp$MPoint)
    }
    TargetLocation <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$MDSmat[, 1] - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$MDSmat[, 2] - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        CoordsChange(as.numeric(xPlotCoord), as.numeric(yPlotCoord))
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", OnPlotLeftClick)
        }
        tkdestroy(MGcomp$LPoint)
        tkdestroy(MGcomp$MPoint)
    }
    CoordsChange <- function(newx, newy) {
        MGvar$prevMat <<- MGvar$MDSmat
        MGvar$MDSmat[MGvar$Movingindex.Closest, 1] <<- newx
        MGvar$MDSmat[MGvar$Movingindex.Closest, 2] <<- newy
        ActiveCoordMat()
        tabplot()
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
    }
    MovePoint <- function() {
        MGcomp$MPoint <<- tktoplevel()
        tkwm.title(MGcomp$MPoint, "Choose Point")
        tkwm.geometry(MGcomp$MPoint, "280x90")
        MPointcanvas = tkcanvas(MGcomp$MPoint, width = 300, height = 120, 
            bg = col.sec)
        tkplace(MPointcanvas, relx = 0, rely = 0, `in` = MGcomp$MPoint)
        frameMP <- tkwidget(MGcomp$MPoint, "TitleFrame", text = "Active Cursor", 
            background = "white")
        tkplace(frameMP, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = MGcomp$MPoint)
        tkplace(tklabel(frameMP, text = "The cursor is now active. Please\nuse the cursor to select the point\nwhich you would like to relocate.", 
            background = "white"), relx = 0.1, rely = 0.2, `in` = frameMP)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", PointtoMove)
            tkconfigure(img, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", PointtoMove)
            tkconfigure(img2, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", PointtoMove)
            tkconfigure(img3, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", PointtoMove)
            tkconfigure(img4, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", PointtoMove)
            tkconfigure(img5, cursor = "crosshair")
        }
        tkwait.window(MGcomp$MPoint)
        ActiveArrowCursor()
    }
    GetClosestPointMovePoint <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        MGvar$Movingindex.Closest <<- which.min(squared.Distance)
        MGcomp$LPoint <<- tktoplevel()
        tkwm.title(MGcomp$LPoint, "Choose Point")
        tkwm.geometry(MGcomp$LPoint, "280x90")
        LPointcanvas = tkcanvas(MGcomp$LPoint, width = 300, height = 120, 
            bg = col.sec)
        tkplace(LPointcanvas, relx = 0, rely = 0, `in` = MGcomp$LPoint)
        frameLP <- tkwidget(MGcomp$LPoint, "TitleFrame", text = "Active Cursor", 
            background = "white")
        tkplace(frameLP, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = MGcomp$LPoint)
        tkplace(tklabel(frameLP, text = "The cursor is now active. Please\nuse the cursor to select the location\nwhich you would like to move the point to.", 
            background = "white"), relx = 0.1, rely = 0.2, `in` = frameLP)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", TargetLocation)
            tkconfigure(img, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", TargetLocation)
            tkconfigure(img2, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", TargetLocation)
            tkconfigure(img3, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", TargetLocation)
            tkconfigure(img4, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", TargetLocation)
            tkconfigure(img5, cursor = "crosshair")
        }
        tkwait.window(MGcomp$MPoint)
        ActiveArrowCursor()
    }
    NewStartConfig <- function() {
        MGvar$CSinitConf <<- "no"
        if (MGvar$GenSet.CalcScree == "yes") {
            MGvar$GenSet.CalcScree <<- "no"
            screetracker = 1
        }
        if (MGvar$MDStype == "ClasScal") {
            tkmessageBox(message = "Classical Scaling does not use a starting configuration")
        }
        if (MGvar$MDStype == "Sammon") {
            SammonMapping(initial = FALSE)
        }
        if (MGvar$MDStype == "Kruskal") {
            MyKruskalsMDS(initial = FALSE)
        }
        if (MGvar$MDStype == "M.Smac.Sym") {
            MyMetricSmacofSym(initial = FALSE)
        }
        if (MGvar$MDStype == "NM.Smac.Sym") {
            MyNonMetricSmacofSym(initial = FALSE)
        }
        if (screetracker == 1) {
            MGvar$GenSet.CalcScree <<- "yes"
            screetrackerr = 0
        }
    }
    angle3Dup <- function() {
        MGvar$sthreeDplot.angle <<- MGvar$sthreeDplot.angle + 
            5
        tkrreplot(img3Dstat)
    }
    angle3Ddown <- function() {
        MGvar$sthreeDplot.angle <<- MGvar$sthreeDplot.angle - 
            5
        tkrreplot(img3Dstat)
    }
    scale3Dup <- function() {
        MGvar$sthreeDplot.yscale <<- MGvar$sthreeDplot.yscale + 
            0.05
        tkrreplot(img3Dstat)
    }
    scale3Ddown <- function() {
        MGvar$sthreeDplot.yscale <<- MGvar$sthreeDplot.yscale - 
            0.05
        tkrreplot(img3Dstat)
    }
    orig3D <- function() {
        MGvar$sthreeDplot.angle <<- 40
        MGvar$sthreeDplot.yscale <<- 1
        tkrreplot(img3Dstat)
    }
    RemovePoint <- function() {
        MGcomp$RmPoint <<- tktoplevel()
        tkwm.title(MGcomp$RmPoint, "Choose Point")
        tkwm.geometry(MGcomp$RmPoint, "280x90")
        RmPointcanvas = tkcanvas(MGcomp$RmPoint, width = 300, 
            height = 120, bg = col.sec)
        tkplace(RmPointcanvas, relx = 0, rely = 0, `in` = MGcomp$RmPoint)
        frameRmp <- tkwidget(MGcomp$RmPoint, "TitleFrame", text = "Active Cursor", 
            background = "white")
        tkplace(frameRmp, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = MGcomp$RmPoint)
        tkplace(tklabel(frameRmp, text = "The cursor is now active. Please\nuse the cursor to select the point\nwhich you would like to remove.", 
            background = "white"), relx = 0.1, rely = 0.2, `in` = frameRmp)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", GetCoordsLeftClick.RemovePoint)
            tkconfigure(img, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", GetCoordsLeftClick.RemovePoint)
            tkconfigure(img2, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", GetCoordsLeftClick.RemovePoint)
            tkconfigure(img3, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", GetCoordsLeftClick.RemovePoint)
            tkconfigure(img4, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", GetCoordsLeftClick.RemovePoint)
            tkconfigure(img5, cursor = "crosshair")
        }
        tkwait.window(MGcomp$RmPoint)
        ActiveArrowCursor()
    }
    GetClosestPointRemovePoint <- function(xClick, yClick, imgXcoords, 
        imgYcoords) {
        squared.Distance <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        Closest <- which.min(squared.Distance)
        compnams = rownames(MGvar$activedata)
        MDSnams = rownames(MGvar$MDSmat)
        MGvar$GenSet.ClearIL <<- "yes"
        ActiveIndex()
        ptm = Closest
        MGvar$remindex <<- MGvar$remindex + 1
        Activeremindex()
        MGvar$RemovedPoints[MGvar$remindex] <<- rownames(MGvar$MDSmat)[ptm]
        MGvar$remPcompindex <<- c()
        for (i in 1:length(MGvar$RemovedPoints)) {
            for (j in 1:length(compnams)) {
                if (MGvar$RemovedPoints[i] == compnams[j]) {
                  MGvar$remPcompindex <<- append(MGvar$remPcompindex, 
                    j)
                }
            }
        }
        ActiveremPcompindex()
        ActiveRemovePoints()
        MGvar$MDSmat <<- MGvar$MDSmat[-ptm, ]
        MGvar$distmat <<- MGvar$distmat[-ptm, -ptm]
        ActiveCoordMat()
        ActivedistMat()
        MGvar$removedpoints <<- TRUE
        ActiveRemPts()
        MGvar$removedpointsactivedata <<- MGvar$removedpointsactivedata[-ptm, 
            ]
        ActiveRemPtsData()
        MGvar$tShepx <<- as.vector(0)
        tabplot()
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
        tableupdate.Remp()
        StressUpdate()
    }
    ActiveremPcompindex <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$remPcompindex.T1 <<- MGvar$remPcompindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$remPcompindex.T2 <<- MGvar$remPcompindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$remPcompindex.T3 <<- MGvar$remPcompindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$remPcompindex.T4 <<- MGvar$remPcompindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$remPcompindex.T5 <<- MGvar$remPcompindex
        }
    }
    Activeremindex <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$remindex.T1 <<- MGvar$remindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$remindex.T2 <<- MGvar$remindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$remindex.T3 <<- MGvar$remindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$remindex.T4 <<- MGvar$remindex
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$remindex.T5 <<- MGvar$remindex
        }
    }
    ActiveRemovePoints <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$RemovedPoints.T1 <<- MGvar$RemovedPoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$RemovedPoints.T2 <<- MGvar$RemovedPoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$RemovedPoints.T3 <<- MGvar$RemovedPoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$RemovedPoints.T4 <<- MGvar$RemovedPoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$RemovedPoints.T5 <<- MGvar$RemovedPoints
        }
    }
    ActiveRemPtsData <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$removedpointsactivedata.T1 <<- MGvar$removedpointsactivedata
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$removedpointsactivedata.T2 <<- MGvar$removedpointsactivedata
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$removedpointsactivedata.T3 <<- MGvar$removedpointsactivedata
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$removedpointsactivedata.T4 <<- MGvar$removedpointsactivedata
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$removedpointsactivedata.T5 <<- MGvar$removedpointsactivedata
        }
    }
    ActiveRemPts <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$removedpoints.T1 <<- MGvar$removedpoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$removedpoints.T2 <<- MGvar$removedpoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$removedpoints.T3 <<- MGvar$removedpoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$removedpoints.T4 <<- MGvar$removedpoints
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$removedpoints.T5 <<- MGvar$removedpoints
        }
    }
    GetCoordsLeftClick.RemovePoint <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$MDSmat[, 1] - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$MDSmat[, 2] - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        xPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        yPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        GetClosestPointRemovePoint(xClick, yClick, imgXcoords, 
            imgYcoords)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", OnPlotLeftClick)
        }
        tkdestroy(MGcomp$RmPoint)
    }
    freshstart <- function() {
        MGvar$removedpoints <<- FALSE
        MGvar$removedpointsactivedata <<- MGvar$activedata
        ActiveRemPtsData()
    }
    BrushingShep <- function(x, y) {
        if (MGvar$GenSet.CalcShep == "yes") {
            MGvar$Shep.indexLabeled <<- c()
            xClick <- x
            yClick <- y
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                imgshep)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                imgshep)))
            xMin <- MGvar$Shep.parPlotSize[1] * width
            xMax <- MGvar$Shep.parPlotSize[2] * width
            yMin <- MGvar$Shep.parPlotSize[3] * height
            yMax <- MGvar$Shep.parPlotSize[4] * height
            rangeX <- MGvar$Shep.usrCoords[2] - MGvar$Shep.usrCoords[1]
            rangeY <- MGvar$Shep.usrCoords[4] - MGvar$Shep.usrCoords[3]
            imgXcoords <- (MGvar$Shepx - MGvar$Shep.usrCoords[1]) * 
                (xMax - xMin)/rangeX + xMin
            imgYcoords <- (MGvar$Shepy - MGvar$Shep.usrCoords[3]) * 
                (yMax - yMin)/rangeY + yMin
            xClick <- as.numeric(xClick) + 0.5
            yClick <- as.numeric(yClick) + 0.5
            yClick <- height - yClick
            MGvar$latest.xPlotCoord <<- MGvar$Shep.usrCoords[1] + 
                (xClick - xMin) * rangeX/(xMax - xMin)
            MGvar$latest.yPlotCoord <<- MGvar$Shep.usrCoords[3] + 
                (yClick - yMin) * rangeY/(yMax - yMin)
            if (MGvar$shep.initpos) {
                MGvar$first.xPlotCoord <<- MGvar$latest.xPlotCoord
                MGvar$first.yPlotCoord <<- MGvar$latest.yPlotCoord
                MGvar$shep.initpos <<- FALSE
            }
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
        }
    }
    OnShepRelease <- function() {
        if (MGvar$GenSet.CalcShep == "yes") {
            count = length(MGvar$Shep.indexLabeled) + 1
            for (i in 1:length(MGvar$Shepx)) {
                if (min(MGvar$first.xPlotCoord, MGvar$latest.xPlotCoord) < 
                  MGvar$Shepx[i] && max(MGvar$first.xPlotCoord, 
                  MGvar$latest.xPlotCoord) > MGvar$Shepx[i] && 
                  min(MGvar$first.yPlotCoord, MGvar$latest.yPlotCoord) < 
                    MGvar$Shepy[i] && max(MGvar$first.yPlotCoord, 
                  MGvar$latest.yPlotCoord) > MGvar$Shepy[i]) {
                  MGvar$Shep.indexLabeled <<- c(MGvar$Shep.indexLabeled, 
                    i)
                  count = count + 1
                }
            }
            tabplot()
            MGvar$latest.xPlotCoord <<- 0
            MGvar$latest.yPlotCoord <<- 0
            MGvar$first.xPlotCoord <<- 0
            MGvar$first.yPlotCoord <<- 0
            MGvar$shep.initpos <<- TRUE
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
        }
    }
    BrushingPointMove <- function(x, y) {
        if (MGvar$main.initpos) {
            Enableundo()
        }
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$X.Coords - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$Y.Coords - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        MGvar$latest.xPlotCoord <- MGvar$usrCoords[1] + (xClick - 
            xMin) * rangeX/(xMax - xMin)
        MGvar$latest.yPlotCoord <- MGvar$usrCoords[3] + (yClick - 
            yMin) * rangeY/(yMax - yMin)
        if (MGvar$main.initpos) {
            MGvar$main.initpos <<- FALSE
            GetActivePoint(xClick, yClick, imgXcoords, imgYcoords)
        }
        MGvar$MDSmat[MGvar$ActivePoint, 1] <<- MGvar$latest.xPlotCoord
        MGvar$MDSmat[MGvar$ActivePoint, 2] <<- MGvar$latest.yPlotCoord
        ActiveCoordMat()
        tabplot()
        StressUpdate()
        if (MGvar$GenSet.UpShep == "yes" && MGvar$GenSet.CalcShep == 
            "yes") {
            tkrreplot(imgshep)
        }
        MGvar$movingpoint <<- TRUE
    }
    GetActivePoint <- function(xClick, yClick, imgXcoords, imgYcoords) {
        squared.Distance.P <- (xClick - imgXcoords)^2 + (yClick - 
            imgYcoords)^2
        MGvar$ActivePoint <<- which.min(squared.Distance.P)
    }
    OnRelease.Main <- function() {
        if (MGvar$movingpoint) {
            lIL = length(MGvar$indexLabeled)
            MGvar$indexLabeled <<- MGvar$indexLabeled[-lIL]
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$indexLabeled.T1 <<- MGvar$indexLabeled
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$indexLabeled.T2 <<- MGvar$indexLabeled
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$indexLabeled.T3 <<- MGvar$indexLabeled
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$indexLabeled.T4 <<- MGvar$indexLabeled
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$indexLabeled.T5 <<- MGvar$indexLabeled
            }
            MGvar$ActivePoint <<- 0
        }
        MGvar$movingpoint <<- FALSE
        tabplot()
        MGvar$main.initpos <<- TRUE
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
    }
    UndoMove <- function() {
        if (nrow(MGvar$prevMat) > 1) {
            MGvar$MDSmat <<- MGvar$prevMat
            ActiveCoordMat()
            tabplot()
            StressUpdate()
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
        }
        ResetUndo()
    }
    Enableundo <- function() {
        MGvar$prevMat <<- MGvar$MDSmat
        tkentryconfigure(GeneralMenu, 0, state = "active")
    }
    ResetUndo <- function() {
        MGvar$prevMat <<- as.matrix(0)
        tkentryconfigure(GeneralMenu, 0, state = "disabled")
    }
    BrushingMainPlot <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$X.Coords - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$Y.Coords - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        MGvar$latest.xCoord <<- MGvar$usrCoords[1] + (xClick - 
            xMin) * rangeX/(xMax - xMin)
        MGvar$latest.yCoord <<- MGvar$usrCoords[3] + (yClick - 
            yMin) * rangeY/(yMax - yMin)
        if (MGvar$main.initpos) {
            MGvar$first.xCoord <<- MGvar$latest.xCoord
            MGvar$first.yCoord <<- MGvar$latest.yCoord
            MGvar$main.initpos <<- FALSE
        }
        tabplot()
    }
    OnRelease.Main.Col <- function() {
        lIL = length(MGvar$indexLabeled)
        MGvar$indexLabeled <<- MGvar$indexLabeled[-lIL]
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$indexLabeled.T1 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$indexLabeled.T2 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$indexLabeled.T3 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$indexLabeled.T4 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$indexLabeled.T5 <<- MGvar$indexLabeled
        }
        PCol <- tktoplevel()
        tkwm.resizable(PCol, "0", "0")
        tkwm.title(PCol, "Point Colour")
        tkwm.geometry(PCol, "250x90")
        PColcanvas = tkcanvas(PCol, width = 300, height = 120, 
            bg = col.sec)
        tkplace(PColcanvas, relx = 0, rely = 0, `in` = PCol)
        framePC <- tkwidget(PCol, "TitleFrame", text = "Choose Point Colour", 
            background = "white")
        tkplace(framePC, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = PCol)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            ptcol <- MGvar$activeplot.pointcol.T1
            ptcol.temp <- MGvar$activeplot.pointcol.T1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            ptcol <- MGvar$activeplot.pointcol.T2
            ptcol.temp <- MGvar$activeplot.pointcol.T2
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            ptcol <- MGvar$activeplot.pointcol.T3
            ptcol.temp <- MGvar$activeplot.pointcol.T3
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            ptcol <- MGvar$activeplot.pointcol.T4
            ptcol.temp <- MGvar$activeplot.pointcol.T4
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            ptcol <- MGvar$activeplot.pointcol.T5
            ptcol.temp <- MGvar$activeplot.pointcol.T5
        }
        tkplace(tklabel(framePC, text = "Selected Points", background = "white"), 
            relx = 0.1, rely = 0.35, `in` = framePC)
        ChangePtCol <- function() {
            ptcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = ptcol.temp, title = "Choose a Colour"))))
            if (nchar(ptcol.temp) > 0) 
                tkconfigure(Ptcol.but, bg = ptcol.temp)
            ptcol <<- ptcol.temp
            tkdestroy(PCol)
        }
        Ptcol.but <- tkbutton(PCol, text = "", width = 2, height = 1, 
            bg = ptcol, command = function() ChangePtCol())
        tkplace(Ptcol.but, relx = 0.8, rely = 0.35, `in` = framePC)
        tkwait.window(PCol)
        colvec = c()
        for (i in 1:nrow(MGvar$MDSmat)) {
            if (min(MGvar$first.xCoord, MGvar$latest.xCoord) < 
                MGvar$MDSmat[i, 1] && max(MGvar$first.xCoord, 
                MGvar$latest.xCoord) > MGvar$MDSmat[i, 1] && 
                min(MGvar$first.yCoord, MGvar$latest.yCoord) < 
                  MGvar$MDSmat[i, 2] && max(MGvar$first.yCoord, 
                MGvar$latest.yCoord) > MGvar$MDSmat[i, 2]) {
                for (j in 1:nrow(MGvar$originaldistmat)) {
                  if (rownames(MGvar$MDSmat)[i] == rownames(MGvar$originaldistmat)[j]) {
                    colvec = c(colvec, j)
                  }
                }
            }
        }
        for (j in 1:length(colvec)) {
            MGvar$MDSmat.Cols[colvec[j]] <<- ptcol
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$MDSmat.Cols.T1[colvec[j]] <<- ptcol
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$MDSmat.Cols.T2[colvec[j]] <<- ptcol
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$MDSmat.Cols.T3[colvec[j]] <<- ptcol
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$MDSmat.Cols.T4[colvec[j]] <<- ptcol
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$MDSmat.Cols.T5[colvec[j]] <<- ptcol
            }
        }
        MGvar$main.initpos <<- TRUE
        MGvar$first.xCoord <<- 0
        MGvar$latest.xCoord <<- 0
        MGvar$first.yCoord <<- 0
        MGvar$latest.yCoord <<- 0
        tabplot()
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<B1-Motion>", BrushingPointMove)
            tkbind(img, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<B1-Motion>", BrushingPointMove)
            tkbind(img2, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<B1-Motion>", BrushingPointMove)
            tkbind(img3, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<B1-Motion>", BrushingPointMove)
            tkbind(img4, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<B1-Motion>", BrushingPointMove)
            tkbind(img5, "<ButtonRelease-1>", OnRelease.Main)
        }
    }
    ColButtonInitiate <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<B1-Motion>", BrushingMainPlot)
            tkbind(img, "<ButtonRelease-1>", OnRelease.Main.Col)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<B1-Motion>", BrushingMainPlot)
            tkbind(img2, "<ButtonRelease-1>", OnRelease.Main.Col)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<B1-Motion>", BrushingMainPlot)
            tkbind(img3, "<ButtonRelease-1>", OnRelease.Main.Col)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<B1-Motion>", BrushingMainPlot)
            tkbind(img4, "<ButtonRelease-1>", OnRelease.Main.Col)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<B1-Motion>", BrushingMainPlot)
            tkbind(img5, "<ButtonRelease-1>", OnRelease.Main.Col)
        }
    }
    OnRelease.Main.Move <- function() {
        lIL = length(MGvar$indexLabeled)
        MGvar$indexLabeled <<- MGvar$indexLabeled[-lIL]
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$indexLabeled.T1 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$indexLabeled.T2 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$indexLabeled.T3 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$indexLabeled.T4 <<- MGvar$indexLabeled
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$indexLabeled.T5 <<- MGvar$indexLabeled
        }
        midX = max(MGvar$first.xCoord, MGvar$latest.xCoord) - 
            (max(MGvar$first.xCoord, MGvar$latest.xCoord) - min(MGvar$first.xCoord, 
                MGvar$latest.xCoord))/2
        midY = max(MGvar$first.yCoord, MGvar$latest.yCoord) - 
            (max(MGvar$first.yCoord, MGvar$latest.yCoord) - min(MGvar$first.yCoord, 
                MGvar$latest.yCoord))/2
        MGvar$relIndex <<- as.vector(0)
        MGvar$relXCoords <<- as.vector(0)
        MGvar$relYCoords <<- as.vector(0)
        count = 1
        for (i in 1:nrow(MGvar$MDSmat)) {
            if (min(MGvar$first.xCoord, MGvar$latest.xCoord) < 
                MGvar$MDSmat[i, 1] && max(MGvar$first.xCoord, 
                MGvar$latest.xCoord) > MGvar$MDSmat[i, 1] && 
                min(MGvar$first.yCoord, MGvar$latest.yCoord) < 
                  MGvar$MDSmat[i, 2] && max(MGvar$first.yCoord, 
                MGvar$latest.yCoord) > MGvar$MDSmat[i, 2]) {
                MGvar$relIndex[count] <<- i
                MGvar$relXCoords[count] <<- MGvar$MDSmat[i, 1] - 
                  midX
                MGvar$relYCoords[count] <<- MGvar$MDSmat[i, 2] - 
                  midY
                count = count + 1
            }
        }
        MGvar$main.initpos <<- TRUE
        MGvar$first.xCoord <<- 0
        MGvar$latest.xCoord <<- 0
        MGvar$first.yCoord <<- 0
        MGvar$latest.yCoord <<- 0
        MoveGroup()
    }
    MoveButtonInitiate <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkconfigure(img, cursor = "crosshair")
            tkbind(img, "<B1-Motion>", BrushingMainPlot)
            tkbind(img, "<ButtonRelease-1>", OnRelease.Main.Move)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkconfigure(img2, cursor = "crosshair")
            tkbind(img2, "<B1-Motion>", BrushingMainPlot)
            tkbind(img2, "<ButtonRelease-1>", OnRelease.Main.Move)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkconfigure(img3, cursor = "crosshair")
            tkbind(img3, "<B1-Motion>", BrushingMainPlot)
            tkbind(img3, "<ButtonRelease-1>", OnRelease.Main.Move)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkconfigure(img4, cursor = "crosshair")
            tkbind(img4, "<B1-Motion>", BrushingMainPlot)
            tkbind(img4, "<ButtonRelease-1>", OnRelease.Main.Move)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkconfigure(img5, cursor = "crosshair")
            tkbind(img5, "<B1-Motion>", BrushingMainPlot)
            tkbind(img5, "<ButtonRelease-1>", OnRelease.Main.Move)
        }
    }
    MoveGroupFunction <- function(TargetX, TargetY) {
        MGvar$prevMat <<- MGvar$MDSmat
        for (i in 1:length(MGvar$relIndex)) {
            MGvar$MDSmat[MGvar$relIndex[i], 1] <<- TargetX + 
                MGvar$relXCoords[i]
            MGvar$MDSmat[MGvar$relIndex[i], 2] <<- TargetY + 
                MGvar$relYCoords[i]
        }
        ActiveCoordMat()
        tkdestroy(MGcomp$MGroup)
        tabplot()
        if (MGvar$GenSet.CalcShep == "yes") {
            tkrreplot(imgshep)
        }
        StressUpdate()
    }
    GroupMoveTargetLocation <- function(x, y) {
        xClick <- x
        yClick <- y
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img2)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img2)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img3)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img3)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img4)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img4)))
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            width <- as.numeric(tclvalue(tkwinfo("reqwidth", 
                img5)))
            height <- as.numeric(tclvalue(tkwinfo("reqheight", 
                img5)))
        }
        xMin <- MGvar$parPlotSize[1] * width
        xMax <- MGvar$parPlotSize[2] * width
        yMin <- MGvar$parPlotSize[3] * height
        yMax <- MGvar$parPlotSize[4] * height
        rangeX <- MGvar$usrCoords[2] - MGvar$usrCoords[1]
        rangeY <- MGvar$usrCoords[4] - MGvar$usrCoords[3]
        imgXcoords <- (MGvar$X.Coords - MGvar$usrCoords[1]) * 
            (xMax - xMin)/rangeX + xMin
        imgYcoords <- (MGvar$Y.Coords - MGvar$usrCoords[3]) * 
            (yMax - yMin)/rangeY + yMin
        xClick <- as.numeric(xClick) + 0.5
        yClick <- as.numeric(yClick) + 0.5
        yClick <- height - yClick
        TargetxPlotCoord <- MGvar$usrCoords[1] + (xClick - xMin) * 
            rangeX/(xMax - xMin)
        TargetyPlotCoord <- MGvar$usrCoords[3] + (yClick - yMin) * 
            rangeY/(yMax - yMin)
        MoveGroupFunction(TargetxPlotCoord, TargetyPlotCoord)
    }
    MoveGroup <- function() {
        MGcomp$MGroup <<- tktoplevel()
        tkwm.title(MGcomp$MGroup, "Choose Point")
        tkwm.geometry(MGcomp$MGroup, "280x90")
        MGroupcanvas = tkcanvas(MGcomp$MGroup, width = 300, height = 120, 
            bg = col.sec)
        tkplace(MGroupcanvas, relx = 0, rely = 0, `in` = MGcomp$MGroup)
        frameMG <- tkwidget(MGcomp$MGroup, "TitleFrame", text = "Active Cursor", 
            background = "white")
        tkplace(frameMG, relx = 0.05, rely = 0.05, relheight = 0.9, 
            relwidth = 0.9, `in` = MGcomp$MGroup)
        tkplace(tklabel(frameMG, text = "The cursor is now active. Please\nuse the cursor to select the central\npoint of the groups' location.", 
            background = "white"), relx = 0.1, rely = 0.2, `in` = frameMG)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", GroupMoveTargetLocation)
            tkconfigure(img, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", GroupMoveTargetLocation)
            tkconfigure(img2, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", GroupMoveTargetLocation)
            tkconfigure(img3, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", GroupMoveTargetLocation)
            tkconfigure(img4, cursor = "crosshair")
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", GroupMoveTargetLocation)
            tkconfigure(img5, cursor = "crosshair")
        }
        tkwait.window(MGcomp$MGroup)
        ActiveArrowCursor()
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<Button-1>", OnPlotLeftClick)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tkbind(img, "<B1-Motion>", BrushingPointMove)
            tkbind(img, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tkbind(img2, "<B1-Motion>", BrushingPointMove)
            tkbind(img2, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tkbind(img3, "<B1-Motion>", BrushingPointMove)
            tkbind(img3, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tkbind(img4, "<B1-Motion>", BrushingPointMove)
            tkbind(img4, "<ButtonRelease-1>", OnRelease.Main)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tkbind(img5, "<B1-Motion>", BrushingPointMove)
            tkbind(img5, "<ButtonRelease-1>", OnRelease.Main)
        }
    }
    ColoursDataOption <- function() {
        MGvar$MDSmat.Cols <<- as.matrix(MGvar$MDSmat.Cols)
        rownames(MGvar$MDSmat.Cols) <<- rownames(MGvar$activedata)
        MGvar$MDSmat.Cols <<- my.tk2edit(MGvar$MDSmat.Cols)
        MGvar$MDSmat.Cols <<- as.vector(MGvar$MDSmat.Cols)
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
        }
        tabplot()
    }
    plottingTabLeftClick <- function() {
        AcTab = tk2notetab.text(myPlottingNB)
        if (AcTab == "Plot1") {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab1"
            ActiveTabChanges()
        }
        if (AcTab == "Plot2") {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab2"
            ActiveTabChanges()
        }
        if (AcTab == "Plot3") {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab3"
            ActiveTabChanges()
        }
        if (AcTab == "Plot4") {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab4"
            ActiveTabChanges()
        }
        if (AcTab == "Plot5") {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab5"
            ActiveTabChanges()
        }
    }
    StressCalculation <- function(MDSmat, distmat, isMet = TRUE) {
        Calc <- Shepard(as.dist(distmat), MDSmat)
        proxdist <- Calc$x
        distfunc <- Calc$yf
        MDSdist <- Calc$y
        STRESS1 <- function(dhat, d) {
            S1 = sqrt(sum((d - dhat)^2)/sum(d^2))
            S1 = min(S1, 1)
            return(S1)
        }
        STRESS2 <- function(dhat, d) {
            S2 = sqrt((sum((d - dhat)^2))/(sum((d - mean(d))^2)))
            S2 = min(S2, 1)
            return(S2)
        }
        NormRawSTRESS <- function(dhat, d) {
            NRS <- (STRESS1(dhat, d))^2
            NRS = min(NRS, 1)
            return(NRS)
        }
        CorCoef <- function(dhat, d) {
            V1 = as.matrix(as.vector(d))
            V2 = as.matrix(as.vector(dhat))
            CC = corr(cbind(V1, V2))
            return(CC)
        }
        if (isMet) {
            if (MGvar$StCalc == "NormRawStress") {
                Stress = NormRawSTRESS(proxdist, MDSdist)
            }
            if (MGvar$StCalc == "Stress-1") {
                Stress = STRESS1(proxdist, MDSdist)
            }
            if (MGvar$StCalc == "Stress-2") {
                Stress = STRESS2(proxdist, MDSdist)
            }
            if (MGvar$StCalc == "CorCoef") {
                Stress = CorCoef(proxdist, MDSdist)
            }
        }
        if (!isMet) {
            if (MGvar$StCalc == "NormRawStress") {
                Stress = NormRawSTRESS(distfunc, MDSdist)
            }
            if (MGvar$StCalc == "Stress-1") {
                Stress = STRESS1(distfunc, MDSdist)
            }
            if (MGvar$StCalc == "Stress-2") {
                Stress = STRESS2(distfunc, MDSdist)
            }
            if (MGvar$StCalc == "CorCoef") {
                Stress = CorCoef(distfunc, MDSdist)
            }
        }
        return(Stress)
    }
    StressUpdate <- function() {
        if (nrow(MGvar$MDSmat.T1) != 1 && ncol(MGvar$MDSmat.T1) != 
            1) {
            if (MGvar$MDStype.T1 != "INDSCAL" && MGvar$MDStype.T1 != 
                "Gifi") {
                MGvar$MDSStress.T1 <<- round(StressCalculation(MGvar$MDSmat.T1, 
                  MGvar$distmat.T1, MGvar$is.Metric.T1), 3)
            }
            else {
                MGvar$MDSStress.T1 <<- paste("*")
            }
        }
        if (nrow(MGvar$MDSmat.T2) != 1 && ncol(MGvar$MDSmat.T2) != 
            1) {
            if (MGvar$MDStype.T2 != "INDSCAL" && MGvar$MDStype.T2 != 
                "Gifi") {
                MGvar$MDSStress.T2 <<- round(StressCalculation(MGvar$MDSmat.T2, 
                  MGvar$distmat.T2, MGvar$is.Metric.T2), 3)
            }
            else {
                MGvar$MDSStress.T2 <<- paste("*")
            }
        }
        if (nrow(MGvar$MDSmat.T3) != 1 && ncol(MGvar$MDSmat.T3) != 
            1) {
            if (MGvar$MDStype.T3 != "INDSCAL" && MGvar$MDStype.T3 != 
                "Gifi") {
                MGvar$MDSStress.T3 <<- round(StressCalculation(MGvar$MDSmat.T3, 
                  MGvar$distmat.T3, MGvar$is.Metric.T3), 3)
            }
            else {
                MGvar$MDSStress.T3 <<- paste("*")
            }
        }
        if (nrow(MGvar$MDSmat.T4) != 1 && ncol(MGvar$MDSmat.T4) != 
            1) {
            if (MGvar$MDStype.T4 != "INDSCAL" && MGvar$MDStype.T4 != 
                "Gifi") {
                MGvar$MDSStress.T4 <<- round(StressCalculation(MGvar$MDSmat.T4, 
                  MGvar$distmat.T4, MGvar$is.Metric.T4), 3)
            }
            else {
                MGvar$MDSStress.T4 <<- paste("*")
            }
        }
        if (nrow(MGvar$MDSmat.T5) != 1 && ncol(MGvar$MDSmat.T5) != 
            1) {
            if (MGvar$MDStype.T5 != "INDSCAL" && MGvar$MDStype.T5 != 
                "Gifi") {
                MGvar$MDSStress.T5 <<- round(StressCalculation(MGvar$MDSmat.T5, 
                  MGvar$distmat.T5, MGvar$is.Metric.T5), 3)
            }
            else {
                MGvar$MDSStress.T5 <<- paste("*")
            }
        }
        tableupdate()
    }
    CatCols <- function() {
        CCols = tktoplevel()
        tkwm.resizable(CCols, "0", "0")
        tkwm.deiconify(CCols)
        tkwm.title(CCols, "Catagory Colours")
        tkwm.geometry(CCols, "320x150")
        CColscanvas = tkcanvas(CCols, bg = col.sec)
        tkplace(CColscanvas, relx = 0, rely = 0, relheight = 1, 
            relwidth = 1, `in` = CCols)
        frameCC <- tkwidget(CCols, "TitleFrame", text = "Catagory Colours Options", 
            background = "white")
        tkplace(frameCC, relx = 0.01, rely = 0.01, relwidth = 0.98, 
            relheight = 0.98, `in` = CCols)
        fontsmall <- tkfont.create(family = "times", size = 9)
        tkplace(tklabel(frameCC, text = "If a column of object catagories was included in the uploaded\ndata, the related catagory colours can be adjusted here.", 
            font = fontsmall), relx = 0.03, rely = 0.15, `in` = frameCC)
        Cats <- names(MGvar$ClasTab)
        CatVar <- tclVar(Cats[1])
        CatVar.ComboBox <- tkwidget(CCols, "ComboBox", editable = FALSE, 
            values = Cats, width = 15, textvariable = CatVar, 
            command = function() ColChanger())
        tkplace(CatVar.ComboBox, relx = 0.15, rely = 0.48, `in` = frameCC)
        catcol = MGvar$ClasTabCols[1]
        catcol.temp = MGvar$ClasTabCols[1]
        ColChanger <- function() {
            CurCat <- as.character(tclvalue(CatVar))
            for (i in 1:length(MGvar$ClasTabCols)) {
                if (CurCat == names(MGvar$ClasTabCols)[i]) {
                  tkconfigure(CatColBut, bg = MGvar$ClasTabCols[i])
                  break
                }
            }
        }
        ChangeCatCol <- function() {
            catcol.temp <<- tclvalue(.Tcl(paste("tk_chooseColor", 
                .Tcl.args(initialcolor = catcol.temp, title = "Choose a colour"))))
            if (nchar(catcol.temp) > 0) 
                tkconfigure(CatColBut, bg = catcol.temp)
        }
        CatColBut <- tkbutton(CCols, text = "", width = 2, height = 1, 
            bg = catcol, command = function() ChangeCatCol())
        tkplace(CatColBut, relx = 0.7, rely = 0.47, `in` = frameCC)
        ChangeCC <- function() {
            catcol <<- catcol.temp
            CatV = as.character(tclvalue(CatVar))
            for (i in 1:length(MGvar$ClasTabCols)) {
                if (CatV == names(MGvar$ClasTabCols)[i]) {
                  MGvar$ClasTabCols[i] <<- catcol
                }
            }
            for (i in 1:nrow(MGvar$activedata)) {
                if (MGvar$ClasVec[i] == CatV) {
                  MGvar$MDSmat.Cols[i] <<- catcol
                }
            }
            MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
            MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
            MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
            MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
            MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
            tabplot()
        }
        tkplace(tkbutton(CCols, text = "Change", width = 15, 
            command = function() ChangeCC()), relx = 0.325, rely = 0.75, 
            `in` = frameCC)
        tkbind(CatVar.ComboBox, "<Leave>", ColChanger)
        tkbind(CatVar.ComboBox, "<Enter>", ColChanger)
        tkbind(CatColBut, "<Enter>", ColChanger)
        tkbind(CatColBut, "<Enter>", ColChanger)
        tkbind(CCols, "<Enter>", ColChanger)
    }
    tableplaceRemP <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tabtitle = "Plot1"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tabtitle = "Plot2"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tabtitle = "Plot3"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tabtitle = "Plot4"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tabtitle = "Plot5"
        }
        Remp = MGvar$RemovedPoints
        myRemParray <- c(tabtitle, Remp[1], Remp[2], Remp[3], 
            Remp[4], Remp[5], Remp[6], Remp[7], "-", Remp[8], 
            Remp[9], Remp[10], Remp[11], Remp[12], Remp[13], 
            Remp[14], "-", Remp[15], Remp[16], Remp[17], Remp[18], 
            Remp[19], Remp[20], Remp[21], "-", Remp[22], Remp[23], 
            Remp[24], Remp[25], Remp[26], Remp[27], Remp[28], 
            "-", Remp[29], Remp[30], Remp[31], Remp[32], Remp[33], 
            Remp[34], Remp[35], "-", Remp[36], Remp[37], Remp[38], 
            Remp[39], Remp[40], Remp[41], Remp[42])
        MGcomp$Remparray <<- tclArray()
        dim(myRemParray) <- c(8, 6)
        for (i in 0:7) for (j in 0:5) {
            MGcomp$Remparray[[i, j]] <- myRemParray[i + 1, j + 
                1]
        }
        table.Remp <- tkwidget(frameRemPoints, "table", rows = 9, 
            cols = 6, variable = MGcomp$Remparray, titlerows = "1", 
            background = "white", width = 6, height = 9)
        tcl(table.Remp, "width", 0, 12)
        tcl(table.Remp, "width", 1, 12)
        tcl(table.Remp, "width", 2, 12)
        tcl(table.Remp, "width", 3, 12)
        tcl(table.Remp, "width", 4, 11)
        tcl(table.Remp, "width", 5, 11)
        tabpointfont <- tkfont.create(family = "times", size = 10, 
            weight = "bold")
        count = 0
        for (i in 1:42) {
            if (Remp[i] != "-") {
                count = count + 1
            }
        }
        if (count > 0) {
            for (i in 1:count) {
                if (i%%7 == 0) {
                  r = 7
                }
                if (i%%7 != 0) {
                  r = i%%7
                }
                if (i < 8) {
                  RaC = paste(r, ",0", sep = "")
                  cind = 0
                }
                if (i > 7 && i < 15) {
                  RaC = paste(r, ",1", sep = "")
                  cind = 1
                }
                if (i > 14 && i < 22) {
                  RaC = paste(r, ",2", sep = "")
                  cind = 2
                }
                if (i > 21 && i < 29) {
                  RaC = paste(r, ",3", sep = "")
                  cind = 3
                }
                if (i > 28 && i < 36) {
                  RaC = paste(r, ",4", sep = "")
                  cind = 4
                }
                if (i > 35 && i < 43) {
                  RaC = paste(r, ",5", sep = "")
                  cind = 5
                }
                for (j in 1:nrow(MGvar$activedata)) {
                  if (Remp[i] == rownames(MGvar$activedata)[j]) {
                    Col = MGvar$MDSmat.Cols[j]
                  }
                }
                tcl(.Tk.ID(table.Remp), "tag", "celltag", i + 
                  cind, RaC)
                tcl(.Tk.ID(table.Remp), "tag", "configure", i + 
                  cind, fg = Col, font = tabpointfont)
            }
        }
        tkgrid(table.Remp)
        MGcomp$table.Remp <<- table.Remp
    }
    tableupdate.Remp <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tabtitle = "Plot1"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tabtitle = "Plot2"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tabtitle = "Plot3"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tabtitle = "Plot4"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tabtitle = "Plot5"
        }
        Remp = MGvar$RemovedPoints
        myRemParray <- c(tabtitle, Remp[1], Remp[2], Remp[3], 
            Remp[4], Remp[5], Remp[6], Remp[7], "-", Remp[8], 
            Remp[9], Remp[10], Remp[11], Remp[12], Remp[13], 
            Remp[14], "-", Remp[15], Remp[16], Remp[17], Remp[18], 
            Remp[19], Remp[20], Remp[21], "-", Remp[22], Remp[23], 
            Remp[24], Remp[25], Remp[26], Remp[27], Remp[28], 
            "-", Remp[29], Remp[30], Remp[31], Remp[32], Remp[33], 
            Remp[34], Remp[35], "-", Remp[36], Remp[37], Remp[38], 
            Remp[39], Remp[40], Remp[41], Remp[42])
        MGcomp$Remparray <<- tclArray()
        dim(myRemParray) <- c(8, 6)
        for (i in 0:7) for (j in 0:5) {
            MGcomp$Remparray[[i, j]] <- myRemParray[i + 1, j + 
                1]
        }
        tkconfigure(MGcomp$table.Remp, variable = MGcomp$Remparray)
        tabpointfont <- tkfont.create(family = "times", size = 10, 
            weight = "bold")
        count = 0
        for (i in 1:42) {
            if (Remp[i] != "-") {
                count = count + 1
            }
        }
        if (count > 0) {
            for (i in 1:count) {
                if (i%%7 == 0) {
                  r = 7
                }
                if (i%%7 != 0) {
                  r = i%%7
                }
                if (i < 8) {
                  RaC = paste(r, ",0", sep = "")
                  cind = 0
                }
                if (i > 7 && i < 15) {
                  RaC = paste(r, ",1", sep = "")
                  cind = 1
                }
                if (i > 14 && i < 22) {
                  RaC = paste(r, ",2", sep = "")
                  cind = 2
                }
                if (i > 21 && i < 29) {
                  RaC = paste(r, ",3", sep = "")
                  cind = 3
                }
                if (i > 28 && i < 36) {
                  RaC = paste(r, ",4", sep = "")
                  cind = 4
                }
                if (i > 35 && i < 43) {
                  RaC = paste(r, ",5", sep = "")
                  cind = 5
                }
                for (j in 1:nrow(MGvar$activedata)) {
                  if (Remp[i] == rownames(MGvar$activedata)[j]) {
                    Col = MGvar$MDSmat.Cols[j]
                  }
                }
                tcl(.Tk.ID(MGcomp$table.Remp), "tag", "celltag", 
                  i + cind, RaC)
                tcl(.Tk.ID(MGcomp$table.Remp), "tag", "configure", 
                  i + cind, fg = Col, font = tabpointfont)
            }
        }
    }
    PointReplace <- function() {
        activeRow <- as.numeric(tkindex(MGcomp$table.Remp, "active", 
            "row"))
        activeCol <- as.numeric(tkindex(MGcomp$table.Remp, "active", 
            "col"))
        tabind <- 1 + 7 * activeCol + (activeRow - 1)
        if (MGvar$RemovedPoints[tabind] == "-") {
            tkmessageBox(message = "The Active Cell is empty! Please select an occupied cell.")
        }
        if (MGvar$RemovedPoints[tabind] != "-") {
            MGvar$RemovedPoints <<- MGvar$RemovedPoints[-tabind]
            MGvar$RemovedPoints <<- append(MGvar$RemovedPoints, 
                "-")
            ActiveRemovePoints()
            MGvar$remindex <<- MGvar$remindex - 1
            Activeremindex()
            for (i in 1:nrow(MGvar$activedata)) {
                if (MGvar$remPcompindex[tabind] == i) {
                  ptm = i
                }
            }
            MGvar$remPcompindex <<- MGvar$remPcompindex[-tabind]
            ActiveremPcompindex()
            if (length(MGvar$remPcompindex) > 0) {
                MGvar$distmat <<- MGvar$originaldistmat[-MGvar$remPcompindex, 
                  -MGvar$remPcompindex]
            }
            if (length(MGvar$remPcompindex) == 0) {
                MGvar$distmat <<- MGvar$originaldistmat
            }
            ActivedistMat()
            ExistPoints <- c()
            for (i in 1:nrow(MGvar$activedata)) {
                for (j in 1:nrow(MGvar$MDSmat)) if (rownames(MGvar$activedata)[i] == 
                  rownames(MGvar$MDSmat)[j]) {
                  ExistPoints = c(ExistPoints, i)
                }
            }
            newMDSmat <- matrix(0, nrow = nrow(MGvar$MDSmat) + 
                1, ncol = 2)
            newMDSmat <- as.matrix(newMDSmat)
            if (length(MGvar$remPcompindex) > 0) {
                rownames(newMDSmat) = rownames(MGvar$activedata)[-MGvar$remPcompindex]
            }
            if (length(MGvar$remPcompindex) == 0) {
                rownames(newMDSmat) = rownames(MGvar$activedata)
            }
            counter = 1
            for (i in 1:nrow(newMDSmat)) {
                indic = TRUE
                for (j in 1:nrow(MGvar$MDSmat)) {
                  if (rownames(newMDSmat)[i] == rownames(MGvar$MDSmat)[j]) {
                    newMDSmat[i, ] = MGvar$MDSmat[j, ]
                    indic = FALSE
                  }
                }
                if (indic) {
                  Xrange = max(MGvar$MDSmat[, 1]) - min(MGvar$MDSmat[, 
                    1])
                  newMDSmat[i, 1] = min(MGvar$MDSmat[, 1]) - 
                    0.05 * Xrange
                  newMDSmat[i, 2] = max(MGvar$MDSmat[, 1])
                }
            }
            MGvar$MDSmat <<- newMDSmat
            ActiveCoordMat()
            MGvar$tShepx <<- as.vector(0)
            tabplot()
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            tableupdate.Remp()
            StressUpdate()
        }
    }
    tableplaceRemAx <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tabtitle = "Plot1"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tabtitle = "Plot2"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tabtitle = "Plot3"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tabtitle = "Plot4"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tabtitle = "Plot5"
        }
        RemAx = MGvar$RemovedAxes
        myRemAxarray <- c(tabtitle, RemAx[1], RemAx[2], RemAx[3], 
            RemAx[4], RemAx[5], RemAx[6], RemAx[7], "-", RemAx[8], 
            RemAx[9], RemAx[10], RemAx[11], RemAx[12], RemAx[13], 
            RemAx[14], "-", RemAx[15], RemAx[16], RemAx[17], 
            RemAx[18], RemAx[19], RemAx[20], RemAx[21], "-", 
            RemAx[22], RemAx[23], RemAx[24], RemAx[25], RemAx[26], 
            RemAx[27], RemAx[28], "-", RemAx[29], RemAx[30], 
            RemAx[31], RemAx[32], RemAx[33], RemAx[34], RemAx[35], 
            "-", RemAx[36], RemAx[37], RemAx[38], RemAx[39], 
            RemAx[40], RemAx[41], RemAx[42])
        MGcomp$RemAxarray <<- tclArray()
        dim(myRemAxarray) <- c(8, 6)
        for (i in 0:7) for (j in 0:5) {
            MGcomp$RemAxarray[[i, j]] <- myRemAxarray[i + 1, 
                j + 1]
        }
        table.RemAx <- tkwidget(frameRemAx, "table", rows = 9, 
            cols = 6, variable = MGcomp$RemAxarray, titlerows = "1", 
            background = "white", width = 6, height = 9)
        tcl(table.RemAx, "width", 0, 12)
        tcl(table.RemAx, "width", 1, 12)
        tcl(table.RemAx, "width", 2, 12)
        tcl(table.RemAx, "width", 3, 12)
        tcl(table.RemAx, "width", 4, 11)
        tcl(table.RemAx, "width", 5, 11)
        clrs = c()
        for (i in 1:ncol(MGvar$activedata)) {
            if (i%%7 == 0) {
                num = 7
            }
            if (i%%7 != 0) {
                num = i%%7
            }
            clrs[i] = brewer.pal(7, "Dark2")[num]
        }
        tabpointfont <- tkfont.create(family = "times", size = 10, 
            weight = "bold")
        count = 0
        for (i in 1:42) {
            if (RemAx[i] != "-") {
                count = count + 1
            }
        }
        if (count > 0) {
            for (i in 1:count) {
                if (i%%7 == 0) {
                  r = 7
                }
                if (i%%7 != 0) {
                  r = i%%7
                }
                if (i < 8) {
                  RaC = paste(r, ",0", sep = "")
                  cind = 0
                }
                if (i > 7 && i < 15) {
                  RaC = paste(r, ",1", sep = "")
                  cind = 1
                }
                if (i > 14 && i < 22) {
                  RaC = paste(r, ",2", sep = "")
                  cind = 2
                }
                if (i > 21 && i < 29) {
                  RaC = paste(r, ",3", sep = "")
                  cind = 3
                }
                if (i > 28 && i < 36) {
                  RaC = paste(r, ",4", sep = "")
                  cind = 4
                }
                if (i > 35 && i < 43) {
                  RaC = paste(r, ",5", sep = "")
                  cind = 5
                }
                for (j in 1:ncol(MGvar$activedata)) {
                  if (RemAx[i] == colnames(MGvar$activedata)[j]) {
                    Col = clrs[j]
                  }
                }
                print(Col)
                tcl(.Tk.ID(table.RemAx), "tag", "celltag", i + 
                  cind, RaC)
                tcl(.Tk.ID(table.RemAx), "tag", "configure", 
                  i + cind, fg = Col, font = tabpointfont)
            }
        }
        tkgrid(table.RemAx)
        MGcomp$table.RemAx <<- table.RemAx
    }
    tableupdate.RemAx <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tabtitle = "Plot1"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tabtitle = "Plot2"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tabtitle = "Plot3"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tabtitle = "Plot4"
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tabtitle = "Plot5"
        }
        RemAx = MGvar$RemovedAxes
        myRemAxarray <- c(tabtitle, RemAx[1], RemAx[2], RemAx[3], 
            RemAx[4], RemAx[5], RemAx[6], RemAx[7], "-", RemAx[8], 
            RemAx[9], RemAx[10], RemAx[11], RemAx[12], RemAx[13], 
            RemAx[14], "-", RemAx[15], RemAx[16], RemAx[17], 
            RemAx[18], RemAx[19], RemAx[20], RemAx[21], "-", 
            RemAx[22], RemAx[23], RemAx[24], RemAx[25], RemAx[26], 
            RemAx[27], RemAx[28], "-", RemAx[29], RemAx[30], 
            RemAx[31], RemAx[32], RemAx[33], RemAx[34], RemAx[35], 
            "-", RemAx[36], RemAx[37], RemAx[38], RemAx[39], 
            RemAx[40], RemAx[41], RemAx[42])
        MGcomp$RemAxarray <<- tclArray()
        dim(myRemAxarray) <- c(8, 6)
        for (i in 0:7) for (j in 0:5) {
            MGcomp$RemAxarray[[i, j]] <- myRemAxarray[i + 1, 
                j + 1]
        }
        tkconfigure(MGcomp$table.RemAx, variable = MGcomp$RemAxarray)
        clrs = c()
        for (i in 1:ncol(MGvar$activedata)) {
            if (i%%7 == 0) {
                num = 7
            }
            if (i%%7 != 0) {
                num = i%%7
            }
            clrs[i] = brewer.pal(7, "Dark2")[num]
        }
        tabpointfont <- tkfont.create(family = "times", size = 10, 
            weight = "bold")
        count = 0
        for (i in 1:42) {
            if (RemAx[i] != "-") {
                count = count + 1
            }
        }
        if (count > 0) {
            for (i in 1:count) {
                if (i%%7 == 0) {
                  r = 7
                }
                if (i%%7 != 0) {
                  r = i%%7
                }
                if (i < 8) {
                  RaC = paste(r, ",0", sep = "")
                  cind = 0
                }
                if (i > 7 && i < 15) {
                  RaC = paste(r, ",1", sep = "")
                  cind = 1
                }
                if (i > 14 && i < 22) {
                  RaC = paste(r, ",2", sep = "")
                  cind = 2
                }
                if (i > 21 && i < 29) {
                  RaC = paste(r, ",3", sep = "")
                  cind = 3
                }
                if (i > 28 && i < 36) {
                  RaC = paste(r, ",4", sep = "")
                  cind = 4
                }
                if (i > 35 && i < 43) {
                  RaC = paste(r, ",5", sep = "")
                  cind = 5
                }
                for (j in 1:ncol(MGvar$activedata)) {
                  if (RemAx[i] == colnames(MGvar$activedata)[j]) {
                    Col = clrs[j]
                  }
                }
                tcl(.Tk.ID(MGcomp$table.RemAx), "tag", "celltag", 
                  i + cind, RaC)
                tcl(.Tk.ID(MGcomp$table.RemAx), "tag", "configure", 
                  i + cind, fg = Col, font = tabpointfont)
            }
        }
    }
    ShowregfromMenu <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            DAx = as.character(tclvalue(MenuDispAx1))
            if (DAx == "1") {
                MGvar$activeplot.showreg.T1 <<- "yes"
                MGvar$activeplot.showreg <<- "yes"
                tkentryconfigure(MainPlotMenu1, 16, state = "active")
            }
            if (DAx == "0") {
                MGvar$activeplot.showreg.T1 <<- "no"
                MGvar$activeplot.showreg <<- "no"
                tkentryconfigure(MainPlotMenu1, 16, state = "disabled")
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            DAx = as.character(tclvalue(MenuDispAx2))
            if (DAx == "1") {
                MGvar$activeplot.showreg.T2 <<- "yes"
                MGvar$activeplot.showreg <<- "yes"
                tkentryconfigure(MainPlotMenu2, 16, state = "active")
            }
            if (DAx == "0") {
                MGvar$activeplot.showreg.T2 <<- "no"
                MGvar$activeplot.showreg <<- "no"
                tkentryconfigure(MainPlotMenu2, 16, state = "disabled")
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            DAx = as.character(tclvalue(MenuDispAx3))
            if (DAx == "1") {
                MGvar$activeplot.showreg.T3 <<- "yes"
                MGvar$activeplot.showreg <<- "yes"
                tkentryconfigure(MainPlotMenu3, 16, state = "active")
            }
            if (DAx == "0") {
                MGvar$activeplot.showreg.T3 <<- "no"
                MGvar$activeplot.showreg <<- "no"
                tkentryconfigure(MainPlotMenu3, 16, state = "disabled")
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            DAx = as.character(tclvalue(MenuDispAx4))
            if (DAx == "1") {
                MGvar$activeplot.showreg.T4 <<- "yes"
                MGvar$activeplot.showreg <<- "yes"
                tkentryconfigure(MainPlotMenu4, 16, state = "active")
            }
            if (DAx == "0") {
                MGvar$activeplot.showreg.T4 <<- "no"
                MGvar$activeplot.showreg <<- "no"
                tkentryconfigure(MainPlotMenu4, 16, state = "disabled")
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            DAx = as.character(tclvalue(MenuDispAx5))
            if (DAx == "1") {
                MGvar$activeplot.showreg.T5 <<- "yes"
                MGvar$activeplot.showreg <<- "yes"
                tkentryconfigure(MainPlotMenu5, 16, state = "active")
            }
            if (DAx == "0") {
                MGvar$activeplot.showreg.T5 <<- "no"
                MGvar$activeplot.showreg <<- "no"
                tkentryconfigure(MainPlotMenu5, 16, state = "disabled")
            }
        }
        tabplot()
    }
    AxesRemove <- function() {
        if (nrow(MGvar$activedata) == 1 && ncol(MGvar$activedata) == 
            1) 
            tkmessageBox(message = "No Data. Please Upload!")
        else {
            RemAx = tktoplevel()
            tkwm.resizable(RemAx, "0", "0")
            tkwm.deiconify(RemAx)
            tkwm.title(RemAx, "Remove Axes")
            tkwm.geometry(RemAx, "320x170")
            RemAxcanvas = tkcanvas(RemAx, bg = col.sec)
            tkplace(RemAxcanvas, relx = 0, rely = 0, relheight = 1, 
                relwidth = 1, `in` = RemAx)
            frameRA <- tkwidget(RemAx, "TitleFrame", text = "Remove Variable Regression Axes", 
                background = "white")
            tkplace(frameRA, relx = 0.01, rely = 0.01, relwidth = 0.98, 
                relheight = 0.98, `in` = RemAx)
            fontsmall <- tkfont.create(family = "times", size = 9)
            tkplace(tklabel(frameRA, text = "Select from the list below the variable you would like to\nremove from the plot. Repeating this process will\npreserve previous axes' removal.", 
                font = fontsmall), relx = 0.06, rely = 0.15, 
                `in` = frameRA)
            DatVarNames <- colnames(MGvar$activedata)
            DatVarNamesTemp <- DatVarNames
            if (length(MGvar$remAxcompindex) > 0) {
                DatVarNamesTemp <- DatVarNames[-MGvar$remAxcompindex]
            }
            DatVar.Var <- tclVar(DatVarNamesTemp[1])
            DatVar.ComboBox <- tkwidget(RemAx, "ComboBox", editable = FALSE, 
                values = DatVarNamesTemp, width = 15, textvariable = DatVar.Var)
            tkplace(DatVar.ComboBox, relx = 0.33, rely = 0.55, 
                `in` = frameRA)
            RemoveAxesBut <- function() {
                count = length(MGvar$remAxcompindex) + 1
                RA <- as.character(tclvalue(DatVar.Var))
                for (i in 1:length(DatVarNames)) {
                  if (RA == DatVarNames[i]) {
                    MGvar$remAxcompindex <<- append(MGvar$remAxcompindex, 
                      i)
                    MGvar$RemovedAxes[count] <<- DatVarNames[i]
                  }
                }
                DatVarNamesTemp = DatVarNames[-MGvar$remAxcompindex]
                tclvalue(DatVar.Var) = (DatVarNamesTemp[1])
                tkconfigure(DatVar.ComboBox, values = DatVarNamesTemp, 
                  textvariable = DatVar.Var)
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                  MGvar$remAxcompindex.T1 <<- MGvar$remAxcompindex
                  MGvar$RemovedAxes.T1 <<- MGvar$RemovedAxes
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                  MGvar$remAxcompindex.T2 <<- MGvar$remAxcompindex
                  MGvar$RemovedAxes.T2 <<- MGvar$RemovedAxes
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                  MGvar$remAxcompindex.T3 <<- MGvar$remAxcompindex
                  MGvar$RemovedAxes.T3 <<- MGvar$RemovedAxes
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                  MGvar$remAxcompindex.T4 <<- MGvar$remAxcompindex
                  MGvar$RemovedAxes.T4 <<- MGvar$RemovedAxes
                }
                if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                  MGvar$remAxcompindex.T5 <<- MGvar$remAxcompindex
                  MGvar$RemovedAxes.T5 <<- MGvar$RemovedAxes
                }
                tableupdate.RemAx()
                tabplot()
            }
            tkplace(tkbutton(RemAx, text = "Remove", width = 15, 
                command = function() RemoveAxesBut()), relx = 0.325, 
                rely = 0.75, `in` = RemAx)
        }
    }
    AxesReplace <- function() {
        activeRow <- as.numeric(tkindex(MGcomp$table.RemAx, "active", 
            "row"))
        activeCol <- as.numeric(tkindex(MGcomp$table.RemAx, "active", 
            "col"))
        tabind <- 1 + 7 * activeCol + (activeRow - 1)
        if (MGvar$RemovedAxes[tabind] == "-") {
            tkmessageBox(message = "The Active Cell is empty! Please select an occupied cell.")
        }
        if (MGvar$RemovedAxes[tabind] != "-") {
            MGvar$RemovedAxes <<- MGvar$RemovedAxes[-tabind]
            MGvar$RemovedAxes <<- append(MGvar$RemovedAxes, "-")
            MGvar$remAxcompindex <<- MGvar$remAxcompindex[-tabind]
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$remAxcompindex.T1 <<- MGvar$remAxcompindex
                MGvar$RemovedAxes.T1 <<- MGvar$RemovedAxes
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$remAxcompindex.T2 <<- MGvar$remAxcompindex
                MGvar$RemovedAxes.T2 <<- MGvar$RemovedAxes
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$remAxcompindex.T3 <<- MGvar$remAxcompindex
                MGvar$RemovedAxes.T3 <<- MGvar$RemovedAxes
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$remAxcompindex.T4 <<- MGvar$remAxcompindex
                MGvar$RemovedAxes.T4 <<- MGvar$RemovedAxes
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$remAxcompindex.T5 <<- MGvar$remAxcompindex
                MGvar$RemovedAxes.T5 <<- MGvar$RemovedAxes
            }
            tabplot()
            tableupdate.RemAx()
        }
    }
    ClearRemindex <- function() {
        MGvar$RemovedPoints <<- as.vector(0)
        for (i in 1:42) {
            MGvar$RemovedPoints[i] <<- "-"
        }
        MGvar$RemovedPoints.T1 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T2 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T3 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T4 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T5 <<- MGvar$RemovedPoints
        MGvar$remindex <<- 0
        MGvar$remindex.T1 <<- 0
        MGvar$remindex.T2 <<- 0
        MGvar$remindex.T3 <<- 0
        MGvar$remindex.T4 <<- 0
        MGvar$remindex.T5 <<- 0
        MGvar$remPcompindex <<- c()
        MGvar$remPcompindex.T1 <<- c()
        MGvar$remPcompindex.T2 <<- c()
        MGvar$remPcompindex.T3 <<- c()
        MGvar$remPcompindex.T4 <<- c()
        MGvar$remPcompindex.T5 <<- c()
        MGvar$RemovedAxes <<- as.vector(0)
        for (i in 1:42) {
            MGvar$RemovedAxes[i] <<- "-"
        }
        MGvar$RemovedAxes.T1 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T2 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T3 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T4 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T5 <<- MGvar$RemovedAxes
        MGvar$remAxindex <<- 0
        MGvar$remAxindex.T1 <<- 0
        MGvar$remAxindex.T2 <<- 0
        MGvar$remAxindex.T3 <<- 0
        MGvar$remAxindex.T4 <<- 0
        MGvar$remAxindex.T5 <<- 0
        MGvar$remAxcompindex <<- c()
        MGvar$remAxcompindex.T1 <<- c()
        MGvar$remAxcompindex.T2 <<- c()
        MGvar$remAxcompindex.T3 <<- c()
        MGvar$remAxcompindex.T4 <<- c()
        MGvar$remAxcompindex.T5 <<- c()
        tableupdate.RemAx()
        tableupdate.Remp()
    }
    tempdistnameassign <- function(name) {
        MGvar$dMeas.T1.temp <<- name
        MGvar$dMeas.T2.temp <<- name
        MGvar$dMeas.T3.temp <<- name
        MGvar$dMeas.T4.temp <<- name
        MGvar$dMeas.T5.temp <<- name
    }
    distnameassign <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$dMeas <<- MGvar$dMeas.T1.temp
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$dMeas <<- MGvar$dMeas.T2.temp
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$dMeas <<- MGvar$dMeas.T3.temp
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$dMeas <<- MGvar$dMeas.T4.temp
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$dMeas <<- MGvar$dMeas.T5.temp
        }
    }
    new.dissim.meas <- function(adata, rbVal = tclvalue(MGvar$dMeasVar)) {
        if (nrow(adata) == 1 && ncol(adata) == 1) {
            tkmessageBox(message = "Please Upload Data!", icon = "error")
        }
        else {
            ActiveWatchCursor()
            if (rbVal == "Euc") {
                MGvar$distmat <<- as.matrix(dist(adata))
                MGvar$dMeas <<- "Euclidean.Distance"
                tempdistnameassign("Euclidean.Distance")
            }
            if (rbVal == "WEuc") {
                WtEuctt <- tktoplevel()
                tkwm.deiconify(WtEuctt)
                tkwm.title(WtEuctt, "Weighted-Euclidean.Distance")
                tkwm.geometry(WtEuctt, "290x140")
                WEuccanvas = tkcanvas(WtEuctt, width = "1128", 
                  height = "756", bg = col.sec)
                tkplace(WEuccanvas, relx = 0, rely = 0, relwidth = 1, 
                  relheight = 1, `in` = WtEuctt)
                frameWE <- tkwidget(WtEuctt, "TitleFrame", text = "Weights Vector", 
                  background = "white")
                tkplace(frameWE, relx = 0.01, rely = 0.01, relwidth = 0.98, 
                  relheight = 0.98, `in` = WtEuctt)
                fontsmall <- tkfont.create(family = "times", 
                  size = 9)
                tkplace(tklabel(frameWE, text = "The Weighted Euclidean Distance requires a   \nvector of weights for computation. This vector \nmust have a length equal to the number of \nvariables in your data.", 
                  font = fontsmall), relx = 0.08, rely = 0.15, 
                  `in` = frameWE)
                On.WEUpload <- function() {
                  fileName <- tclvalue(tkgetOpenFile())
                  if (!nchar(fileName)) 
                    tkmessageBox(message = "No file was selected!")
                  else {
                    MGvar$wts <<- read.table(fileName)
                    if (isTRUE(ncol(MGvar$wts) == nrow(adata)) && 
                      isTRUE(nrow(MGvar$wts) == 1)) {
                      MGvar$distmat <<- Wtd.Euc(adata, MGvar$wts)
                      MGvar$dMeas <<- "Weighted.Euclidean.Distance"
                      tempdistnameassign("Weighted.Euclidean.Distance")
                      MGvar$weights.vec <<- MGvar$wts
                    }
                    if (isTRUE(nrow(MGvar$wts) == nrow(adata)) && 
                      isTRUE(ncol(MGvar$wts) == 1)) {
                      MGvar$distmat <<- Wtd.Euc(adata, t(MGvar$wts))
                      MGvar$dMeas <<- "Weighted Euclidean Distance"
                      tempdistnameassign("Weighted.Euclidean.Distance")
                      MGvar$weights.vec <<- MGvar$wts
                    }
                    else {
                      tkmessageBox(message = "Weights Vector has incorrect dimensions!")
                    }
                  }
                  tkdestroy(WtEuctt)
                }
                tkplace(tkbutton(WtEuctt, text = "Upload", width = 15, 
                  command = function() On.WEUpload()), relx = 0.3, 
                  rely = 0.7, `in` = frameWE)
                tkbind(WtEuctt, "<Return>", On.WEUpload)
            }
            if (rbVal == "Mah") {
                Mahtt <- tktoplevel()
                tkwm.deiconify(Mahtt)
                tkwm.title(Mahtt, "Mahalanobis.Distance")
                tkwm.geometry(Mahtt, "290x200")
                Mahcanvas = tkcanvas(Mahtt, width = "1128", height = "756", 
                  bg = col.sec)
                tkplace(Mahcanvas, relx = 0, rely = 0, relwidth = 1, 
                  relheight = 1, `in` = Mahtt)
                frameMah <- tkwidget(Mahtt, "TitleFrame", text = "Sigma Options", 
                  background = "white")
                tkplace(frameMah, relx = 0.01, rely = 0.01, relwidth = 0.98, 
                  relheight = 0.98, `in` = Mahtt)
                fontsmall <- tkfont.create(family = "times", 
                  size = 9)
                tkplace(tklabel(frameMah, text = "The Mahalanobis Distance requires a sigma   \nmatrix for computation. You may use the\n covariance matrix of your data as a default,\n or upload your own.", 
                  font = fontsmall), relx = 0.1, rely = 0.15, 
                  `in` = frameMah)
                tkplace(tklabel(frameMah, text = "Upload Sigma Matrix", 
                  background = "white"), relx = 0.1, rely = 0.6, 
                  `in` = frameMah)
                MahUp <- tk2checkbutton(Mahtt)
                MahUp.val <- tclVar("0")
                tkconfigure(MahUp, variable = MahUp.val)
                tkplace(MahUp, relx = 0.8, rely = 0.6, `in` = frameMah)
                OnOK.Mah <- function() {
                  MahV = as.character(tclvalue(MahUp.val))
                  if (MahV == "0") {
                    tkdestroy(Mahtt)
                    X <- scale(adata, scale = F)
                    MGvar$distmat <<- Mahalanobis.dist(X, cov(X))
                    MGvar$dMeas <<- "Mahalanobis.Distance"
                    tempdistnameassign("Mahalanobis.Distance")
                  }
                  if (MahV == "1") {
                    tkdestroy(Mahtt)
                    fileName <- tclvalue(tkgetOpenFile())
                    if (!nchar(fileName)) 
                      tkmessageBox(message = "No file was selected!")
                    else {
                      sig <- read.table(fileName)
                      if (isTRUE(nrow(sig) == ncol(adata)) && 
                        isTRUE(ncol(sig) == ncol(adata))) {
                        X <- scale(adata, scale = F)
                        MGvar$distmat <<- Mahalanobis.dist(X, 
                          sig)
                        MGvar$dMeas <<- "Mahalanobis.Distance"
                        tempdistnameassign("Mahalanobis.Distance")
                        MGvar$sigma <<- sig
                      }
                    }
                  }
                }
                tkplace(tkbutton(frameMah, text = "Enter", width = 15, 
                  command = function() OnOK.Mah()), relx = 0.3, 
                  rely = 0.8, `in` = frameMah)
                tkbind(Mahtt, "<Return>", OnOK.Mah)
            }
            if (rbVal == "CBM") {
                MGvar$distmat <<- CityBlock.dist(adata)
                MGvar$dMeas <<- "CityBlock.Metric"
                tempdistnameassign("CityBlock.Metric")
            }
            if (rbVal == "Mink") {
                MGvar$distmat <<- Minkowski.dist(adata)
                MGvar$dMeas <<- "Minkowski.Metric"
                tempdistnameassign("Minkowski.Metric")
            }
            if (rbVal == "Can") {
                MGvar$distmat <<- Canberra.dist(adata)
                MGvar$dMeas <<- "Canberra.Metric"
                tempdistnameassign("Canberra.Metric")
            }
            if (rbVal == "Div") {
                MGvar$distmat <<- Div.dist(adata)
                MGvar$dMeas <<- "Divergence"
                tempdistnameassign("Divergence")
            }
            if (rbVal == "BC") {
                MGvar$distmat <<- BC.dist(adata)
                MGvar$dMeas <<- "Bray-Curtis.Distance"
                tempdistnameassign("Bray-Curtis.Distance")
            }
            if (rbVal == "Soe") {
                MGvar$distmat <<- Soergel.dist(adata)
                MGvar$dMeas <<- "Soergel.Distance"
                tempdistnameassign("Soergel.Metric")
            }
            if (rbVal == "Bhat") {
                MGvar$distmat <<- Bhatt.dist(adata)
                MGvar$dMeas <<- "Bhattacharyya.Distance"
                tempdistnameassign("Bhattacharyya.Distance")
            }
            if (rbVal == "WH") {
                MGvar$distmat <<- WH.dist(adata)
                MGvar$dMeas <<- "Wave-Hedges"
                tempdistnameassign("Wave-Hedges")
            }
            if (rbVal == "Ang") {
                MGvar$distmat <<- AngSep.dist(adata)
                MGvar$dMeas <<- "Angular.Seperation"
                tempdistnameassign("Angular.Seperation")
                MGvar$GenSet.CalcScree <<- "no"
            }
            if (rbVal == "Corr") {
                MGvar$distmat <<- Corr.dist2(adata)
                MGvar$dMeas <<- "Correlation"
                tempdistnameassign("Correlation")
                MGvar$GenSet.CalcScree <<- "no"
            }
            rownames(MGvar$distmat) <<- rownames(adata)
            colnames(MGvar$distmat) <<- rownames(adata)
            MGvar$originaldistmat <<- MGvar$distmat
            ActivedistMat()
            ActivedMeas()
            ActiveArrowCursor()
            tkentryconfigure(EditDataMenu, 1, state = "active")
            distmatcheck()
        }
    }
    distmatcheck <- function() {
        if (any(is.na(MGvar$distmat))) {
            tkmessageBox(message = "Dissimilarity matrix contains at least one NA value. All MDS methods have been disabled. Please Select another dissimilarity metric.")
            tkentryconfigure(MDSMenu, 3, state = "disabled")
            tkentryconfigure(MDSMenu, 4, state = "disabled")
            tkentryconfigure(MDSMenu, 5, state = "disabled")
            tkentryconfigure(MDSMenu, 9, state = "disabled")
            tkentryconfigure(MDSMenu, 10, state = "disabled")
            tkentryconfigure(MDSMenu, 11, state = "disabled")
        }
        else {
            tkentryconfigure(MDSMenu, 3, state = "active")
            tkentryconfigure(MDSMenu, 4, state = "active")
            tkentryconfigure(MDSMenu, 5, state = "active")
            tkentryconfigure(MDSMenu, 9, state = "active")
            tkentryconfigure(MDSMenu, 10, state = "active")
            tkentryconfigure(MDSMenu, 11, state = "active")
            if (any(MGvar$distmat < 0)) {
                tkmessageBox(message = "Dissimilarity matrix contains at least one less than zero element. Some MDS methods have been disabled.")
                tkentryconfigure(MDSMenu, 10, state = "disabled")
                tkentryconfigure(MDSMenu, 11, state = "disabled")
            }
        }
    }
    HomePage.Help <- function() {
        shell.exec("http://mdsgui.r-forge.r-project.org/")
    }
    Documents.Help <- function() {
	  shell.exec("https://www.dropbox.com/sh/mvqn4dwdoi7edac/AAB2uuGELTW3pKtWvT3bpUR5a")	
}


    FuncCode.Display <- function() {
        CodeDisp = tktoplevel()
        tkwm.resizable(CodeDisp, "0", "0")
        tkwm.deiconify(CodeDisp)
        tkwm.title(CodeDisp, "Display Function Code")
        tkwm.geometry(CodeDisp, "320x200")
        CodeDispCanvas = tkcanvas(CodeDisp, bg = col.sec)
        tkplace(CodeDispCanvas, relx = 0, rely = 0, relheight = 1, 
            relwidth = 1, `in` = CodeDisp)
        FrameCD <- tkwidget(CodeDisp, "TitleFrame", text = "Display Function Code", 
            background = "white")
        tkplace(FrameCD, relx = 0.02, relwidth = 0.96, rely = 0.02, 
            relheight = 0.96, `in` = CodeDisp)
        tkplace(tklabel(FrameCD, text = "The following list is a select group of functions used\nby the MDS-GUI. Selecting an item will display the R\ncode for the respective function for the users interest.", 
            background = "white"), relx = 0.05, rely = 0.15, 
            `in` = FrameCD)
        funclist <- c("Classical Scaling", "Metric SMACOF", "2D Plot", 
            "3D Plot", "Plotting Options", "Label Point", "Shepard Plot", 
            "Scree Plot", "Procrustes Calc.", "Load Data", "etc")
        DispFunc.Var <- tclVar(funclist[1])
        DispFunc.ComboBox <- tkwidget(CodeDisp, "ComboBox", editable = FALSE, 
            values = funclist, width = 15, textvariable = DispFunc.Var)
        tkplace(DispFunc.ComboBox, relx = 0.35, rely = 0.5, `in` = CodeDisp)
        DispCode <- function() {
            CodeTxtTL = tktoplevel()
            Codetxtscr <- tkscrollbar(CodeTxtTL, repeatinterval = 5, 
                command = function(...) tkyview(CodeTxt, ...))
            CodetxtXscr <- tkscrollbar(CodeTxtTL, orient = "horizontal", 
                repeatinterval = 5, command = function(...) tkxview(CodeTxt, 
                  ...))
            CodeTxt <- tk2text(CodeTxtTL, yscrollcommand = function(...) tkset(Codetxtscr, 
                ...), xscrollcommand = function(...) tkset(CodetxtXscr, 
                ...), wrap = "none")
            tkgrid(CodeTxt, Codetxtscr)
            tkgrid(CodetxtXscr)
            tkgrid.configure(Codetxtscr, sticky = "ns")
            tkgrid.configure(CodetxtXscr, sticky = "ew")
            copyCText <- function() tcl("event", "generate", 
                .Tk.ID(CodeTxt), "<<Copy>>")
            PasteCText <- function() tcl("event", "generate", 
                .Tk.ID(CodeTxt), "<<Paste>>")
            TexteditPopupMenuC <- tkmenu(CodeTxt, tearoff = FALSE)
            tkadd(TexteditPopupMenuC, "command", label = "Paste <Ctrl-V>", 
                command = PasteCText)
            tkadd(TexteditPopupMenuC, "command", label = "Copy <Ctrl-C>", 
                command = copyCText)
            RightClickNoteTextC <- function(x, y) {
                rootx <- as.integer(tkwinfo("rootx", CodeTxt))
                rooty <- as.integer(tkwinfo("rooty", CodeTxt))
                xTxt <- as.integer(x) + rootx
                yTxt <- as.integer(y) + rooty
                tcl("tk_popup", TexteditPopupMenuC, xTxt, yTxt)
            }
            tkbind(CodeTxt, "<Button-3>", RightClickNoteTextC)
            if (tclvalue(DispFunc.Var) == "Classical Scaling") {
                activecode = capture.output(ClasScal)
                tkwm.title(CodeTxtTL, "Classical Scaling Code")
            }
            if (tclvalue(DispFunc.Var) == "Metric SMACOF") {
                activecode = capture.output(smacof.UCT)
                tkwm.title(CodeTxtTL, "Metric SMACOF Code")
            }
            if (tclvalue(DispFunc.Var) == "2D plot") {
                activecode = capture.output(plotting2D)
                tkwm.title(CodeTxtTL, "2D plot Code")
            }
            if (tclvalue(DispFunc.Var) == "3D plot") {
                activecode = capture.output(plotting3D)
                tkwm.title(CodeTxtTL, "3D plot Code")
            }
            if (tclvalue(DispFunc.Var) == "Plotting Options") {
                activecode = capture.output(ConfPlotOptions)
                tkwm.title(CodeTxtTL, "Plotting Options Code")
            }
            if (tclvalue(DispFunc.Var) == "Label Point") {
                activecode = capture.output(LabelSpecificPoint)
                tkwm.title(CodeTxtTL, "Label Point Code")
            }
            if (tclvalue(DispFunc.Var) == "Shepard Plot") {
                activecode = capture.output(plotShepard)
                tkwm.title(CodeTxtTL, "Shepard Plot Code")
            }
            if (tclvalue(DispFunc.Var) == "Scree Plot") {
                activecode = capture.output(plotScree)
                tkwm.title(CodeTxtTL, "Scree Plot Code")
            }
            if (tclvalue(DispFunc.Var) == "Stress Plot") {
                activecode = capture.output(StressPlot)
                tkwm.title(CodeTxtTL, "Stress Plot Code")
            }
            if (tclvalue(DispFunc.Var) == "Procrustes Calc.") {
                activecode = capture.output(Procrust)
                tkwm.title(CodeTxtTL, "Procrustes Calc. Code")
            }
            if (tclvalue(DispFunc.Var) == "Load Data") {
                activecode = capture.output(LoadDataSettxt)
                tkwm.title(CodeTxtTL, "Load Data Code")
            }
            for (i in 1:length(activecode)) tkinsert(CodeTxt, 
                "end", paste(activecode[i]), sep = "\n")
        }
        ExitBut <- function() {
            tkdestroy(CodeDisp)
        }
        DispCodeTab <- function() {
            if (tclvalue(DispFunc.Var) == "Classical Scaling") {
                activecode = capture.output(ClasScal)
            }
            if (tclvalue(DispFunc.Var) == "Metric SMACOF") {
                activecode = capture.output(smacof.UCT)
            }
            if (tclvalue(DispFunc.Var) == "2D plot") {
                activecode = capture.output(plotting2D)
            }
            if (tclvalue(DispFunc.Var) == "3D plot") {
                activecode = capture.output(plotting3D)
            }
            if (tclvalue(DispFunc.Var) == "Plotting Options") {
                activecode = capture.output(ConfPlotOptions)
            }
            if (tclvalue(DispFunc.Var) == "Label Point") {
                activecode = capture.output(LabelSpecificPoint)
            }
            if (tclvalue(DispFunc.Var) == "Shepard Plot") {
                activecode = capture.output(plotShepard)
            }
            if (tclvalue(DispFunc.Var) == "Scree Plot") {
                activecode = capture.output(plotScree)
            }
            if (tclvalue(DispFunc.Var) == "Stress Plot") {
                activecode = capture.output(StressPlot)
            }
            if (tclvalue(DispFunc.Var) == "Procrustes Calc.") {
                activecode = capture.output(Procrust)
            }
            if (tclvalue(DispFunc.Var) == "Load Data") {
                activecode = capture.output(LoadDataSettxt)
            }
            for (i in 1:length(activecode)) tkinsert(Notetxt, 
                "end", paste(activecode[i]), sep = "\n")
        }
        tkplace(tkbutton(CodeDisp, text = "External Window", 
            width = 15, command = function() DispCode()), relx = 0.12, 
            rely = 0.78, `in` = CodeDisp)
        tkplace(tkbutton(CodeDisp, text = "In Notes/Script Tab", 
            width = 15, command = function() DispCodeTab()), 
            relx = 0.55, rely = 0.78, `in` = CodeDisp)
    }
    ReusePOW <- function() {
        if (as.numeric(tclvalue(MGvar$GS.Win.val)) == 1) {
            if (MGvar$EnActivePlot.switch.T1 == "on") {
                tkdestroy(MGcomp$EActive)
                MGvar$EnActivePlot.switch.T1 <<- "off"
            }
            if (MGvar$EnActivePlot.switch.T2 == "on") {
                tkdestroy(MGcomp$EActive)
                MGvar$EnActivePlot.switch.T2 <<- "off"
            }
            if (MGvar$EnActivePlot.switch.T3 == "on") {
                tkdestroy(MGcomp$EActive)
                MGvar$EnActivePlot.switch.T3 <<- "off"
            }
            if (MGvar$EnActivePlot.switch.T4 == "on") {
                tkdestroy(MGcomp$EActive)
                MGvar$EnActivePlot.switch.T4 <<- "off"
            }
            if (MGvar$EnActivePlot.switch.T5 == "on") {
                tkdestroy(MGcomp$EActive)
                MGvar$EnActivePlot.switch.T5 <<- "off"
            }
            if (MGvar$EStat.switch == "on") {
                tkdestroy(MGcomp$Estat)
                MGvar$EStat.switch <<- "off"
            }
            if (MGvar$EnStress.switch == "on") {
                tkdestroy(MGcomp$EStress)
                MGvar$EnStress.switch <<- "off"
            }
            if (MGvar$EnStress2.switch == "on") {
                tkdestroy(MGcomp$EStress2)
                MGvar$EnStress2.switch <<- "off"
            }
            if (MGvar$EnProcPlot.switch == "on") {
                tkdestroy(MGcomp$EProc)
                MGvar$EnProcPlot.switch <<- "off"
            }
            if (MGvar$EnShep.switch == "on") {
                tkdestroy(MGcomp$EShep)
                MGvar$EnShep.switch <<- "off"
            }
            if (EnScree.switch == "on") {
                tkdestroy(MGcomp$EScree)
                EnScree.switch <<- "off"
            }
            if (popzoomswitch == "on") {
                tkdestroy(MGcomp$zoomplottt)
                popzoomswitch <<- "off"
            }
        }
    }
    exportsweave <- function(All = TRUE, P1 = FALSE, P2 = FALSE, 
        P3 = FALSE, P4 = FALSE, P5 = FALSE) {
        tkmessageBox(message = "The following allows you to export your MDS-GUI output into a pdf file using Sweave and latex. In order to use this, you must have some latex software installed on your computer (e.g. Miktex). In addition to this, the Sweave.sty file must be placed in your R home directory. This file can be copied from your R system folder. e.g. \nProgram Files/R/R-2.13.0/share/textmf/tex/latex/Sweave.sty\nA Tex file will be created in your home folder which may be run to obtain the pdf version.", 
            type = "ok")
        remTab <- tclvalue(MGvar$ActivePlottingTab)
        remStCalc <- MGvar$StCalc
        usrname = tclvalue(MGvar$Active.UserName)
        errorplot <- c()
        if (P1) {
            filename <- paste("MDSGUIPlot1.Rnw")
        }
        if (P2) {
            filename <- paste("MDSGUIPlot2.Rnw")
        }
        if (P3) {
            filename <- paste("MDSGUIPlot3.Rnw")
        }
        if (P4) {
            filename <- paste("MDSGUIPlot4.Rnw")
        }
        if (P5) {
            filename <- paste("MDSGUIPlot5.Rnw")
        }
        if (All) {
            filename <- paste("MDSGUIAllPlots.Rnw")
        }
        cat("\\documentclass[a4paper]{article}\n", file = filename, 
            fill = TRUE)
        cat("\\title{MDS-GUI Output}\n", file = filename, append = TRUE)
        cat("\\author{", usrname, "}\n", file = filename, append = TRUE)
        cat("\\begin{document}\n", file = filename, append = TRUE)
        cat("\\maketitle\n", file = filename, append = TRUE)
        if ((All || P1) && ncol(MGvar$MDSmat.T1 > 1)) {
            if (nrow(MGvar$MDSmat.T1) > 1 && ncol(MGvar$MDSmat.T1) > 
                1) {
                tclvalue(MGvar$ActivePlottingTab) <<- "Tab1"
                ActiveTabChanges()
                cat("\\section{Plot 1}\n", file = filename, append = TRUE)
                cat(paste("Method:                  ", MGvar$MDStype.T1), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                cat(paste("Distance Measure:        ", MGvar$dMeas.T1), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "NormRawStress"
                cat(paste("Normalised Raw Stress:   ", StressCalculation(MGvar$MDSmat.T1, 
                  MGvar$distmat.T1)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-1"
                cat(paste("Stress-1:                ", StressCalculation(MGvar$MDSmat.T1, 
                  MGvar$distmat.T1)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-2"
                cat(paste("Stress-2:                ", StressCalculation(MGvar$MDSmat.T1, 
                  MGvar$distmat.T1)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "CorCoef"
                cat(paste("Correlation Coefficient: ", StressCalculation(MGvar$MDSmat.T1, 
                  MGvar$distmat.T1)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\subsection{MDS Configuration}\n", file = filename, 
                  append = TRUE)
                cat("\\vspace{10mm}\n", file = filename, append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=8,height=7.5>>=\n", 
                  file = filename, append = TRUE)
                cat("plotting2D(MGvar$MDSmat.T1,title=MGvar$activeplot.title.T1,Measure=MGvar$dMeas.T1, showtitle=MGvar$activeplot.title.show.T1, \nshowmeas=MGvar$activeplot.distmeas.T1,xlabel=MGvar$activeplot.xlab.T1, ylabel=MGvar$activeplot.ylab.T1, bgcol=MGvar$activeplot.bg.T1,\npointcex=MGvar$activeplot.cex.T1, showlabs=MGvar$activeplot.labs.T1,showpoints=MGvar$activeplot.showpoints.T1, pointcol=MGvar$activeplot.pointcol.T1,\npointshape=MGvar$activeplot.type.T1, ymeas = MGvar$activeplot.yaxt.T1, xmeas = MGvar$activeplot.xaxt.T1,axcol =  MGvar$activeplot.axescol.T1,\nindexLabeled=MGvar$indexLabeled.T1,zoomedcoords=MGvar$newCoords.T1,showreg = MGvar$activeplot.showreg.T1,regcol=MGvar$activeplot.regcol.T1,\nshowleg=MGvar$activeplot.showleg.T1,Zrat=MGvar$zoominrat.T1,Mvup=MGvar$moveup.T1,Mvdn=MGvar$movedown.T1,Mvlt=MGvar$moveleft.T1,Mvrt=MGvar$moveright.T1,\nPTcolsindex=MGvar$MDSmat.Cols.T1,showdist=MGvar$activeplot.showdist.T1,distcol=MGvar$activeplot.distcol.T1)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Shepard Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotShepard(MGvar$distmat.T1,MGvar$MDSmat.T1,showtitle=MGvar$shepplot.title.show,showlabs=MGvar$shepplot.labs.show,showleg=MGvar$shepplot.leg.show\n,bgcol=MGvar$shepplot.bg,showpoints=MGvar$shepplot.showpoints,pointcex=MGvar$shepplot.cex,pointshape=MGvar$shepplot.type,pointcol=MGvar$shepplot.pointcol\n,showcurve=MGvar$shepplot.curve.show,curvetype=MGvar$shepplot.curve.type,curvecol=MGvar$shepplot.curvecol,xmeas=MGvar$shepplot.Axes.xaxt\n,ymeas=MGvar$shepplot.Axes.yaxt,axcol =  MGvar$shepplot.Axescol,Tab=tclvalue(MGvar$ActivePlottingTab ),showplabs=MGvar$shepplot.showlabels)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Scree Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotScree(MGvar$scree.stress.T1,MGvar$screepoints.current.T1,MGvar$screepoints.best.T1,MGvar$MDS.dimensions,MGvar$Opt.dim.T1,showtitle=MGvar$screeplot.title.show,showlabs=MGvar$screeplot.labs.show,showleg=MGvar$screeplot.leg.show\n,bgcol=MGvar$screeplot.bg,pointtype=MGvar$screeplot.points.show,showCdim=MGvar$screeplot.Cdim.show,showOdim=MGvar$screeplot.Odim.show,CPcol=MGvar$screeplot.Ccol,\nOPcol=MGvar$screeplot.Ocol,showcurve=MGvar$screeplot.curve.show,curvetype=MGvar$screeplot.curve.type,curvecol=MGvar$screeplot.curvecol,Cline=MGvar$screeplot.Cline.show,\nOline=MGvar$screeplot.Oline.show,xmeas=MGvar$screeplot.Axes.xaxt,ymeas=MGvar$screeplot.Axes.yaxt,axcol=MGvar$screeplot.Axescol,Tab=tclvalue(MGvar$ActivePlottingTab ))\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\clearpage\n", file = filename, append = TRUE)
            }
            else {
                errorplot = c(errorplot, "Plot 1")
            }
        }
        if (All || P2) {
            if (nrow(MGvar$MDSmat.T2) > 1 && ncol(MGvar$MDSmat.T2) > 
                1) {
                tclvalue(MGvar$ActivePlottingTab) <<- "Tab2"
                ActiveTabChanges()
                cat("\\section{Plot 2}\n", file = filename, append = TRUE)
                cat(paste("Method:                  ", MGvar$MDStype.T2), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                cat(paste("Distance Measure:        ", MGvar$dMeas.T2), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "NormRawStress"
                cat(paste("Normalised Raw Stress:   ", StressCalculation(MGvar$MDSmat.T2, 
                  MGvar$distmat.T2)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-1"
                cat(paste("Stress-1:                ", StressCalculation(MGvar$MDSmat.T2, 
                  MGvar$distmat.T2)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-2"
                cat(paste("Stress-2:                ", StressCalculation(MGvar$MDSmat.T2, 
                  MGvar$distmat.T2)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "CorCoef"
                cat(paste("Correlation Coefficient: ", StressCalculation(MGvar$MDSmat.T2, 
                  MGvar$distmat.T2)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\subsection{MDS Configuration}\n", file = filename, 
                  append = TRUE)
                cat("\\vspace{35mm}\n", file = filename, append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=8,height=7.5>>=\n", 
                  file = filename, append = TRUE)
                cat("plotting2D(MGvar$MDSmat.T2,title=MGvar$activeplot.title.T2,Measure=MGvar$dMeas.T2, showtitle=MGvar$activeplot.title.show.T2, \nshowmeas=MGvar$activeplot.distmeas.T2,xlabel=MGvar$activeplot.xlab.T2, ylabel=MGvar$activeplot.ylab.T2, bgcol=MGvar$activeplot.bg.T2,\npointcex=MGvar$activeplot.cex.T2, showlabs=MGvar$activeplot.labs.T2,showpoints=MGvar$activeplot.showpoints.T2, pointcol=MGvar$activeplot.pointcol.T2,\npointshape=MGvar$activeplot.type.T2, ymeas = MGvar$activeplot.yaxt.T2, xmeas = MGvar$activeplot.xaxt.T2,axcol =  MGvar$activeplot.axescol.T2,\nindexLabeled=MGvar$indexLabeled.T2,zoomedcoords=MGvar$newCoords.T2,showreg = MGvar$activeplot.showreg.T2,regcol=MGvar$activeplot.regcol.T2,\nshowleg=MGvar$activeplot.showleg.T2,Zrat=MGvar$zoominrat.T2,Mvup=MGvar$moveup.T2,Mvdn=MGvar$movedown.T2,Mvlt=MGvar$moveleft.T2,Mvrt=MGvar$moveright.T2,\nPTcolsindex=MGvar$MDSmat.Cols.T2,showdist=MGvar$activeplot.showdist.T2,distcol=MGvar$activeplot.distcol.T2)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Shepard Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotShepard(MGvar$distmat.T2,MGvar$MDSmat.T2)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Scree Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotScree(MGvar$scree.stress.T2,MGvar$screepoints.current.T2,MGvar$screepoints.best.T2,MGvar$MDS.dimensions,MGvar$Opt.dim.T2,showtitle=MGvar$screeplot.title.show,showlabs=MGvar$screeplot.labs.show,showleg=MGvar$screeplot.leg.show\n,bgcol=MGvar$screeplot.bg,pointtype=MGvar$screeplot.points.show,showCdim=MGvar$screeplot.Cdim.show,showOdim=MGvar$screeplot.Odim.show,CPcol=MGvar$screeplot.Ccol,\nOPcol=MGvar$screeplot.Ocol,showcurve=MGvar$screeplot.curve.show,curvetype=MGvar$screeplot.curve.type,curvecol=MGvar$screeplot.curvecol,Cline=MGvar$screeplot.Cline.show,\nOline=MGvar$screeplot.Oline.show,xmeas=MGvar$screeplot.Axes.xaxt,ymeas=MGvar$screeplot.Axes.yaxt,axcol=MGvar$screeplot.Axescol,Tab=tclvalue(MGvar$ActivePlottingTab ))\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\clearpage\n", file = filename, append = TRUE)
            }
            else {
                errorplot = c(errorplot, "Plot 2")
            }
        }
        if (All || P3) {
            if (nrow(MGvar$MDSmat.T3) > 1 && ncol(MGvar$MDSmat.T3) > 
                1) {
                tclvalue(MGvar$ActivePlottingTab) <<- "Tab3"
                ActiveTabChanges()
                cat("\\section{Plot 3}\n", file = filename, append = TRUE)
                cat(paste("Method:                  ", MGvar$MDStype.T3), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                cat(paste("Distance Measure:        ", MGvar$dMeas.T3), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "NormRawStress"
                cat(paste("Normalised Raw Stress:   ", StressCalculation(MGvar$MDSmat.T3, 
                  MGvar$distmat.T3)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-1"
                cat(paste("Stress-1:                ", StressCalculation(MGvar$MDSmat.T3, 
                  MGvar$distmat.T3)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-2"
                cat(paste("Stress-2:                ", StressCalculation(MGvar$MDSmat.T3, 
                  MGvar$distmat.T3)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "CorCoef"
                cat(paste("Correlation Coefficient: ", StressCalculation(MGvar$MDSmat.T3, 
                  MGvar$distmat.T3)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\subsection{MDS Configuration}\n", file = filename, 
                  append = TRUE)
                cat("\\vspace{35mm}\n", file = filename, append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=8,height=7.5>>=\n", 
                  file = filename, append = TRUE)
                cat("plotting2D(MGvar$MDSmat.T3,title=MGvar$activeplot.title.T3,Measure=MGvar$dMeas.T3, showtitle=MGvar$activeplot.title.show.T3, \nshowmeas=MGvar$activeplot.distmeas.T3,xlabel=MGvar$activeplot.xlab.T3, ylabel=MGvar$activeplot.ylab.T3, bgcol=MGvar$activeplot.bg.T3,\npointcex=MGvar$activeplot.cex.T3, showlabs=MGvar$activeplot.labs.T3,showpoints=MGvar$activeplot.showpoints.T3, pointcol=MGvar$activeplot.pointcol.T3,\npointshape=MGvar$activeplot.type.T3, ymeas = MGvar$activeplot.yaxt.T3, xmeas = MGvar$activeplot.xaxt.T3,axcol =  MGvar$activeplot.axescol.T3,\nindexLabeled=MGvar$indexLabeled.T3,zoomedcoords=MGvar$newCoords.T3,showreg = MGvar$activeplot.showreg.T3,regcol=MGvar$activeplot.regcol.T3,\nshowleg=MGvar$activeplot.showleg.T3,Zrat=MGvar$zoominrat.T3,Mvup=MGvar$moveup.T3,Mvdn=MGvar$movedown.T3,Mvlt=MGvar$moveleft.T3,Mvrt=MGvar$moveright.T3,\nPTcolsindex=MGvar$MDSmat.Cols.T3,showdist=MGvar$activeplot.showdist.T3,distcol=MGvar$activeplot.distcol.T3)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Shepard Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotShepard(MGvar$distmat.T3,MGvar$MDSmat.T3)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Scree Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotScree(MGvar$scree.stress.T3,MGvar$screepoints.current.T3,MGvar$screepoints.best.T3,MGvar$MDS.dimensions,MGvar$Opt.dim.T3)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\clearpage\n", file = filename, append = TRUE)
            }
            else {
                errorplot = c(errorplot, "Plot 3")
            }
        }
        if (All || P4) {
            if (nrow(MGvar$MDSmat.T4) > 1 && ncol(MGvar$MDSmat.T4) > 
                1) {
                tclvalue(MGvar$ActivePlottingTab) <<- "Tab4"
                ActiveTabChanges()
                cat("\\section{Plot 4}\n", file = filename, append = TRUE)
                cat(paste("Method:                  ", MGvar$MDStype.T4), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                cat(paste("Distance Measure:        ", MGvar$dMeas.T4), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "NormRawStress"
                cat(paste("Normalised Raw Stress:   ", StressCalculation(MGvar$MDSmat.T4, 
                  MGvar$distmat.T4)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-1"
                cat(paste("Stress-1:                ", StressCalculation(MGvar$MDSmat.T4, 
                  MGvar$distmat.T4)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-2"
                cat(paste("Stress-2:                ", StressCalculation(MGvar$MDSmat.T4, 
                  MGvar$distmat.T4)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "CorCoef"
                cat(paste("Correlation Coefficient: ", StressCalculation(MGvar$MDSmat.T4, 
                  MGvar$distmat.T4)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\subsection{MDS Configuration}\n", file = filename, 
                  append = TRUE)
                cat("\\vspace{35mm}\n", file = filename, append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=8,height=7.5>>=\n", 
                  file = filename, append = TRUE)
                cat("plotting2D(MGvar$MDSmat.T4,title=MGvar$activeplot.title.T4,Measure=MGvar$dMeas.T4, showtitle=MGvar$activeplot.title.show.T4, \nshowmeas=MGvar$activeplot.distmeas.T4,xlabel=MGvar$activeplot.xlab.T4, ylabel=MGvar$activeplot.ylab.T4, bgcol=MGvar$activeplot.bg.T4,\npointcex=MGvar$activeplot.cex.T4, showlabs=MGvar$activeplot.labs.T4,showpoints=MGvar$activeplot.showpoints.T4, pointcol=MGvar$activeplot.pointcol.T4,\npointshape=MGvar$activeplot.type.T4, ymeas = MGvar$activeplot.yaxt.T4, xmeas = MGvar$activeplot.xaxt.T4,axcol =  MGvar$activeplot.axescol.T4,\nindexLabeled=MGvar$indexLabeled.T4,zoomedcoords=MGvar$newCoords.T4,showreg = MGvar$activeplot.showreg.T4,regcol=MGvar$activeplot.regcol.T4,\nshowleg=MGvar$activeplot.showleg.T4,Zrat=MGvar$zoominrat.T4,Mvup=MGvar$moveup.T4,Mvdn=MGvar$movedown.T4,Mvlt=MGvar$moveleft.T4,Mvrt=MGvar$moveright.T4,\nPTcolsindex=MGvar$MDSmat.Cols.T4,showdist=MGvar$activeplot.showdist.T4,distcol=MGvar$activeplot.distcol.T4)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Shepard Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotShepard(MGvar$distmat.T4,MGvar$MDSmat.T4)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Scree Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotScree(MGvar$scree.stress.T4,MGvar$screepoints.current.T4,MGvar$screepoints.best.T4,MGvar$MDS.dimensions,MGvar$Opt.dim.T4)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\clearpage\n", file = filename, append = TRUE)
            }
            else {
                errorplot = c(errorplot, "Plot 4")
            }
        }
        if (All || P5) {
            if (nrow(MGvar$MDSmat.T5) > 1 && ncol(MGvar$MDSmat.T5) > 
                1) {
                tclvalue(MGvar$ActivePlottingTab) <<- "Tab5"
                ActiveTabChanges()
                cat("\\section{Plot 5}\n", file = filename, append = TRUE)
                cat(paste("Method:                  ", MGvar$MDStype.T5), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                cat(paste("Distance Measure:        ", MGvar$dMeas.T5), 
                  sep = "\n", file = filename, append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "NormRawStress"
                cat(paste("Normalised Raw Stress:   ", StressCalculation(MGvar$MDSmat.T5, 
                  MGvar$distmat.T5)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-1"
                cat(paste("Stress-1:                ", StressCalculation(MGvar$MDSmat.T5, 
                  MGvar$distmat.T5)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "Stress-2"
                cat(paste("Stress-2:                ", StressCalculation(MGvar$MDSmat.T5, 
                  MGvar$distmat.T5)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\\\\n", file = filename, append = TRUE)
                MGvar$StCalc <<- "CorCoef"
                cat(paste("Correlation Coefficient: ", StressCalculation(MGvar$MDSmat.T5, 
                  MGvar$distmat.T5)), sep = "\n", file = filename, 
                  append = TRUE)
                cat("\\subsection{MDS Configuration}\n", file = filename, 
                  append = TRUE)
                cat("\\vspace{35mm}\n", file = filename, append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=8,height=7.5>>=\n", 
                  file = filename, append = TRUE)
                cat("plotting2D(MGvar$MDSmat.T5,title=MGvar$activeplot.title.T5,Measure=MGvar$dMeas.T5, showtitle=MGvar$activeplot.title.show.T5, \nshowmeas=MGvar$activeplot.distmeas.T5,xlabel=MGvar$activeplot.xlab.T5, ylabel=MGvar$activeplot.ylab.T5, bgcol=MGvar$activeplot.bg.T5,\npointcex=MGvar$activeplot.cex.T5, showlabs=MGvar$activeplot.labs.T5,showpoints=MGvar$activeplot.showpoints.T5, pointcol=MGvar$activeplot.pointcol.T5,\npointshape=MGvar$activeplot.type.T5, ymeas = MGvar$activeplot.yaxt.T5, xmeas = MGvar$activeplot.xaxt.T5,axcol =  MGvar$activeplot.axescol.T5,\nindexLabeled=MGvar$indexLabeled.T5,zoomedcoords=MGvar$newCoords.T5,showreg = MGvar$activeplot.showreg.T5,regcol=MGvar$activeplot.regcol.T5,\nshowleg=MGvar$activeplot.showleg.T5,Zrat=MGvar$zoominrat.T5,Mvup=MGvar$moveup.T5,Mvdn=MGvar$movedown.T5,Mvlt=MGvar$moveleft.T5,Mvrt=MGvar$moveright.T5,\nPTcolsindex=MGvar$MDSmat.Cols.T5,showdist=MGvar$activeplot.showdist.T5,distcol=MGvar$activeplot.distcol.T5)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Shepard Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotShepard(MGvar$distmat.T5,MGvar$MDSmat.T5)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\subsection{Scree Plot}\n", file = filename, 
                  append = TRUE)
                cat("\\begin{center}\n", file = filename, append = TRUE)
                cat("<<fig=TRUE,echo=FALSE,width=5,height=4>>=\n", 
                  file = filename, append = TRUE)
                cat("plotScree(MGvar$scree.stress.T5,MGvar$screepoints.current.T5,MGvar$screepoints.best.T5,MGvar$MDS.dimensions,MGvar$Opt.dim.T5)\n", 
                  file = filename, append = TRUE)
                cat("@\n", file = filename, append = TRUE)
                cat("\\end{center}\n", file = filename, append = TRUE)
                cat("\\clearpage\n", file = filename, append = TRUE)
            }
            else {
                errorplot = c(errorplot, "Plot 5")
            }
        }
        cat("\\end{document}\n", file = filename, append = TRUE)
        Sweave(filename)
        tclvalue(MGvar$ActivePlottingTab) <<- remTab
        MGvar$StCalc <<- remStCalc
        ActiveTabChanges()
        if (length(errorplot) > 0) {
            errmessage <- paste("The following plot could not be output:", 
                errorplot[1], if (length(errorplot) > 1) {
                  errorplot[2]
                }, if (length(errorplot) > 2) {
                  errorplot[3]
                }, if (length(errorplot) > 3) {
                  errorplot[4]
                }, if (length(errorplot) > 4) {
                  errorplot[5]
                })
            tkmessageBox(message = errmessage, type = "ok")
        }
    }
    completeplot <- function() {
        curtab <- tclvalue(MGvar$ActivePlottingTab)
        allplot <- function() {
            MGvar$tShepx <<- as.vector(0)
            tabplot()
            if (MGvar$GenSet.CalcShep == "yes") {
                tkrreplot(imgshep)
            }
            tkrreplot(imgscree)
            tableupdate.Remp()
            tableupdate.RemAx()
            tableupdate()
        }
        tclvalue(MGvar$ActivePlottingTab) <<- "Tab1"
        ActiveTabChanges()
        allplot()
        tclvalue(MGvar$ActivePlottingTab) <<- "Tab2"
        ActiveTabChanges()
        allplot()
        tclvalue(MGvar$ActivePlottingTab) <<- "Tab3"
        ActiveTabChanges()
        allplot()
        tclvalue(MGvar$ActivePlottingTab) <<- "Tab4"
        ActiveTabChanges()
        allplot()
        tclvalue(MGvar$ActivePlottingTab) <<- "Tab5"
        ActiveTabChanges()
        tclvalue(MGvar$ActivePlottingTab) <<- curtab
        ActiveTabChanges()
    }
    DeactivateAll <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            tabin = 0
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            tabin = 1
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            tabin = 2
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            tabin = 3
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            tabin = 4
        }
        if (tabin != 0) {
            tcl(myPlottingNB, "tab", 0, state = "disabled")
        }
        if (tabin != 1) {
            tcl(myPlottingNB, "tab", 1, state = "disabled")
        }
        if (tabin != 2) {
            tcl(myPlottingNB, "tab", 2, state = "disabled")
        }
        if (tabin != 3) {
            tcl(myPlottingNB, "tab", 3, state = "disabled")
        }
        if (tabin != 4) {
            tcl(myPlottingNB, "tab", 4, state = "disabled")
        }
        tkentryconfigure(topMenu, 1, state = "disabled")
        tkentryconfigure(topMenu, 2, state = "disabled")
        tkentryconfigure(topMenu, 3, state = "disabled")
        tkentryconfigure(topMenu, 4, state = "disabled")
        tkentryconfigure(topMenu, 5, state = "disabled")
        tkentryconfigure(topMenu, 6, state = "disabled")
    }
    ActivateAll <- function() {
        tcl(myPlottingNB, "tab", 0, state = "normal")
        tcl(myPlottingNB, "tab", 1, state = "normal")
        tcl(myPlottingNB, "tab", 2, state = "normal")
        tcl(myPlottingNB, "tab", 3, state = "normal")
        tcl(myPlottingNB, "tab", 4, state = "normal")
        tkentryconfigure(topMenu, 1, state = "active")
        tkentryconfigure(topMenu, 2, state = "active")
        tkentryconfigure(topMenu, 3, state = "active")
        tkentryconfigure(topMenu, 4, state = "active")
        tkentryconfigure(topMenu, 5, state = "active")
        tkentryconfigure(topMenu, 6, state = "active")
        if (MGvar$MDS.dimensions == 2) {
            ActiveTabChanges()
        }
    }
    MDSGUI.print <- function() {
        try(win.print(), silent = TRUE)
        if (geterrmessage() != "Error in win.print() : unable to start device devWindows\n") {
            printplot()
            dev.off()
        }
    }
    printplot <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            plotting2D(MGvar$MDSmat.T1, title = MGvar$activeplot.title.T1, 
                Measure = MGvar$dMeas.T1, showtitle = MGvar$activeplot.title.show.T1, 
                showmeas = MGvar$activeplot.distmeas.T1, xlabel = MGvar$activeplot.xlab.T1, 
                ylabel = MGvar$activeplot.ylab.T1, bgcol = MGvar$activeplot.bg.T1, 
                pointcex = MGvar$activeplot.cex.T1, showlabs = MGvar$activeplot.labs.T1, 
                showpoints = MGvar$activeplot.showpoints.T1, 
                pointcol = MGvar$activeplot.pointcol.T1, pointshape = MGvar$activeplot.type.T1, 
                ymeas = MGvar$activeplot.yaxt.T1, xmeas = MGvar$activeplot.xaxt.T1, 
                axcol = MGvar$activeplot.axescol.T1, indexLabeled = MGvar$indexLabeled.T1, 
                zoomedcoords = MGvar$newCoords.T1, showreg = MGvar$activeplot.showreg.T1, 
                regcol = MGvar$activeplot.regcol.T1, showleg = MGvar$activeplot.showleg.T1, 
                Zrat = MGvar$zoominrat.T1, Mvup = MGvar$moveup.T1, 
                Mvdn = MGvar$movedown.T1, Mvlt = MGvar$moveleft.T1, 
                Mvrt = MGvar$moveright.T1, PTcolsindex = MGvar$MDSmat.Cols.T1, 
                showdist = MGvar$activeplot.showdist.T1, distcol = MGvar$activeplot.distcol.T1)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            plotting2D(MGvar$MDSmat.T2, title = MGvar$activeplot.title.T2, 
                Measure = MGvar$dMeas.T2, showtitle = MGvar$activeplot.title.show.T2, 
                showmeas = MGvar$activeplot.distmeas.T2, xlabel = MGvar$activeplot.xlab.T2, 
                ylabel = MGvar$activeplot.ylab.T2, bgcol = MGvar$activeplot.bg.T2, 
                pointcex = MGvar$activeplot.cex.T2, showlabs = MGvar$activeplot.labs.T2, 
                showpoints = MGvar$activeplot.showpoints.T2, 
                pointcol = MGvar$activeplot.pointcol.T2, pointshape = MGvar$activeplot.type.T2, 
                ymeas = MGvar$activeplot.yaxt.T2, xmeas = MGvar$activeplot.xaxt.T2, 
                axcol = MGvar$activeplot.axescol.T2, indexLabeled = MGvar$indexLabeled.T2, 
                zoomedcoords = MGvar$newCoords.T2, showreg = MGvar$activeplot.showreg.T2, 
                regcol = MGvar$activeplot.regcol.T2, showleg = MGvar$activeplot.showleg.T2, 
                Zrat = MGvar$zoominrat.T2, Mvup = MGvar$moveup.T2, 
                Mvdn = MGvar$movedown.T2, Mvlt = MGvar$moveleft.T2, 
                Mvrt = MGvar$moveright.T2, PTcolsindex = MGvar$MDSmat.Cols.T2, 
                showdist = MGvar$activeplot.showdist.T2, distcol = MGvar$activeplot.distcol.T2)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            plotting2D(MGvar$MDSmat.T3, title = MGvar$activeplot.title.T3, 
                Measure = MGvar$dMeas.T3, showtitle = MGvar$activeplot.title.show.T3, 
                showmeas = MGvar$activeplot.distmeas.T3, xlabel = MGvar$activeplot.xlab.T3, 
                ylabel = MGvar$activeplot.ylab.T3, bgcol = MGvar$activeplot.bg.T3, 
                pointcex = MGvar$activeplot.cex.T3, showlabs = MGvar$activeplot.labs.T3, 
                showpoints = MGvar$activeplot.showpoints.T3, 
                pointcol = MGvar$activeplot.pointcol.T3, pointshape = MGvar$activeplot.type.T3, 
                ymeas = MGvar$activeplot.yaxt.T3, xmeas = MGvar$activeplot.xaxt.T3, 
                axcol = MGvar$activeplot.axescol.T3, indexLabeled = MGvar$indexLabeled.T3, 
                zoomedcoords = MGvar$newCoords.T3, showreg = MGvar$activeplot.showreg.T3, 
                regcol = MGvar$activeplot.regcol.T3, showleg = MGvar$activeplot.showleg.T3, 
                Zrat = MGvar$zoominrat.T3, Mvup = MGvar$moveup.T3, 
                Mvdn = MGvar$movedown.T3, Mvlt = MGvar$moveleft.T3, 
                Mvrt = MGvar$moveright.T3, PTcolsindex = MGvar$MDSmat.Cols.T3, 
                showdist = MGvar$activeplot.showdist.T3, distcol = MGvar$activeplot.distcol.T3)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            plotting2D(MGvar$MDSmat.T4, title = MGvar$activeplot.title.T4, 
                Measure = MGvar$dMeas.T4, showtitle = MGvar$activeplot.title.show.T4, 
                showmeas = MGvar$activeplot.distmeas.T4, xlabel = MGvar$activeplot.xlab.T4, 
                ylabel = MGvar$activeplot.ylab.T4, bgcol = MGvar$activeplot.bg.T4, 
                pointcex = MGvar$activeplot.cex.T4, showlabs = MGvar$activeplot.labs.T4, 
                showpoints = MGvar$activeplot.showpoints.T4, 
                pointcol = MGvar$activeplot.pointcol.T4, pointshape = MGvar$activeplot.type.T4, 
                ymeas = MGvar$activeplot.yaxt.T4, xmeas = MGvar$activeplot.xaxt.T4, 
                axcol = MGvar$activeplot.axescol.T4, indexLabeled = MGvar$indexLabeled.T4, 
                zoomedcoords = MGvar$newCoords.T4, showreg = MGvar$activeplot.showreg.T4, 
                regcol = MGvar$activeplot.regcol.T4, showleg = MGvar$activeplot.showleg.T4, 
                Zrat = MGvar$zoominrat.T4, Mvup = MGvar$moveup.T4, 
                Mvdn = MGvar$movedown.T4, Mvlt = MGvar$moveleft.T4, 
                Mvrt = MGvar$moveright.T4, PTcolsindex = MGvar$MDSmat.Cols.T4, 
                showdist = MGvar$activeplot.showdist.T4, distcol = MGvar$activeplot.distcol.T4)
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            plotting2D(MGvar$MDSmat.T5, title = MGvar$activeplot.title.T5, 
                Measure = MGvar$dMeas.T5, showtitle = MGvar$activeplot.title.show.T5, 
                showmeas = MGvar$activeplot.distmeas.T5, xlabel = MGvar$activeplot.xlab.T5, 
                ylabel = MGvar$activeplot.ylab.T5, bgcol = MGvar$activeplot.bg.T5, 
                pointcex = MGvar$activeplot.cex.T5, showlabs = MGvar$activeplot.labs.T5, 
                showpoints = MGvar$activeplot.showpoints.T5, 
                pointcol = MGvar$activeplot.pointcol.T5, pointshape = MGvar$activeplot.type.T5, 
                ymeas = MGvar$activeplot.yaxt.T5, xmeas = MGvar$activeplot.xaxt.T5, 
                axcol = MGvar$activeplot.axescol.T5, indexLabeled = MGvar$indexLabeled.T5, 
                zoomedcoords = MGvar$newCoords.T5, showreg = MGvar$activeplot.showreg.T5, 
                regcol = MGvar$activeplot.regcol.T5, showleg = MGvar$activeplot.showleg.T5, 
                Zrat = MGvar$zoominrat.T5, Mvup = MGvar$moveup.T5, 
                Mvdn = MGvar$movedown.T5, Mvlt = MGvar$moveleft.T5, 
                Mvrt = MGvar$moveright.T5, PTcolsindex = MGvar$MDSmat.Cols.T5, 
                showdist = MGvar$activeplot.showdist.T5, distcol = MGvar$activeplot.distcol.T5)
        }
    }
    DispPopUpHelp <- function() {
        if (as.character(tclvalue(MGvar$DispHelp.Var)) == 1) {
            tk2tip(img, "Click any point on the plot to identify\nthe closest point and label it")
            tk2tip(img2, "Click any point on the plot to identify\nthe closest point and label it")
            tk2tip(img3, "Click any point on the plot to identify\nthe closest point and label it")
            tk2tip(img4, "Click any point on the plot to identify\nthe closest point and label it")
            tk2tip(img5, "Click any point on the plot to identify\nthe closest point and label it")
            tk2tip(procimg, "Result of Procrustes Analysis Shown here.\nRequires two plots in any Plot1-5 with the\nsame object space. Go to Multivariate Tools\n-> Procrustes Analysis.")
            tk2tip(img3Dstat, "Provides 3D result when p=3, and `Static Plot'\nis selected. Uses scatterplot3d package. To\nchange p, go to Multivariate Tools -> MDS Options.")
            tk2tip(imgshep, "Click any point to label it. The corresponding\nline is drawn on the active configuration plot")
            tk2tip(imgstress, "Plots Stress over number of iterations.\nOnly applies to Metric and Non-Metric\nSMACOF.")
            tk2tip(imgstress2, "Plots the logged differences of stress\nbetween iterations. Only applies to\nMetric and Non-Metric SMACOF.")
            tk2tip(imgscree, "Plots Stress over 10 dimensions. If data is very\nlarge, this process may take very long to compute.\nIf so, Click `End Process', go to General Settings\nand deactivate the Scree Plot.")
            tk2tip(imgseczoom, "The Secondary Zoomed Plot is not interactive.\n Use Active Tab to add labels and re-zoom")
            tk2tip(EndButton, "Will end the current MDS process. Use\nis only suggested when process is taking\nmuch longer than expected. May cause\nunexpected runtime errors. Use with caution.")
            tk2tip(MDSGUI.ProgressBar, "Tracks the number of iterations as a percentage\nof the maximum iterations. Bar will progress both\nduring the configuration calculation and\nthroughout the Scree Plot Calculation. Maximum\niterations can be changed via General -> General Settings.")
            tk2tip(frameBottomStrip, "Provides information regarding, name\nof current data set, the active plotting area,\nand number of plotting dimensions.")
            tk2tip(GoToButton, "When either the `Static 3D Plot' or\n`Procrustes' tab are focused, this\nbutton will jump to the active main\narea (Plot1-5)")
            tk2tip(tb3.1, "Provides information of all Plotting areas\n(Plot1-5) for direct comparison.")
            tk2tip(tb3.2, "Provides information of points that have\nbeen removed from each Plotting Area\n(Plot1-5). Table is specific to the active area.\nPoints are removed via the `Remove a Point'\noption of the via accessed via the right click\nmenu of the plotting area.")
            tk2tip(tb3.3, "Provides information of axes that have\nbeen removed from each Plotting Area\n(Plot1-5). Table is specific to the active area.\nVariable Axes are removed via the `Remove Axes\nof Variable(s)' option of the via accessed via the\nright click menu of the plotting area.")
        }
        if (as.character(tclvalue(MGvar$DispHelp.Var)) == 0) {
            tk2tip(img, "")
            tk2tip(img2, "")
            tk2tip(img3, "")
            tk2tip(img4, "")
            tk2tip(img5, "")
            tk2tip(procimg, "")
            tk2tip(img3Dstat, "")
            tk2tip(imgshep, "")
            tk2tip(imgstress, "")
            tk2tip(imgstress2, "")
            tk2tip(imgscree, "")
            tk2tip(imgseczoom, "")
            tk2tip(EndButton, "")
            tk2tip(MDSGUI.ProgressBar, "")
            tk2tip(frameBottomStrip, "")
            tk2tip(GoToButton, "")
            tk2tip(tb3.1, "")
            tk2tip(tb3.2, "")
            tk2tip(tb3.3, "")
        }
    }
    ChangeTabSC <- function(num) {
        if (num == 1) {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab1"
            tk2notetab.select(myPlottingNB, "Plot1")
        }
        if (num == 2) {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab2"
            tk2notetab.select(myPlottingNB, "Plot2")
        }
        if (num == 3) {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab3"
            tk2notetab.select(myPlottingNB, "Plot3")
        }
        if (num == 4) {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab4"
            tk2notetab.select(myPlottingNB, "Plot4")
        }
        if (num == 5) {
            tclvalue(MGvar$ActivePlottingTab) <<- "Tab5"
            tk2notetab.select(myPlottingNB, "Plot5")
        }
        ActiveTabChanges()
    }
    DefaultPtCols <- function() {
        if (length(MGvar$ClasTab) == 0) {
            PointColInitialise()
        }
        if (length(MGvar$ClasTab) > 0) {
            potcols <- c(brewer.pal(9, "Set1"), brewer.pal(12, 
                "Paired"))
            for (i in 1:nrow(MGvar$activedata)) {
                for (j in 1:length(MGvar$ClasTab)) {
                  if (MGvar$ClasVec[i] == names(MGvar$ClasTab)[j]) {
                    if (j == length(potcols)) {
                      num = length(potcols)
                    }
                    if (j != length(potcols)) {
                      num = j%%length(potcols)
                    }
                    MGvar$MDSmat.Cols[i] <<- potcols[num]
                  }
                }
            }
            for (i in 1:length(MGvar$ClasTab)) {
                MGvar$ClasTabCols <<- c(MGvar$ClasTabCols, potcols[i])
                names(MGvar$ClasTabCols)[i] <<- names(MGvar$ClasTab)[i]
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
                MGvar$MDSmat.Cols.T1 <<- MGvar$MDSmat.Cols
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
                MGvar$MDSmat.Cols.T2 <<- MGvar$MDSmat.Cols
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
                MGvar$MDSmat.Cols.T3 <<- MGvar$MDSmat.Cols
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
                MGvar$MDSmat.Cols.T4 <<- MGvar$MDSmat.Cols
            }
            if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
                MGvar$MDSmat.Cols.T5 <<- MGvar$MDSmat.Cols
            }
        }
        tabplot()
    }
    OptDim.Scree <- function(scree.stress) {
        ndim = MGvar$MDS.dimensions
        bestchange = 0
        cap = min(MGvar$maxdims, 12)
        for (i in 2:(cap - 2)) {
            m1 = (scree.stress[i - 1] - scree.stress[i])
            thet1 = atan(m1/1) * 57.2957795
            m2 = (MGvar$scree.stress[i] - scree.stress[i + 1])
            thet2 = atan(1/m2) * 57.2957795
            newthet = 360 - (90 + thet1 + thet2)
            if (newthet < bestchange) {
                bestchange = newthet
                MGvar$Opt.dim <<- i
            }
            print(newthet)
        }
        return(MGvar$Opt.dim)
    }
    MGvar <- list()
    Initialise.Variables <- function() {
        MGvar$activedata <<- as.matrix(0)
        MGvar$distmat <<- as.matrix(0)
        MGvar$distmat.T1 <<- as.matrix(0)
        MGvar$distmat.T2 <<- as.matrix(0)
        MGvar$distmat.T3 <<- as.matrix(0)
        MGvar$distmat.T4 <<- as.matrix(0)
        MGvar$distmat.T5 <<- as.matrix(0)
        MGvar$simmat <<- as.matrix(0)
        MGvar$simmat.T1 <<- as.matrix(0)
        MGvar$simmat.T2 <<- as.matrix(0)
        MGvar$simmat.T3 <<- as.matrix(0)
        MGvar$simmat.T4 <<- as.matrix(0)
        MGvar$simmat.T5 <<- as.matrix(0)
        MGvar$cormat <<- as.matrix(0)
        MGvar$cormat.T1 <<- as.matrix(0)
        MGvar$cormat.T2 <<- as.matrix(0)
        MGvar$cormat.T3 <<- as.matrix(0)
        MGvar$cormat.T4 <<- as.matrix(0)
        MGvar$cormat.T5 <<- as.matrix(0)
        MGvar$MDSmat <<- as.matrix(0)
        MGvar$MDSmat.T1 <<- as.matrix(0)
        MGvar$MDSmat.T2 <<- as.matrix(0)
        MGvar$MDSmat.T3 <<- as.matrix(0)
        MGvar$MDSmat.T4 <<- as.matrix(0)
        MGvar$MDSmat.T5 <<- as.matrix(0)
        MGvar$CSinitConf <<- "yes"
        MGvar$MDStype <<- ""
        MGvar$MDStype.T1 <<- "-"
        MGvar$MDStype.T2 <<- "-"
        MGvar$MDStype.T3 <<- "-"
        MGvar$MDStype.T4 <<- "-"
        MGvar$MDStype.T5 <<- "-"
        MGvar$MDStype.3S <<- "-"
        MGvar$MDStype.3R <<- "-"
        MGvar$weights.vec <<- as.matrix(0)
        MGvar$sigma <<- as.matrix(0)
        MGvar$threeDconf <<- as.matrix(0)
        MGvar$temp.dims <<- 2
        MGvar$datatitle <<- ""
        MGvar$activeplot.title <<- ""
        MGvar$activeplot.title.T1 <<- ""
        MGvar$activeplot.title.T2 <<- ""
        MGvar$activeplot.title.T3 <<- ""
        MGvar$activeplot.title.T4 <<- ""
        MGvar$activeplot.title.T5 <<- ""
        MGvar$activeshepplot.title <<- " "
        MGvar$activescreeplot.title <<- " "
        MGvar$dMeasVar <<- tclVar("Euc")
        MGvar$dMeas <<- " "
        MGvar$dMeas.T1 <<- "-"
        MGvar$dMeas.T2 <<- "-"
        MGvar$dMeas.T3 <<- "-"
        MGvar$dMeas.T4 <<- "-"
        MGvar$dMeas.T5 <<- "-"
        MGvar$dMeas.3S <<- "-"
        MGvar$dMeas.3R <<- "-"
        MGvar$ActivePlottingTab <<- tclVar("Tab1")
        MGvar$rb.PT.Value <<- tclVar("PTab1")
        MGvar$UserName <<- tclVar("")
        MGvar$Active.UserName <<- tclVar("None")
        MGvar$distdat <<- NULL
        MGvar$MDS.dimensions <<- 2
        MGvar$MDS.dimensions.T1 <<- "-"
        MGvar$MDS.dimensions.T2 <<- "-"
        MGvar$MDS.dimensions.T3 <<- "-"
        MGvar$MDS.dimensions.T4 <<- "-"
        MGvar$MDS.dimensions.T5 <<- "-"
        MGvar$MDS.dimensions.3S <<- "-"
        MGvar$MDS.dimensions.3R <<- "-"
        MGvar$maxdims <<- 10
        MGvar$tempdims <<- 2
        MGvar$Opt.dim <<- 1
        MGvar$Opt.dim.T1 <<- 1
        MGvar$Opt.dim.T2 <<- 1
        MGvar$Opt.dim.T3 <<- 1
        MGvar$Opt.dim.T4 <<- 1
        MGvar$Opt.dim.T5 <<- 1
        MGvar$scree.stress <<- as.vector(0)
        MGvar$scree.stress.T1 <<- as.vector(0)
        MGvar$scree.stress.T2 <<- as.vector(0)
        MGvar$scree.stress.T3 <<- as.vector(0)
        MGvar$scree.stress.T4 <<- as.vector(0)
        MGvar$scree.stress.T5 <<- as.vector(0)
        MGvar$screepoints.current <<- as.vector(0)
        MGvar$screepoints.current.T1 <<- as.vector(0)
        MGvar$screepoints.current.T2 <<- as.vector(0)
        MGvar$screepoints.current.T3 <<- as.vector(0)
        MGvar$screepoints.current.T4 <<- as.vector(0)
        MGvar$screepoints.current.T5 <<- as.vector(0)
        MGvar$screepoints.best <<- as.vector(0)
        MGvar$screepoints.best.T1 <<- as.vector(0)
        MGvar$screepoints.best.T2 <<- as.vector(0)
        MGvar$screepoints.best.T3 <<- as.vector(0)
        MGvar$screepoints.best.T4 <<- as.vector(0)
        MGvar$screepoints.best.T5 <<- as.vector(0)
        MGvar$MDS.iter.max <<- 1000
        MGvar$MDS.iter.T1 <<- "-"
        MGvar$MDS.iter.T2 <<- "-"
        MGvar$MDS.iter.T3 <<- "-"
        MGvar$MDS.iter.T4 <<- "-"
        MGvar$MDS.iter.T5 <<- "-"
        MGvar$MDS.iter.3S <<- "-"
        MGvar$MDS.iter.3R <<- "-"
        MGvar$MDS.tol <<- 1e-05
        MGvar$MDS.tol.T1 <<- "-"
        MGvar$MDS.tol.T2 <<- "-"
        MGvar$MDS.tol.T3 <<- "-"
        MGvar$MDS.tol.T4 <<- "-"
        MGvar$MDS.tol.T5 <<- "-"
        MGvar$MDS.tol.3S <<- "-"
        MGvar$MDS.tol.3R <<- "-"
        MGvar$zoomedplot.title.show <<- "yes"
        MGvar$zoomedplot.distmeas <<- "no"
        MGvar$zoomedplot.xlab <<- ""
        MGvar$zoomedplot.ylab <<- ""
        MGvar$zoomedplot.bg <<- "white"
        MGvar$zoomedplot.bg.temp <<- "white"
        MGvar$zoomedplot.cex <<- 1
        MGvar$zoomedplot.labs <<- "yes"
        MGvar$zoomedplot.showpoints <<- "no"
        MGvar$zoomedplot.pointcol <<- "black"
        MGvar$zoomedplot.pointcol.temp <<- "black"
        MGvar$zoomedplot.type <<- 1
        MGvar$zoomedplot.xaxt <<- "n"
        MGvar$zoomedplot.yaxt <<- "n"
        MGvar$zoomedplot.axescol <<- "black"
        MGvar$zoomedplot.axescol.temp <<- "black"
        MGvar$zoomedplot.showzoom <<- "no"
        MGvar$zoomedplot.showleg <<- "no"
        MGvar$Zoom.Main.val <<- tclVar("1")
        MGvar$Zoom.Dist.val <<- tclVar("1")
        MGvar$Zoom.Ylab.val <<- tclVar("0")
        MGvar$Zoom.Xlab.val <<- tclVar("0")
        MGvar$Zoom.Points.val <<- tclVar("0")
        MGvar$Zoom.Labels.val <<- tclVar("1")
        MGvar$Zoom.PT.var <<- tclVar("Empty Circles")
        MGvar$Zoom.ShowZoom.val <<- tclVar("0")
        MGvar$Zoom.AxesMeas.val <<- tclVar("0")
        MGvar$MDSStress <<- ""
        MGvar$MDSStress.T1 <<- "-"
        MGvar$MDSStress.T2 <<- "-"
        MGvar$MDSStress.T3 <<- "-"
        MGvar$MDSStress.T4 <<- "-"
        MGvar$MDSStress.T5 <<- "-"
        MGvar$MDSStress.3S <<- "-"
        MGvar$MDSStress.3R <<- "-"
        MGvar$GenSet.CalcMDS <<- "yes"
        MGvar$GenSet.CalcShep <<- "yes"
        MGvar$GenSet.CalcStress <<- "yes"
        MGvar$GenSet.CalcScree <<- "yes"
        MGvar$GenSet.RUWin <<- "no"
        MGvar$GenSet.DefLabs <<- "yes"
        MGvar$GenSet.DefPts <<- "no"
        MGvar$GenSet.DefReg <<- "no"
        MGvar$GenSet.KeepShep <<- "no"
        MGvar$GenSet.ClearCols <<- "no"
        MGvar$GenSet.ClearIL <<- "yes"
        MGvar$GenSet.ILPos <<- 3
        MGvar$GenSet.UpConf <<- "yes"
        MGvar$GenSet.UpShep <<- "yes"
        MGvar$GenSet.UpStress <<- "yes"
        MGvar$GenSet.UpProg <<- "yes"
        MGvar$MDSops.startconfig <<- tclVar("ClasScal")
        MGvar$MDSops.UseEC.Plot <<- tclVar("Plot1")
        MGvar$startconfig.UseEC.plot <<- "Plot1"
        MGvar$MDSops.StressCalc <<- tclVar("Norm. Raw Stress")
        MGvar$StCalc <<- "NormRawStress"
        MGvar$PlottingDimX <<- 1
        MGvar$PlottingDimY <<- 2
        MGvar$PlottingDimX.3D <<- 1
        MGvar$PlottingDimY.3D <<- 2
        MGvar$PlottingDimZ.3D <<- 3
        MGvar$TabDims.T1 <<- "-"
        MGvar$TabDims.T2 <<- "-"
        MGvar$TabDims.T3 <<- "-"
        MGvar$TabDims.T4 <<- "-"
        MGvar$TabDims.T5 <<- "-"
        MGvar$TabDims.3S <<- "-"
        MGvar$TabDims.3R <<- "-"
        MGvar$MDSmat.LD <<- as.matrix(0)
        MGvar$fromscree <<- FALSE
        MGvar$breakscree <<- FALSE
        MGvar$stressitervec <<- as.vector(0)
        MGvar$endprocess <<- "no"
        MGvar$proctime <<- 0
        MGvar$removedpoints <<- FALSE
        MGvar$removedpoints.T1 <<- FALSE
        MGvar$removedpoints.T2 <<- FALSE
        MGvar$removedpoints.T3 <<- FALSE
        MGvar$removedpoints.T4 <<- FALSE
        MGvar$removedpoints.T5 <<- FALSE
        MGvar$removedpointsactivedata <<- as.matrix(0)
        MGvar$removedpointsactivedata.T1 <<- as.matrix(0)
        MGvar$removedpointsactivedata.T2 <<- as.matrix(0)
        MGvar$removedpointsactivedata.T3 <<- as.matrix(0)
        MGvar$removedpointsactivedata.T4 <<- as.matrix(0)
        MGvar$removedpointsactivedata.T5 <<- as.matrix(0)
        MGvar$ClasVec <<- c()
        MGvar$ClasTabCols <<- c()
        MGvar$ClasTab <<- c()
        MGvar$RemovedPoints <<- as.vector(0)
        for (i in 1:42) {
            MGvar$RemovedPoints[i] <<- "-"
        }
        MGvar$RemovedPoints.T1 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T2 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T3 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T4 <<- MGvar$RemovedPoints
        MGvar$RemovedPoints.T5 <<- MGvar$RemovedPoints
        MGvar$remindex <<- 0
        MGvar$remindex.T1 <<- 0
        MGvar$remindex.T2 <<- 0
        MGvar$remindex.T3 <<- 0
        MGvar$remindex.T4 <<- 0
        MGvar$remindex.T5 <<- 0
        MGvar$remPcompindex <<- c()
        MGvar$remPcompindex.T1 <<- c()
        MGvar$remPcompindex.T2 <<- c()
        MGvar$remPcompindex.T3 <<- c()
        MGvar$remPcompindex.T4 <<- c()
        MGvar$remPcompindex.T5 <<- c()
        MGvar$RemovedAxes <<- as.vector(0)
        for (i in 1:42) {
            MGvar$RemovedAxes[i] <<- "-"
        }
        MGvar$RemovedAxes.T1 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T2 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T3 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T4 <<- MGvar$RemovedAxes
        MGvar$RemovedAxes.T5 <<- MGvar$RemovedAxes
        MGvar$remAxindex <<- 0
        MGvar$remAxindex.T1 <<- 0
        MGvar$remAxindex.T2 <<- 0
        MGvar$remAxindex.T3 <<- 0
        MGvar$remAxindex.T4 <<- 0
        MGvar$remAxindex.T5 <<- 0
        MGvar$remAxcompindex <<- c()
        MGvar$remAxcompindex.T1 <<- c()
        MGvar$remAxcompindex.T2 <<- c()
        MGvar$remAxcompindex.T3 <<- c()
        MGvar$remAxcompindex.T4 <<- c()
        MGvar$remAxcompindex.T5 <<- c()
        MGvar$DistFunc <<- as.vector(0)
        MGvar$DistFunc.T1 <<- as.vector(0)
        MGvar$DistFunc.T2 <<- as.vector(0)
        MGvar$DistFunc.T3 <<- as.vector(0)
        MGvar$DistFunc.T4 <<- as.vector(0)
        MGvar$DistFunc.T5 <<- as.vector(0)
        MGvar$is.Metric <<- ""
        MGvar$is.Metric.T1 <<- ""
        MGvar$is.Metric.T2 <<- ""
        MGvar$is.Metric.T3 <<- ""
        MGvar$is.Metric.T4 <<- ""
        MGvar$is.Metric.T5 <<- ""
        MGvar$ESt1hscale <<- 1.8
        MGvar$ESt1vscale <<- 1.3
        MGvar$ESt2hscale <<- 1.8
        MGvar$ESt2vscale <<- 1.3
        MGvar$main.initpos <<- TRUE
        MGvar$first.xCoord <<- 0
        MGvar$first.yCoord <<- 0
        MGvar$latest.xCoord <<- 0
        MGvar$latest.yCoord <<- 0
        MGvar$ActivePoint <<- 0
        MGvar$movingpoint <<- TRUE
        MGvar$DispHelp.Var <<- tclVar(FALSE)
        MGvar$prevMat <<- as.matrix(0)
        MGvar$activeplot.title.show <<- "yes"
        MGvar$activeplot.title.show.T1 <<- "yes"
        MGvar$activeplot.title.show.T2 <<- "yes"
        MGvar$activeplot.title.show.T3 <<- "yes"
        MGvar$activeplot.title.show.T4 <<- "yes"
        MGvar$activeplot.title.show.T5 <<- "yes"
        MGvar$activeplot.distmeas <<- "yes"
        MGvar$activeplot.distmeas.T1 <<- "yes"
        MGvar$activeplot.distmeas.T2 <<- "yes"
        MGvar$activeplot.distmeas.T3 <<- "yes"
        MGvar$activeplot.distmeas.T4 <<- "yes"
        MGvar$activeplot.distmeas.T5 <<- "yes"
        MGvar$activeplot.xlab <<- ""
        MGvar$activeplot.xlab.T1 <<- ""
        MGvar$activeplot.xlab.T2 <<- ""
        MGvar$activeplot.xlab.T3 <<- ""
        MGvar$activeplot.xlab.T4 <<- ""
        MGvar$activeplot.xlab.T5 <<- ""
        MGvar$activeplot.ylab <<- ""
        MGvar$activeplot.ylab.T1 <<- ""
        MGvar$activeplot.ylab.T2 <<- ""
        MGvar$activeplot.ylab.T3 <<- ""
        MGvar$activeplot.ylab.T4 <<- ""
        MGvar$activeplot.ylab.T5 <<- ""
        MGvar$activeplot.bg <<- "white"
        MGvar$activeplot.bg.temp <<- "white"
        MGvar$activeplot.bg.T1 <<- "white"
        MGvar$activeplot.bg.T2 <<- "white"
        MGvar$activeplot.bg.T3 <<- "white"
        MGvar$activeplot.bg.T4 <<- "white"
        MGvar$activeplot.bg.T5 <<- "white"
        MGvar$activeplot.cex <<- 0.6
        MGvar$activeplot.cex.T1 <<- 0.6
        MGvar$activeplot.cex.T2 <<- 0.6
        MGvar$activeplot.cex.T3 <<- 0.6
        MGvar$activeplot.cex.T4 <<- 0.6
        MGvar$activeplot.cex.T5 <<- 0.6
        MGvar$activeplot.labs <<- "yes"
        MGvar$activeplot.labs.T1 <<- "yes"
        MGvar$activeplot.labs.T2 <<- "yes"
        MGvar$activeplot.labs.T3 <<- "yes"
        MGvar$activeplot.labs.T4 <<- "yes"
        MGvar$activeplot.labs.T5 <<- "yes"
        MGvar$activeplot.showpoints <<- "no"
        MGvar$activeplot.showpoints.T1 <<- "no"
        MGvar$activeplot.showpoints.T2 <<- "no"
        MGvar$activeplot.showpoints.T3 <<- "no"
        MGvar$activeplot.showpoints.T4 <<- "no"
        MGvar$activeplot.showpoints.T5 <<- "no"
        MGvar$activeplot.pointcol <<- "black"
        MGvar$activeplot.pointcol.temp <<- "black"
        MGvar$activeplot.pointcol.T1 <<- "black"
        MGvar$activeplot.pointcol.T2 <<- "black"
        MGvar$activeplot.pointcol.T3 <<- "black"
        MGvar$activeplot.pointcol.T4 <<- "black"
        MGvar$activeplot.pointcol.T5 <<- "black"
        MGvar$activeplot.type <<- 1
        MGvar$activeplot.type.T1 <<- 1
        MGvar$activeplot.type.T2 <<- 1
        MGvar$activeplot.type.T3 <<- 1
        MGvar$activeplot.type.T4 <<- 1
        MGvar$activeplot.type.T5 <<- 1
        MGvar$activeplot.xaxt <<- "n"
        MGvar$activeplot.xaxt.T1 <<- "n"
        MGvar$activeplot.xaxt.T2 <<- "n"
        MGvar$activeplot.xaxt.T3 <<- "n"
        MGvar$activeplot.xaxt.T4 <<- "n"
        MGvar$activeplot.xaxt.T5 <<- "n"
        MGvar$activeplot.yaxt <<- "n"
        MGvar$activeplot.yaxt.T1 <<- "n"
        MGvar$activeplot.yaxt.T2 <<- "n"
        MGvar$activeplot.yaxt.T3 <<- "n"
        MGvar$activeplot.yaxt.T4 <<- "n"
        MGvar$activeplot.yaxt.T5 <<- "n"
        MGvar$activeplot.axescol <<- "black"
        MGvar$activeplot.axescol.temp <<- "black"
        MGvar$activeplot.axescol.T1 <<- "black"
        MGvar$activeplot.axescol.T2 <<- "black"
        MGvar$activeplot.axescol.T3 <<- "black"
        MGvar$activeplot.axescol.T4 <<- "black"
        MGvar$activeplot.axescol.T5 <<- "black"
        MGvar$activeplot.showreg <<- "no"
        MGvar$activeplot.showreg.T1 <<- "no"
        MGvar$activeplot.showreg.T2 <<- "no"
        MGvar$activeplot.showreg.T3 <<- "no"
        MGvar$activeplot.showreg.T4 <<- "no"
        MGvar$activeplot.showreg.T5 <<- "no"
        MGvar$activeplot.regcol <<- "red"
        MGvar$activeplot.regcol.temp <<- "red"
        MGvar$activeplot.regcol.T1 <<- "red"
        MGvar$activeplot.regcol.T2 <<- "red"
        MGvar$activeplot.regcol.T3 <<- "red"
        MGvar$activeplot.regcol.T4 <<- "red"
        MGvar$activeplot.regcol.T5 <<- "red"
        MGvar$activeplot.showleg <<- "no"
        MGvar$activeplot.showleg.T1 <<- "no"
        MGvar$activeplot.showleg.T2 <<- "no"
        MGvar$activeplot.showleg.T3 <<- "no"
        MGvar$activeplot.showleg.T4 <<- "no"
        MGvar$activeplot.showleg.T5 <<- "no"
        MGvar$activeplot.showdist <<- "yes"
        MGvar$activeplot.showdist.T1 <<- "yes"
        MGvar$activeplot.showdist.T2 <<- "yes"
        MGvar$activeplot.showdist.T3 <<- "yes"
        MGvar$activeplot.showdist.T4 <<- "yes"
        MGvar$activeplot.showdist.T5 <<- "yes"
        MGvar$activeplot.distcol <<- "red"
        MGvar$activeplot.distcol.temp <<- "red"
        MGvar$activeplot.distcol.T1 <<- "red"
        MGvar$activeplot.distcol.T2 <<- "red"
        MGvar$activeplot.distcol.T3 <<- "red"
        MGvar$activeplot.distcol.T4 <<- "red"
        MGvar$activeplot.distcol.T5 <<- "red"
        MGvar$Tab1.zoomedswitch <<- "off"
        MGvar$Tab2.zoomedswitch <<- "off"
        MGvar$Tab3.zoomedswitch <<- "off"
        MGvar$Tab4.zoomedswitch <<- "off"
        MGvar$Tab5.zoomedswitch <<- "off"
        MGvar$newCoords <<- as.matrix(0)
        MGvar$newCoords.T1 <<- as.matrix(0)
        MGvar$newCoords.T2 <<- as.matrix(0)
        MGvar$newCoords.T3 <<- as.matrix(0)
        MGvar$newCoords.T4 <<- as.matrix(0)
        MGvar$newCoords.T5 <<- as.matrix(0)
        MGvar$indexLabeled <<- c()
        MGvar$indexLabeled.T1 <<- c()
        MGvar$indexLabeled.T2 <<- c()
        MGvar$indexLabeled.T3 <<- c()
        MGvar$indexLabeled.T4 <<- c()
        MGvar$indexLabeled.T5 <<- c()
        MGvar$labeledPoints <<- list()
        MGvar$MDSmat.Cols <<- as.vector(0)
        MGvar$MDSmat.Cols.T1 <<- as.vector(0)
        MGvar$MDSmat.Cols.T2 <<- as.vector(0)
        MGvar$MDSmat.Cols.T3 <<- as.vector(0)
        MGvar$MDSmat.Cols.T4 <<- as.vector(0)
        MGvar$MDSmat.Cols.T5 <<- as.vector(0)
        MGvar$EnActivePlot.switch.T1 <<- "off"
        MGvar$EnActivePlot.switch.T2 <<- "off"
        MGvar$EnActivePlot.switch.T3 <<- "off"
        MGvar$EnActivePlot.switch.T4 <<- "off"
        MGvar$EnActivePlot.switch.T5 <<- "off"
        MGvar$sthreeDplot.title <<- ""
        MGvar$sthreeDplot.title.show <<- "yes"
        MGvar$sthreeDplot.distmeas <<- "yes"
        MGvar$sthreeDplot.xlab <<- paste("Dim", MGvar$PlottingDimX.3D)
        MGvar$sthreeDplot.ylab <<- paste("Dim", MGvar$PlottingDimY.3D)
        MGvar$sthreeDplot.zlab <<- paste("Dim", MGvar$PlottingDimZ.3D)
        MGvar$sthreeDplot.leg.show <<- "yes"
        MGvar$sthreeDplot.bg <<- "white"
        MGvar$sthreeDplot.bg.temp <<- "white"
        MGvar$sthreeDplot.showpoints <<- "no"
        MGvar$sthreeDplot.showlabels <<- "yes"
        MGvar$sthreeDplot.cex <<- 0.6
        MGvar$sthreeDplot.type <<- 1
        MGvar$sthreeDplot.pointcol <<- "black"
        MGvar$sthreeDplot.pointcol.temp <<- "black"
        MGvar$sthreeDplot.showregX <<- "no"
        MGvar$sthreeDplot.regXcol <<- "red"
        MGvar$sthreeDplot.regXcol.temp <<- "red"
        MGvar$sthreeDplot.showregY <<- "no"
        MGvar$sthreeDplot.regYcol <<- "blue"
        MGvar$sthreeDplot.showregZ <<- "no"
        MGvar$sthreeDplot.regZcol <<- "green"
        MGvar$sthreeDplot.showaxes <<- "yes"
        MGvar$sthreeDplot.axescol <<- "black"
        MGvar$sthreeDplot.axescol.temp <<- "black"
        MGvar$sthreeDplot.showgrid <<- "yes"
        MGvar$sthreeDplot.angle <<- 40
        MGvar$sthreeDplot.HL <<- "no"
        MGvar$sthreeDplot.yscale <<- 1
        MGvar$rglplot.title <<- ""
        MGvar$rglplot.title.show <<- "yes"
        MGvar$rglplot.xlab <<- paste("Dim", MGvar$PlottingDimX.3D)
        MGvar$rglplot.ylab <<- paste("Dim", MGvar$PlottingDimY.3D)
        MGvar$rglplot.zlab <<- paste("Dim", MGvar$PlottingDimZ.3D)
        MGvar$rglplot.showpoints <<- "no"
        MGvar$rglplot.showlabels <<- "yes"
        MGvar$rglplot.ptcol <<- "black"
        MGvar$rglplot.ptcol.temp <<- "black"
        MGvar$rglplot.ptsize <<- "2"
        MGvar$rglplot.axmeas <<- "no"
        MGvar$shep.firstrun <<- "yes"
        MGvar$shepplot.title.show <<- "yes"
        MGvar$shepplot.labs.show <<- "yes"
        MGvar$shepplot.leg.show <<- "no"
        MGvar$shepplot.bg <<- "white"
        MGvar$shepplot.bg.temp <<- "white"
        MGvar$shepplot.showpoints <<- "yes"
        MGvar$shepplot.showlabels <<- "no"
        MGvar$shepplot.cex <<- 0.6
        MGvar$shepplot.type <<- 1
        MGvar$shepplot.pointcol.temp <<- "black"
        MGvar$shepplot.pointcol <<- "black"
        MGvar$shepplot.curve.show <<- "yes"
        MGvar$shepplot.curve.type <<- 1
        MGvar$shepplot.curvecol.temp <<- "red"
        MGvar$shepplot.curvecol <<- "red"
        MGvar$shepplot.Axes.xaxt <<- "s"
        MGvar$shepplot.Axes.yaxt <<- "s"
        MGvar$shepplot.Axescol.temp <<- "black"
        MGvar$shepplot.Axescol <<- "black"
        MGvar$shepplot.bold <<- FALSE
        MGvar$Shep.indexLabeled <<- c()
        MGvar$Shepx <<- as.vector(0)
        MGvar$Shepy <<- as.vector(0)
        MGvar$tShepx <<- as.vector(0)
        MGvar$screeplot.title.show <<- "yes"
        MGvar$screeplot.labs.show <<- "yes"
        MGvar$screeplot.leg.show <<- "yes"
        MGvar$screeplot.bg <<- "white"
        MGvar$screeplot.bg.temp <<- "white"
        MGvar$screeplot.points.show <<- "no"
        MGvar$screeplot.Cdim.show <<- "yes"
        MGvar$screeplot.Odim.show <<- "yes"
        MGvar$screeplot.Ccol <<- "red"
        MGvar$screeplot.Ccol.temp <<- "red"
        MGvar$screeplot.Ocol <<- "blue"
        MGvar$screeplot.Ocol.temp <<- "blue"
        MGvar$screeplot.curve.show <<- "yes"
        MGvar$screeplot.curve.type <<- 1
        MGvar$screeplot.curvecol <<- "black"
        MGvar$screeplot.curvecol.temp <<- "black"
        MGvar$screeplot.Cline.show <<- "yes"
        MGvar$screeplot.Oline.show <<- "yes"
        MGvar$screeplot.Axes.xaxt <<- "s"
        MGvar$screeplot.Axes.yaxt <<- "s"
        MGvar$screeplot.Axescol <<- "black"
        MGvar$screeplot.Axescol.temp <<- "black"
        MGvar$stressplot.title.show <<- "yes"
        MGvar$stressplot.time.show <<- "yes"
        MGvar$stressplot.labs.show <<- "yes"
        MGvar$stressplot.bg <<- "white"
        MGvar$stressplot.bg.temp <<- "white"
        MGvar$stressplot.curve.type <<- 1
        MGvar$stressplot.curvecol <<- "black"
        MGvar$stressplot.curvecol.temp <<- "black"
        MGvar$stressplot.Axes.xaxt <<- "s"
        MGvar$stressplot.Axes.yaxt <<- "s"
        MGvar$stressplot.Axescol <<- "black"
        MGvar$stressplot.Axescol.temp <<- "black"
        MGvar$zoominrat <<- 1
        MGvar$zoominrat.T1 <<- 1
        MGvar$zoominrat.T2 <<- 1
        MGvar$zoominrat.T3 <<- 1
        MGvar$zoominrat.T4 <<- 1
        MGvar$zoominrat.T5 <<- 1
        MGvar$moveup <<- 1
        MGvar$moveup.T1 <<- 1
        MGvar$moveup.T2 <<- 1
        MGvar$moveup.T3 <<- 1
        MGvar$moveup.T4 <<- 1
        MGvar$moveup.T5 <<- 1
        MGvar$movedown <<- 1
        MGvar$movedown.T1 <<- 1
        MGvar$movedown.T2 <<- 1
        MGvar$movedown.T3 <<- 1
        MGvar$movedown.T4 <<- 1
        MGvar$movedown.T5 <<- 1
        MGvar$moveleft <<- 1
        MGvar$moveleft.T1 <<- 1
        MGvar$moveleft.T2 <<- 1
        MGvar$moveleft.T3 <<- 1
        MGvar$moveleft.T4 <<- 1
        MGvar$moveleft.T5 <<- 1
        MGvar$moveright <<- 1
        MGvar$moveright.T1 <<- 1
        MGvar$moveright.T2 <<- 1
        MGvar$moveright.T3 <<- 1
        MGvar$moveright.T4 <<- 1
        MGvar$moveright.T5 <<- 1
        MGvar$EStat.switch <<- "off"
        MGvar$EnStress.switch <<- "off"
        MGvar$EnStress2.switch <<- "off"
        MGvar$procplot.title.show <<- "yes"
        MGvar$procplot.leg.show <<- "no"
        MGvar$procplot.ylab <<- ""
        MGvar$procplot.xlab <<- ""
        MGvar$procplot.bg.temp <<- "white"
        MGvar$procplot.bg <<- "white"
        MGvar$procplot.showpoints1 <<- "no"
        MGvar$procplot.labs1 <<- "yes"
        MGvar$procplot.cex1 <<- 0.7
        MGvar$procplot.type1 <<- 1
        MGvar$procplot.point1col.temp <<- "red"
        MGvar$procplot.point1col <<- "red"
        MGvar$procplot.showpoints2 <<- "no"
        MGvar$procplot.labs2 <<- "yes"
        MGvar$procplot.cex2 <<- 0.7
        MGvar$procplot.type2 <<- 1
        MGvar$procplot.point2col.temp <<- "blue"
        MGvar$procplot.point2col <<- "blue"
        MGvar$procplot.yaxt <<- "n"
        MGvar$procplot.xaxt <<- "n"
        MGvar$procplot.axescol.temp <<- "black"
        MGvar$procplot.axescol <<- "black"
        MGvar$procplot.showreg1 <<- "no"
        MGvar$procplot.showreg2 <<- "no"
        MGvar$procplot.regcol1 <<- "red"
        MGvar$procplot.regcol1.temp <<- "red"
        MGvar$procplot.regcol2 <<- "blue"
        MGvar$procplot.regcol2.temp <<- "blue"
        MGvar$Proc.indexLabeled <<- c()
        MGvar$EnProcPlot.switch <<- "off"
        MGvar$proc.zoominrat <<- 1
        MGvar$proc.moveup <<- 1
        MGvar$proc.movedown <<- 1
        MGvar$proc.moveleft <<- 1
        MGvar$proc.moveright <<- 1
        MGvar$Conf.Main.val <<- tclVar("1")
        MGvar$Conf.Leg.val <<- tclVar("0")
        MGvar$Conf.Dist.val <<- tclVar("1")
        MGvar$Conf.Ylab.val <<- tclVar("0")
        MGvar$Conf.Xlab.val <<- tclVar("0")
        if (MGvar$GenSet.DefPts == "yes") {
            MGvar$Conf.Points.val <<- tclVar("0")
        }
        if (MGvar$GenSet.DefPts == "no") {
            MGvar$Conf.Points.val <<- tclVar("0")
        }
        if (MGvar$GenSet.DefLabs == "yes") {
            MGvar$Conf.Labels.val <<- tclVar("1")
        }
        if (MGvar$GenSet.DefLabs == "no") {
            MGvar$Conf.Labels.val <<- tclVar("0")
        }
        MGvar$Conf.PT.var <<- tclVar("Empty Circles")
        MGvar$Conf.AxesMeas.val <<- tclVar("0")
        if (MGvar$GenSet.DefReg == "yes") {
            MGvar$Conf.RegLine.val <<- tclVar("1")
        }
        if (MGvar$GenSet.DefReg == "no") {
            MGvar$Conf.RegLine.val <<- tclVar("0")
        }
        MGvar$Conf.DistLine.val <<- tclVar("1")
        MGvar$Shep.Main.val <<- tclVar("1")
        MGvar$Shep.Lab.val <<- tclVar("1")
        MGvar$Shep.Leg.val <<- tclVar("0")
        MGvar$Shep.Points.val <<- tclVar("1")
        MGvar$Shep.Labels.val <<- tclVar("0")
        MGvar$Shep.PS.var <<- tclVar(MGvar$shepplot.cex)
        MGvar$Shep.PT.var <<- tclVar("Empty Circles")
        MGvar$Shep.Curve.val <<- tclVar("1")
        MGvar$Shep.LineT.val <<- tclVar("Solid Line")
        MGvar$Shep.Bold.val <<- tclVar("0")
        MGvar$Shep.Axes.val <<- tclVar("1")
        MGvar$EnShep.switch <<- "off"
        MGvar$shep.initpos <<- TRUE
        MGvar$latest.xPlotCoord <<- 0
        MGvar$latest.yPlotCoord <<- 0
        MGvar$first.xPlotCoord <<- 0
        MGvar$first.yPlotCoord <<- 0
        MGvar$Scree.Main.val <<- tclVar("1")
        MGvar$Scree.Lab.val <<- tclVar("1")
        MGvar$Scree.Leg.val <<- tclVar("1")
        MGvar$Scree.Points.val <<- tclVar("0")
        MGvar$Scree.CDim.val <<- tclVar("1")
        MGvar$Scree.ODim.val <<- tclVar("1")
        MGvar$Scree.Curve.val <<- tclVar("1")
        MGvar$Scree.LineT.val <<- tclVar("Solid Line")
        MGvar$Scree.CLine.val <<- tclVar("1")
        MGvar$Scree.OLine.val <<- tclVar("1")
        MGvar$Scree.Axes.val <<- tclVar("1")
        MGvar$Stress.Main.val <<- tclVar("1")
        MGvar$Stress.Lab.val <<- tclVar("1")
        MGvar$Stress.Time.val <<- tclVar("1")
        MGvar$Stress.LineT.val <<- tclVar("Solid Line")
        MGvar$Stress.Axes.val <<- tclVar("1")
        MGvar$Zoom.Leg.val <<- tclVar("0")
        MGvar$Zoom.RegLine.val <<- tclVar("0")
        MGvar$Proc.Main.val <<- tclVar("1")
        MGvar$Proc.Leg.val <<- tclVar("0")
        MGvar$Proc.Ylab.val <<- tclVar("0")
        MGvar$Proc.Xlab.val <<- tclVar("0")
        MGvar$Proc.Points1.val <<- tclVar("0")
        MGvar$Proc.Labels1.val <<- tclVar("1")
        MGvar$Proc.PT1.val <<- tclVar("Empty Circles")
        MGvar$Proc.Points2.val <<- tclVar("0")
        MGvar$Proc.Labels2.val <<- tclVar("1")
        MGvar$Proc.PT2.val <<- tclVar("Empty Circles")
        MGvar$Proc.AxesMeas.val <<- tclVar("0")
        MGvar$Proc.RegLine1.val <<- tclVar("0")
        MGvar$Proc.RegLine2.val <<- tclVar("0")
        MGvar$Stat.Main.val <<- tclVar("1")
        MGvar$Stat.Dist.val <<- tclVar("1")
        MGvar$Stat.Leg.val <<- tclVar("0")
        MGvar$Stat.Labs.val <<- tclVar("1")
        MGvar$Stat.Points.val <<- tclVar("0")
        MGvar$Stat.Labels.val <<- tclVar("1")
        MGvar$Stat.PT.var <<- tclVar("Empty Circles")
        MGvar$Stat.HL.var <<- tclVar("0")
        MGvar$Stat.RegLine.val <<- tclVar("0")
        MGvar$Stat.AxesMeas.val <<- tclVar("1")
        MGvar$Stat.Grid.val <<- tclVar("1")
        MGvar$RGL.Main.val <<- tclVar("1")
        MGvar$RGL.axlabs.val <<- tclVar("1")
        MGvar$RGL.points.val <<- tclVar("0")
        MGvar$RGL.ptlabels.val <<- tclVar("1")
        MGvar$RGL.AxesMeas.val <<- tclVar("1")
        MGvar$datnam <<- tclVar("Current name is ")
        MGvar$zoomedplot.showreg <<- "no"
        MGvar$zoomedplot.regcol <<- "red"
        MGvar$zoomedplot.regcol.temp <<- "red"
        MGvar$indexZLabeled <<- c()
        MGvar$labeledZPoints <<- list()
        MGvar$indexZMLabeled <<- c()
        MGvar$seczoomswitch <<- "off"
        MGvar$GS.MDSConfig.val <<- tclVar("1")
        MGvar$GS.Shep.val <<- tclVar("1")
        MGvar$GS.Stress.val <<- tclVar("1")
        MGvar$GS.Scree.val <<- tclVar("1")
        MGvar$GS.Win.val <<- tclVar("0")
        MGvar$GS.Lab.val <<- tclVar("1")
        MGvar$GS.Pt.val <<- tclVar("0")
        MGvar$GS.Reg.val <<- tclVar("0")
        MGvar$GS.ShepP.val <<- tclVar("1")
        MGvar$GS.Pcol.val <<- tclVar("0")
        MGvar$GS.IL.val <<- tclVar("1")
        MGvar$GS.IP.val <<- tclVar("Top")
        MGvar$GS.UpConf.val <<- tclVar("1")
        MGvar$GS.UpShep.val <<- tclVar("1")
        MGvar$GS.UpStress.val <<- tclVar("1")
        MGvar$GS.UpProg.val <<- tclVar("1")
        MGvar$parPlotSize <<- ""
        MGvar$usrCoords <<- ""
        MGvar$labelsVec <<- ""
        MGvar$ZlabelsVec <<- ""
        MGvar$OrderedVals <<- ""
        MGvar$originaldistmat <<- as.matrix(0)
        MGvar$stressval <<- c()
        MGvar$X.Coords <<- ""
        MGvar$Y.Coords <<- ""
        MGvar$Shep.labelsVec <<- ""
        MGvar$Shep.parPlotSize <<- ""
        MGvar$Shep.usrCoords <<- ""
        MGvar$Proc.labelsVec <<- ""
        MGvar$Proc.parPlotSize <<- ""
        MGvar$Proc.usrCoords <<- ""
        MGvar$relIndex <<- ""
        MGvar$RelXCoords <<- ""
        MGvar$RelYCoords <<- ""
        MGvar$diffvec <<- ""
    }
    MGcomp <- list()
    Initialise.GUIcomps <- function() {
        MGcomp$Confarray <<- c()
        MGcomp$RemAxarray <<- c()
        MGcomp$Remparray <<- c()
        MGcomp$table.Conf <<- ""
    }
    Initialise.GUIcomps()
    FirstPointColInitialise <- function() {
        for (i in 1:nrow(MGvar$activedata)) {
            MGvar$MDSmat.Cols.T1[i] <<- MGvar$activeplot.pointcol.T1
            MGvar$MDSmat.Cols.T2[i] <<- MGvar$activeplot.pointcol.T2
            MGvar$MDSmat.Cols.T3[i] <<- MGvar$activeplot.pointcol.T3
            MGvar$MDSmat.Cols.T4[i] <<- MGvar$activeplot.pointcol.T4
            MGvar$MDSmat.Cols.T5[i] <<- MGvar$activeplot.pointcol.T5
        }
    }
    PointColInitialise <- function() {
        for (i in 1:nrow(MGvar$activedata)) {
            MGvar$MDSmat.Cols[i] <<- MGvar$activeplot.pointcol
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            for (i in 1:nrow(MGvar$activedata)) {
                MGvar$MDSmat.Cols.T1[i] <<- MGvar$activeplot.pointcol.T1
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            for (i in 1:nrow(MGvar$activedata)) {
                MGvar$MDSmat.Cols.T2[i] <<- MGvar$activeplot.pointcol.T2
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            for (i in 1:nrow(MGvar$activedata)) {
                MGvar$MDSmat.Cols.T3[i] <<- MGvar$activeplot.pointcol.T3
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            for (i in 1:nrow(MGvar$activedata)) {
                MGvar$MDSmat.Cols.T4[i] <<- MGvar$activeplot.pointcol.T4
            }
        }
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            for (i in 1:nrow(MGvar$activedata)) {
                MGvar$MDSmat.Cols.T5[i] <<- MGvar$activeplot.pointcol.T5
            }
        }
    }
    versioncheck <- function() {
        packs = installed.packages()
        for (i in 1:nrow(packs)) {
            if (rownames(packs)[i] == "tcltk2") {
                ver = packs[i, 3]
                if (ver != "1.1-5") {
                  mes = paste("The MDS-GUI was developed using tcltk2 version 1.1-5. You currently have tcltk2 version ", 
                    ver, " installed. Keeping with your currently installed version may effect the look of the MDS-GUI. It is highly recommened that you install tcltk2 1.1-5.")
                  tkmessageBox(message = mes, icon = "error", 
                    type = "ok")
                }
            }
        }
    }
    mytt <- tktoplevel()
    topMenu <- tkmenu(mytt)
    tkconfigure(mytt, menu = topMenu)
    tkwm.title(mytt, "MDS-GUI")
    GUI.AvailableScreenWidth <- round(as.numeric(tkwinfo("screenwidth", 
        mytt)))
    GUI.AvailableScreenHeight <- round(as.numeric(tkwinfo("screenheight", 
        mytt)))
    if (GUI.AvailableScreenWidth/GUI.AvailableScreenHeight <= 
        1080/660) {
        GUI.ScreenWidth <- min(1080, round(GUI.AvailableScreenWidth * 
            0.9))
        GUI.ScreenHeight <- min(660, round(GUI.AvailableScreenWidth/1080 * 
            660 * 0.9))
    }
    else {
        GUI.ScreenWidth <- min(1080, round(GUI.AvailableScreenHeight/660 * 
            1080 * 0.9))
        GUI.ScreenHeight <- min(660, round(GUI.AvailableScreenHeight * 
            0.9))
    }
    .Tcl(paste("wm geometry ", mytt, " ", GUI.ScreenWidth, "x", 
        GUI.ScreenHeight, "+", round(GUI.AvailableScreenWidth/2 - 
            GUI.ScreenWidth/2, 0), "+", round(GUI.AvailableScreenHeight/2 - 
            GUI.ScreenHeight/2, 0), sep = ""))
    .Tcl(paste("wm minsize ", mytt, " 900 600", sep = ""))
    tkwm.resizable(mytt, "0", "0")
    tkwm.deiconify(mytt)
    col.prim <- "#dae0f1"
    col.prim.temp <- "#dae0f1"
    col.sec <- "#dae0f1"
    col.sec.temp <- "#dae0f1"
    col.frame <- "#9ed3e4"
    backcanvas <- tkcanvas(mytt, width = "1128", height = "756", 
        bg = col.prim)
    tkplace(backcanvas, relx = 0, rely = 0, relheight = 1, relwidth = 1, 
        `in` = mytt)
    Initialise.Variables()
    fontHeading <- tkfont.create(family = "times", size = 24, 
        weight = "bold", slant = "italic")
    tkplace(tklabel(mytt, text = "The MDS-GUI: Version 0.2", 
        foreground = "black", font = fontHeading, background = col.prim), 
        relx = 0.535, rely = 0.01, `in` = mytt)
    UNHeading.Font <- tkfont.create(family = "times", size = 16, 
        weight = "bold")
    tclvalue(MGvar$UserName) <- paste("User: ", tclvalue(MGvar$Active.UserName))
    tkplace(tklabel(mytt, text = tclvalue(MGvar$UserName), textvariable = MGvar$UserName, 
        font = UNHeading.Font, background = col.prim), relx = 0.01, 
        rely = 0.005, `in` = mytt)
    Myhscale <- 1.52
    Myvscale <- 1.36
    myPlottingNB <- tk2notebook(mytt, tabs = NULL)
    tb2.1 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.1, text = "Plot1")
    framePlot1 <- tkframe(tb2.1, relief = "groove", borderwidth = 2, 
        background = "black")
    img <- tkrplot(framePlot1, function() plotting2D(MGvar$MDSmat), 
        hscale = Myhscale, vscale = Myvscale)
    tkplace(img, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot1)
    tkplace(framePlot1, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.1)
    tkbind(img, "<Button-1>", OnPlotLeftClick)
    tkconfigure(img, cursor = "hand2")
    tkbind(img, "<B1-Motion>", BrushingPointMove)
    tkbind(img, "<ButtonRelease-1>", OnRelease.Main)
    tb2.2 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.2, text = "Plot2")
    framePlot2 <- tkframe(tb2.2, relief = "groove", borderwidth = 2, 
        background = "black")
    img2 <- tkrplot(framePlot2, function() plotting2D(MGvar$MDSmat), 
        hscale = Myhscale, vscale = Myvscale)
    tkplace(img2, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot2)
    tkplace(framePlot2, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.2)
    tkbind(img2, "<Button-1>", OnPlotLeftClick)
    tkconfigure(img2, cursor = "arrow")
    tkconfigure(img2, state = "disabled")
    tb2.3 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.3, text = "Plot3")
    framePlot3 <- tkframe(tb2.3, relief = "groove", borderwidth = 2, 
        background = "black")
    img3 <- tkrplot(framePlot3, function() plotting2D(MGvar$MDSmat), 
        hscale = Myhscale, vscale = Myvscale)
    tkplace(img3, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot3)
    tkplace(framePlot3, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.3)
    tkbind(img3, "<Button-1>", OnPlotLeftClick)
    tkconfigure(img3, cursor = "arrow")
    tkconfigure(img3, state = "disabled")
    tb2.4 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.4, text = "Plot4")
    framePlot4 <- tkframe(tb2.4, relief = "groove", borderwidth = 2, 
        background = "black")
    img4 <- tkrplot(framePlot4, function() plotting2D(MGvar$MDSmat), 
        hscale = Myhscale, vscale = Myvscale)
    tkplace(img4, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot4)
    tkplace(framePlot4, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.4)
    tkbind(img4, "<Button-1>", OnPlotLeftClick)
    tkconfigure(img4, cursor = "arrow")
    tkconfigure(img4, state = "disabled")
    tb2.5 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.5, text = "Plot5")
    framePlot5 <- tkframe(tb2.5, relief = "groove", borderwidth = 2, 
        background = "black")
    img5 <- tkrplot(framePlot5, function() plotting2D(MGvar$MDSmat), 
        hscale = Myhscale, vscale = Myvscale)
    tkplace(img5, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot5)
    tkplace(framePlot5, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.5)
    tkbind(img5, "<Button-1>", OnPlotLeftClick)
    tkconfigure(img5, cursor = "arrow")
    tkconfigure(img5, state = "disabled")
    tb2.6 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.6, text = "Static 3D Plot")
    framePlot6 <- tkframe(tb2.6, relief = "groove", borderwidth = 2, 
        background = "black")
    img3Dstat <- tkrplot(framePlot6, function() plotting3Dstatic(MGvar$MDSmat), 
        hscale = Myhscale, vscale = Myvscale)
    tkplace(img3Dstat, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot6)
    tkplace(framePlot6, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.6)
    tkconfigure(img3Dstat, cursor = "arrow")
    ProcConf1 = MGvar$MDSmat.T1
    ProcConf2 = MGvar$MDSmat.T2
    tb2.7 <- tk2frame(myPlottingNB)
    tkadd(myPlottingNB, tb2.7, text = "Procrustes")
    framePlot7 <- tkframe(tb2.7, relief = "groove", borderwidth = 2, 
        background = "black")
    procimg <- tkrplot(framePlot7, function() ProcrustesMDSGUIplot(conf1 = ProcConf1, 
        conf2 = ProcConf2), hscale = Myhscale, vscale = Myvscale)
    tkplace(procimg, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = framePlot7)
    tkplace(framePlot7, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb2.7)
    tkbind(procimg, "<Button-1>", Proc.OnLeftClick)
    tkconfigure(procimg, cursor = "hand2")
    ClearProcPoints <- function() {
        MGvar$Proc.indexLabeled <<- c()
        tkrreplot(procimg)
        if (MGvar$EnProcPlot.switch == "on") {
            tkrreplot(MGcomp$POprocimg)
        }
    }
    ProcPlotMenu <- tkmenu(procimg, tearoff = FALSE)
    tkadd(ProcPlotMenu, "command", label = "Label Specific Point", 
        command = Proc.LabelSpecificPoint)
    tkadd(ProcPlotMenu, "command", label = "Clear Added Point Labels", 
        command = ClearProcPoints)
    tkadd(ProcPlotMenu, "separator")
    tkadd(ProcPlotMenu, "command", label = "Pop-out Plot", command = EnlargedProcPlot)
    tkadd(ProcPlotMenu, "command", label = "Copy to Clipboard", 
        command = function() {
            tkrreplot(procimg)
        })
    tkadd(ProcPlotMenu, "separator")
    tkadd(ProcPlotMenu, "command", label = "Procrustes Plot Options", 
        command = ProcPlotOps)
    RightClickProc <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", procimg))
        rooty <- as.integer(tkwinfo("rooty", procimg))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", ProcPlotMenu, xTxt, yTxt)
    }
    tkbind(procimg, "<Button-3>", RightClickProc)
    tkplace(tkbutton(mytt, text = "+", width = 2, command = function() {
        MGvar$proc.zoominrat <<- MGvar$proc.zoominrat/1.1
        tkrreplot(procimg)
    }), relx = 0.001, rely = 0.105, `in` = procimg)
    tkplace(tkbutton(mytt, text = "-", width = 2, command = function() {
        MGvar$proc.zoominrat <<- MGvar$proc.zoominrat * 1.1
        tkrreplot(procimg)
    }), relx = 0.001, rely = 0.15, `in` = procimg)
    tkplace(tkbutton(mytt, text = "^", width = 2, command = function() {
        MGvar$proc.moveup <<- MGvar$proc.moveup * 1.05
        tkrreplot(procimg)
    }), relx = 0.04, rely = 0.865, `in` = procimg)
    tkplace(tkbutton(mytt, text = "v", width = 2, command = function() {
        MGvar$proc.movedown <<- MGvar$proc.movedown * 1.05
        tkrreplot(procimg)
    }), relx = 0.04, rely = 0.955, `in` = procimg)
    tkplace(tkbutton(mytt, text = "<", width = 2, command = function() {
        MGvar$proc.moveleft <<- MGvar$proc.moveleft * 1.05
        tkrreplot(procimg)
    }), relx = 0.001, rely = 0.91, `in` = procimg)
    tkplace(tkbutton(mytt, text = ">", width = 2, command = function() {
        MGvar$proc.moveright <<- MGvar$proc.moveright * 1.05
        tkrreplot(procimg)
    }), relx = 0.078, rely = 0.91, `in` = procimg)
    tkplace(tkbutton(mytt, text = "C", width = 2, command = function() proc.origpos()), 
        relx = 0.04, rely = 0.91, `in` = procimg)
    ClearMainPoints1 <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab1") {
            MGvar$indexLabeled <<- c()
            MGvar$indexLabeled.T1 <<- c()
            MGvar$labeledPoints <<- list()
            tabplot()
        }
        else {
            tkmessageBox(message = "This is not the active plot", 
                icon = "error")
        }
    }
    fontMenuHeading <- tkfont.create(family = "times", size = 10, 
        weight = "bold", slant = "italic")
    MainPlotMenu1 <- tkmenu(img, tearoff = FALSE)
    tkadd(MainPlotMenu1, "command", label = "                        Plot 1 Menu", 
        font = fontMenuHeading)
    tkadd(MainPlotMenu1, "separator")
    tkadd(MainPlotMenu1, "command", label = "Label Specific Point", 
        command = LabelSpecificPoint)
    tkadd(MainPlotMenu1, "command", label = "Clear Added Point Labels", 
        command = ClearMainPoints1)
    tkadd(MainPlotMenu1, "separator")
    tkadd(MainPlotMenu1, "command", label = "Relocate Group of Points", 
        command = MoveButtonInitiate)
    tkadd(MainPlotMenu1, "command", label = "Remove a Point", 
        command = RemovePoint)
    tkadd(MainPlotMenu1, "command", label = "Use Coordinates as Starting Configuration", 
        command = NewStartConfig)
    tkadd(MainPlotMenu1, "separator")
    tkadd(MainPlotMenu1, "command", label = "Rotate and Reflect", 
        command = RotateandReflect)
    tkadd(MainPlotMenu1, "command", label = "Advanced Zoom", 
        command = ZoomPlot)
    tkadd(MainPlotMenu1, "separator")
    tkadd(MainPlotMenu1, "command", label = "Change Point Colour", 
        command = ColButtonInitiate)
    tkadd(MainPlotMenu1, "command", label = "Default Point Colours", 
        command = DefaultPtCols)
    tkadd(MainPlotMenu1, "separator")
    MenuDispAx1 <- tclVar(0)
    tkadd(MainPlotMenu1, "checkbutton", label = "Display Variable Axes", 
        variable = MenuDispAx1, command = ShowregfromMenu)
    tkadd(MainPlotMenu1, "command", label = "Remove Axes of Variable(s)", 
        command = AxesRemove, state = "disabled")
    tkadd(MainPlotMenu1, "separator")
    tkadd(MainPlotMenu1, "command", label = "Pop-out Plot", command = EnlargedActivePlot)
    tkadd(MainPlotMenu1, "command", label = "Copy Plot1 to Clipboard", 
        command = function() {
            tkrreplot(img, function() plotting2D(MGvar$MDSmat.T1, 
                title = MGvar$activeplot.title.T1, Measure = MGvar$dMeas.T1, 
                showtitle = MGvar$activeplot.title.show.T1, showmeas = MGvar$activeplot.distmeas.T1, 
                xlabel = MGvar$activeplot.xlab.T1, ylabel = MGvar$activeplot.ylab.T1, 
                bgcol = MGvar$activeplot.bg.T1, pointcex = MGvar$activeplot.cex.T1, 
                showlabs = MGvar$activeplot.labs.T1, showpoints = MGvar$activeplot.showpoints.T1, 
                pointcol = MGvar$activeplot.pointcol.T1, pointshape = MGvar$activeplot.type.T1, 
                ymeas = MGvar$activeplot.yaxt.T1, xmeas = MGvar$activeplot.xaxt.T1, 
                axcol = MGvar$activeplot.axescol.T1, indexLabeled = MGvar$indexLabeled.T1, 
                zoomedcoords = MGvar$newCoords.T1, Zrat = MGvar$zoominrat.T1, 
                Mvup = MGvar$moveup.T1, Mvdn = MGvar$movedown.T1, 
                Mvlt = MGvar$moveleft.T1, Mvrt = MGvar$moveright.T1))
        })
    tkadd(MainPlotMenu1, "separator")
    tkadd(MainPlotMenu1, "command", label = "Plot Options", command = ConfPlotOptions)
    RightClickMain1 <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", img))
        rooty <- as.integer(tkwinfo("rooty", img))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", MainPlotMenu1, xTxt, yTxt)
    }
    tkbind(img, "<Button-3>", RightClickMain1)
    z1.1 <- tkbutton(mytt, text = "+", width = 2, command = function() simplezoomin())
    tkplace(z1.1, relx = 0.001, rely = 0.105, `in` = img)
    z1.2 <- tkbutton(mytt, text = "-", width = 2, command = function() simplezoomout())
    tkplace(z1.2, relx = 0.001, rely = 0.15, `in` = img)
    z1.3 <- tkbutton(mytt, text = "^", width = 2, command = function() moveupfunc())
    tkplace(z1.3, relx = 0.04, rely = 0.865, `in` = img)
    z1.4 <- tkbutton(mytt, text = "v", width = 2, command = function() movedownfunc())
    tkplace(z1.4, relx = 0.04, rely = 0.955, `in` = img)
    z1.5 <- tkbutton(mytt, text = "<", width = 2, command = function() moveleftfunc())
    tkplace(z1.5, relx = 0.001, rely = 0.91, `in` = img)
    z1.6 <- tkbutton(mytt, text = ">", width = 2, command = function() moverightfunc())
    tkplace(z1.6, relx = 0.078, rely = 0.91, `in` = img)
    z1.7 <- tkbutton(mytt, text = "C", width = 2, command = function() origpos())
    tkplace(z1.7, relx = 0.04, rely = 0.91, `in` = img)
    ClearMainPoints2 <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab2") {
            MGvar$indexLabeled <<- c()
            MGvar$indexLabeled.T2 <<- c()
            MGvar$labeledPoints <<- list()
            tabplot()
        }
        else {
            tkmessageBox(message = "This is not the active plot", 
                icon = "error")
        }
    }
    MainPlotMenu2 <- tkmenu(img2, tearoff = FALSE)
    tkadd(MainPlotMenu2, "command", label = "                        Plot 2 Menu", 
        font = fontMenuHeading)
    tkadd(MainPlotMenu2, "separator")
    tkadd(MainPlotMenu2, "command", label = "Label Specific Point", 
        command = LabelSpecificPoint)
    tkadd(MainPlotMenu2, "command", label = "Clear Added Point Labels", 
        command = ClearMainPoints2)
    tkadd(MainPlotMenu2, "separator")
    tkadd(MainPlotMenu2, "command", label = "Relocate Group of Points", 
        command = MoveButtonInitiate)
    tkadd(MainPlotMenu2, "command", label = "Remove a Point", 
        command = RemovePoint)
    tkadd(MainPlotMenu2, "command", label = "Use Coordinates as Starting Configuration", 
        command = NewStartConfig)
    tkadd(MainPlotMenu2, "separator")
    tkadd(MainPlotMenu2, "command", label = "Rotate and Reflect", 
        command = RotateandReflect)
    tkadd(MainPlotMenu2, "command", label = "Advanced Zoom", 
        command = ZoomPlot)
    tkadd(MainPlotMenu2, "separator")
    tkadd(MainPlotMenu2, "command", label = "Change Point Colour", 
        command = ColButtonInitiate)
    tkadd(MainPlotMenu2, "command", label = "Default Point Colours", 
        command = DefaultPtCols)
    tkadd(MainPlotMenu2, "separator")
    MenuDispAx2 <- tclVar(0)
    tkadd(MainPlotMenu2, "checkbutton", label = "Display Variable Axes", 
        variable = MenuDispAx2, command = ShowregfromMenu)
    tkadd(MainPlotMenu2, "command", label = "Remove Axes of Variable(s)", 
        command = AxesRemove, state = "disabled")
    tkadd(MainPlotMenu2, "separator")
    tkadd(MainPlotMenu2, "command", label = "Pop-out Plot", command = EnlargedActivePlot)
    tkadd(MainPlotMenu2, "command", label = "Copy Plot2 to Clipboard", 
        command = function() {
            tkrreplot(img2, function() plotting2D(MGvar$MDSmat.T2, 
                title = MGvar$activeplot.title.T2, Measure = MGvar$dMeas.T2, 
                showtitle = MGvar$activeplot.title.show.T2, showmeas = MGvar$activeplot.distmeas.T2, 
                xlabel = MGvar$activeplot.xlab.T2, ylabel = MGvar$activeplot.ylab.T2, 
                bgcol = MGvar$activeplot.bg.T2, pointcex = MGvar$activeplot.cex.T2, 
                showlabs = MGvar$activeplot.labs.T2, showpoints = MGvar$activeplot.showpoints.T2, 
                pointcol = MGvar$activeplot.pointcol.T2, pointshape = MGvar$activeplot.type.T2, 
                ymeas = MGvar$activeplot.yaxt.T2, xmeas = MGvar$activeplot.xaxt.T2, 
                axcol = MGvar$activeplot.axescol.T2, indexLabeled = MGvar$indexLabeled.T2, 
                zoomedcoords = MGvar$newCoords.T2, Zrat = MGvar$zoominrat.T2, 
                Mvup = MGvar$moveup.T2, Mvdn = MGvar$movedown.T2, 
                Mvlt = MGvar$moveleft.T2, Mvrt = MGvar$moveright.T2))
        })
    tkadd(MainPlotMenu2, "separator")
    tkadd(MainPlotMenu2, "command", label = "Plot Options", command = ConfPlotOptions)
    RightClickMain2 <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", img2))
        rooty <- as.integer(tkwinfo("rooty", img2))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", MainPlotMenu2, xTxt, yTxt)
    }
    tkbind(img2, "<Button-3>", RightClickMain2)
    z2.1 <- tkbutton(mytt, text = "+", width = 2, command = function() simplezoomin())
    tkplace(z2.1, relx = 0.001, rely = 0.105, `in` = img2)
    z2.2 <- tkbutton(mytt, text = "-", width = 2, command = function() simplezoomout())
    tkplace(z2.2, relx = 0.001, rely = 0.15, `in` = img2)
    z2.3 <- tkbutton(mytt, text = "^", width = 2, command = function() moveupfunc())
    tkplace(z2.3, relx = 0.04, rely = 0.865, `in` = img2)
    z2.4 <- tkbutton(mytt, text = "v", width = 2, command = function() movedownfunc())
    tkplace(z2.4, relx = 0.04, rely = 0.955, `in` = img2)
    z2.5 <- tkbutton(mytt, text = "<", width = 2, command = function() moveleftfunc())
    tkplace(z2.5, relx = 0.001, rely = 0.91, `in` = img2)
    z2.6 <- tkbutton(mytt, text = ">", width = 2, command = function() moverightfunc())
    tkplace(z2.6, relx = 0.078, rely = 0.91, `in` = img2)
    z2.7 <- tkbutton(mytt, text = "C", width = 2, command = function() origpos())
    tkplace(z2.7, relx = 0.04, rely = 0.91, `in` = img2)
    ClearMainPoints3 <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab3") {
            MGvar$indexLabeled <<- c()
            MGvar$indexLabeled.T3 <<- c()
            MGvar$labeledPoints <<- list()
            tabplot()
        }
        else {
            tkmessageBox(message = "This is not the active plot", 
                icon = "error")
        }
    }
    MainPlotMenu3 <- tkmenu(img3, tearoff = FALSE)
    tkadd(MainPlotMenu3, "command", label = "                        Plot 3 Menu", 
        font = fontMenuHeading)
    tkadd(MainPlotMenu3, "separator")
    tkadd(MainPlotMenu3, "command", label = "Label Specific Point", 
        command = LabelSpecificPoint)
    tkadd(MainPlotMenu3, "command", label = "Clear Added Point Labels", 
        command = ClearMainPoints3)
    tkadd(MainPlotMenu3, "separator")
    tkadd(MainPlotMenu3, "command", label = "Relocate Group of Points", 
        command = MoveButtonInitiate)
    tkadd(MainPlotMenu3, "command", label = "Remove a Point", 
        command = RemovePoint)
    tkadd(MainPlotMenu3, "command", label = "Use Coordinates as Starting Configuration", 
        command = NewStartConfig)
    tkadd(MainPlotMenu3, "separator")
    tkadd(MainPlotMenu3, "command", label = "Rotate and Reflect", 
        command = RotateandReflect)
    tkadd(MainPlotMenu3, "command", label = "Advanced Zoom", 
        command = ZoomPlot)
    tkadd(MainPlotMenu3, "separator")
    tkadd(MainPlotMenu3, "command", label = "Change Point Colour", 
        command = ColButtonInitiate)
    tkadd(MainPlotMenu3, "command", label = "Default Point Colours", 
        command = DefaultPtCols)
    tkadd(MainPlotMenu3, "separator")
    MenuDispAx3 <- tclVar(0)
    tkadd(MainPlotMenu3, "checkbutton", label = "Display Variable Axes", 
        variable = MenuDispAx3, command = ShowregfromMenu)
    tkadd(MainPlotMenu3, "command", label = "Remove Axes of Variable(s)", 
        command = AxesRemove, state = "disabled")
    tkadd(MainPlotMenu3, "separator")
    tkadd(MainPlotMenu3, "command", label = "Pop-out Plot", command = EnlargedActivePlot)
    tkadd(MainPlotMenu3, "command", label = "Copy Plot3 to Clipboard", 
        command = function() {
            tkrreplot(img3, function() plotting2D(MGvar$MDSmat.T3, 
                title = MGvar$activeplot.title.T3, Measure = MGvar$dMeas.T3, 
                showtitle = MGvar$activeplot.title.show.T3, showmeas = MGvar$activeplot.distmeas.T3, 
                xlabel = MGvar$activeplot.xlab.T3, ylabel = MGvar$activeplot.ylab.T3, 
                bgcol = MGvar$activeplot.bg.T3, pointcex = MGvar$activeplot.cex.T3, 
                showlabs = MGvar$activeplot.labs.T3, showpoints = MGvar$activeplot.showpoints.T3, 
                pointcol = MGvar$activeplot.pointcol.T3, pointshape = MGvar$activeplot.type.T3, 
                ymeas = MGvar$activeplot.yaxt.T3, xmeas = MGvar$activeplot.xaxt.T3, 
                axcol = MGvar$activeplot.axescol.T3, indexLabeled = MGvar$indexLabeled.T3, 
                zoomedcoords = MGvar$newCoords.T3, Zrat = MGvar$zoominrat.T3, 
                Mvup = MGvar$moveup.T3, Mvdn = MGvar$movedown.T3, 
                Mvlt = MGvar$moveleft.T3, Mvrt = MGvar$moveright.T3))
        })
    tkadd(MainPlotMenu3, "separator")
    tkadd(MainPlotMenu3, "command", label = "Plot Options", command = ConfPlotOptions)
    RightClickMain3 <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", img3))
        rooty <- as.integer(tkwinfo("rooty", img3))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", MainPlotMenu3, xTxt, yTxt)
    }
    tkbind(img3, "<Button-3>", RightClickMain3)
    z3.1 <- tkbutton(mytt, text = "+", width = 2, command = function() simplezoomin())
    tkplace(z3.1, relx = 0.001, rely = 0.105, `in` = img3)
    z3.2 <- tkbutton(mytt, text = "-", width = 2, command = function() simplezoomout())
    tkplace(z3.2, relx = 0.001, rely = 0.15, `in` = img3)
    z3.3 <- tkbutton(mytt, text = "^", width = 2, command = function() moveupfunc())
    tkplace(z3.3, relx = 0.04, rely = 0.865, `in` = img3)
    z3.4 <- tkbutton(mytt, text = "v", width = 2, command = function() movedownfunc())
    tkplace(z3.4, relx = 0.04, rely = 0.955, `in` = img3)
    z3.5 <- tkbutton(mytt, text = "<", width = 2, command = function() moveleftfunc())
    tkplace(z3.5, relx = 0.001, rely = 0.91, `in` = img3)
    z3.6 <- tkbutton(mytt, text = ">", width = 2, command = function() moverightfunc())
    tkplace(z3.6, relx = 0.078, rely = 0.91, `in` = img3)
    z3.7 <- tkbutton(mytt, text = "C", width = 2, command = function() origpos())
    tkplace(z3.7, relx = 0.04, rely = 0.91, `in` = img3)
    ClearMainPoints4 <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab4") {
            MGvar$indexLabeled <<- c()
            MGvar$indexLabeled.T4 <<- c()
            MGvar$labeledPoints <<- list()
            tabplot()
        }
        else {
            tkmessageBox(message = "This is not the active plot", 
                icon = "error")
        }
    }
    MainPlotMenu4 <- tkmenu(img4, tearoff = FALSE)
    tkadd(MainPlotMenu4, "command", label = "                        Plot 4 Menu", 
        font = fontMenuHeading)
    tkadd(MainPlotMenu4, "separator")
    tkadd(MainPlotMenu4, "command", label = "Label Specific Point", 
        command = LabelSpecificPoint)
    tkadd(MainPlotMenu4, "command", label = "Clear Added Point Labels", 
        command = ClearMainPoints4)
    tkadd(MainPlotMenu4, "separator")
    tkadd(MainPlotMenu4, "command", label = "Relocate Group of Points", 
        command = MoveButtonInitiate)
    tkadd(MainPlotMenu4, "command", label = "Remove a Point", 
        command = RemovePoint)
    tkadd(MainPlotMenu4, "command", label = "Use Coordinates as Starting Configuration", 
        command = NewStartConfig)
    tkadd(MainPlotMenu4, "separator")
    tkadd(MainPlotMenu4, "command", label = "Rotate and Reflect", 
        command = RotateandReflect)
    tkadd(MainPlotMenu4, "command", label = "Advanced Zoom", 
        command = ZoomPlot)
    tkadd(MainPlotMenu4, "separator")
    tkadd(MainPlotMenu4, "command", label = "Change Point Colour", 
        command = ColButtonInitiate)
    tkadd(MainPlotMenu4, "command", label = "Default Point Colours", 
        command = DefaultPtCols)
    tkadd(MainPlotMenu4, "separator")
    MenuDispAx4 <- tclVar(0)
    tkadd(MainPlotMenu4, "checkbutton", label = "Display Variable Axes", 
        variable = MenuDispAx4, command = ShowregfromMenu)
    tkadd(MainPlotMenu4, "command", label = "Remove Axes of Variable(s)", 
        command = AxesRemove, state = "disabled")
    tkadd(MainPlotMenu4, "separator")
    tkadd(MainPlotMenu4, "command", label = "Pop-out Plot", command = EnlargedActivePlot)
    tkadd(MainPlotMenu4, "command", label = "Copy Plot4 to Clipboard", 
        command = function() {
            tkrreplot(img4, function() plotting2D(MGvar$MDSmat.T4, 
                title = MGvar$activeplot.title.T4, Measure = MGvar$dMeas.T4, 
                showtitle = MGvar$activeplot.title.show.T4, showmeas = MGvar$activeplot.distmeas.T4, 
                xlabel = MGvar$activeplot.xlab.T4, ylabel = MGvar$activeplot.ylab.T4, 
                bgcol = MGvar$activeplot.bg.T4, pointcex = MGvar$activeplot.cex.T4, 
                showlabs = MGvar$activeplot.labs.T4, showpoints = MGvar$activeplot.showpoints.T4, 
                pointcol = MGvar$activeplot.pointcol.T4, pointshape = MGvar$activeplot.type.T4, 
                ymeas = MGvar$activeplot.yaxt.T4, xmeas = MGvar$activeplot.xaxt.T4, 
                axcol = MGvar$activeplot.axescol.T4, indexLabeled = MGvar$indexLabeled.T4, 
                zoomedcoords = MGvar$newCoords.T4, Zrat = MGvar$zoominrat.T4, 
                Mvup = MGvar$moveup.T4, Mvdn = MGvar$movedown.T4, 
                Mvlt = MGvar$moveleft.T4, Mvrt = MGvar$moveright.T4))
        })
    tkadd(MainPlotMenu4, "separator")
    tkadd(MainPlotMenu4, "command", label = "Plot Options", command = ConfPlotOptions)
    RightClickMain4 <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", img4))
        rooty <- as.integer(tkwinfo("rooty", img4))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", MainPlotMenu4, xTxt, yTxt)
    }
    tkbind(img4, "<Button-3>", RightClickMain4)
    z4.1 <- tkbutton(mytt, text = "+", width = 2, command = function() simplezoomin())
    tkplace(z4.1, relx = 0.001, rely = 0.105, `in` = img4)
    z4.2 <- tkbutton(mytt, text = "-", width = 2, command = function() simplezoomout())
    tkplace(z4.2, relx = 0.001, rely = 0.15, `in` = img4)
    z4.3 <- tkbutton(mytt, text = "^", width = 2, command = function() moveupfunc())
    tkplace(z4.3, relx = 0.04, rely = 0.865, `in` = img4)
    z4.4 <- tkbutton(mytt, text = "v", width = 2, command = function() movedownfunc())
    tkplace(z4.4, relx = 0.04, rely = 0.955, `in` = img4)
    z4.5 <- tkbutton(mytt, text = "<", width = 2, command = function() moveleftfunc())
    tkplace(z4.5, relx = 0.001, rely = 0.91, `in` = img4)
    z4.6 <- tkbutton(mytt, text = ">", width = 2, command = function() moverightfunc())
    tkplace(z4.6, relx = 0.078, rely = 0.91, `in` = img4)
    z4.7 <- tkbutton(mytt, text = "C", width = 2, command = function() origpos())
    tkplace(z4.7, relx = 0.04, rely = 0.91, `in` = img4)
    ClearMainPoints5 <- function() {
        if (tclvalue(MGvar$ActivePlottingTab) == "Tab5") {
            MGvar$indexLabeled <<- c()
            MGvar$indexLabeled.T5 <<- c()
            MGvar$labeledPoints <<- list()
            tabplot()
        }
        else {
            tkmessageBox(message = "This is not the active plot", 
                icon = "error")
        }
    }
    MainPlotMenu5 <- tkmenu(img5, tearoff = FALSE)
    tkadd(MainPlotMenu5, "command", label = "                        Plot 5 Menu", 
        font = fontMenuHeading)
    tkadd(MainPlotMenu5, "separator")
    tkadd(MainPlotMenu5, "command", label = "Label Specific Point", 
        command = LabelSpecificPoint)
    tkadd(MainPlotMenu5, "command", label = "Clear Added Point Labels", 
        command = ClearMainPoints5)
    tkadd(MainPlotMenu5, "separator")
    tkadd(MainPlotMenu5, "command", label = "Relocate Group of Points", 
        command = MoveButtonInitiate)
    tkadd(MainPlotMenu5, "command", label = "Remove a Point", 
        command = RemovePoint)
    tkadd(MainPlotMenu5, "command", label = "Use Coordinates as Starting Configuration", 
        command = NewStartConfig)
    tkadd(MainPlotMenu5, "separator")
    tkadd(MainPlotMenu5, "command", label = "Rotate and Reflect", 
        command = RotateandReflect)
    tkadd(MainPlotMenu5, "command", label = "Advanced Zoom", 
        command = ZoomPlot)
    tkadd(MainPlotMenu5, "separator")
    tkadd(MainPlotMenu5, "command", label = "Change Point Colour", 
        command = ColButtonInitiate)
    tkadd(MainPlotMenu5, "command", label = "Default Point Colours", 
        command = DefaultPtCols)
    tkadd(MainPlotMenu5, "separator")
    MenuDispAx5 <- tclVar(0)
    tkadd(MainPlotMenu5, "checkbutton", label = "Display Variable Axes", 
        variable = MenuDispAx5, command = ShowregfromMenu)
    tkadd(MainPlotMenu5, "command", label = "Remove Axes of Variable(s)", 
        command = AxesRemove, state = "disabled")
    tkadd(MainPlotMenu5, "separator")
    tkadd(MainPlotMenu5, "command", label = "Pop-out Plot", command = EnlargedActivePlot)
    tkadd(MainPlotMenu5, "command", label = "Copy Plot5 to Clipboard", 
        command = function() {
            tkrreplot(img5, function() plotting2D(MGvar$MDSmat.T5, 
                title = MGvar$activeplot.title.T5, Measure = MGvar$dMeas.T5, 
                showtitle = MGvar$activeplot.title.show.T5, showmeas = MGvar$activeplot.distmeas.T5, 
                xlabel = MGvar$activeplot.xlab.T5, ylabel = MGvar$activeplot.ylab.T5, 
                bgcol = MGvar$activeplot.bg.T5, pointcex = MGvar$activeplot.cex.T5, 
                showlabs = MGvar$activeplot.labs.T5, showpoints = MGvar$activeplot.showpoints.T5, 
                pointcol = MGvar$activeplot.pointcol.T5, pointshape = MGvar$activeplot.type.T5, 
                ymeas = MGvar$activeplot.yaxt.T5, xmeas = MGvar$activeplot.xaxt.T5, 
                axcol = MGvar$activeplot.axescol.T5, indexLabeled = MGvar$indexLabeled.T5, 
                zoomedcoords = MGvar$newCoords.T5, Zrat = MGvar$zoominrat.T5, 
                Mvup = MGvar$moveup.T5, Mvdn = MGvar$movedown.T5, 
                Mvlt = MGvar$moveleft.T5, Mvrt = MGvar$moveright.T5))
        })
    tkadd(MainPlotMenu5, "separator")
    tkadd(MainPlotMenu5, "command", label = "Plot Options", command = ConfPlotOptions)
    RightClickMain5 <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", img5))
        rooty <- as.integer(tkwinfo("rooty", img5))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", MainPlotMenu5, xTxt, yTxt)
    }
    tkbind(img5, "<Button-3>", RightClickMain5)
    z5.1 <- tkbutton(mytt, text = "+", width = 2, command = function() simplezoomin())
    tkplace(z5.1, relx = 0.001, rely = 0.105, `in` = img5)
    z5.2 <- tkbutton(mytt, text = "-", width = 2, command = function() simplezoomout())
    tkplace(z5.2, relx = 0.001, rely = 0.15, `in` = img5)
    z5.3 <- tkbutton(mytt, text = "^", width = 2, command = function() moveupfunc())
    tkplace(z5.3, relx = 0.04, rely = 0.865, `in` = img5)
    z5.4 <- tkbutton(mytt, text = "v", width = 2, command = function() movedownfunc())
    tkplace(z5.4, relx = 0.04, rely = 0.955, `in` = img5)
    z5.5 <- tkbutton(mytt, text = "<", width = 2, command = function() moveleftfunc())
    tkplace(z5.5, relx = 0.001, rely = 0.91, `in` = img5)
    z5.6 <- tkbutton(mytt, text = ">", width = 2, command = function() moverightfunc())
    tkplace(z5.6, relx = 0.078, rely = 0.91, `in` = img5)
    z5.7 <- tkbutton(mytt, text = "C", width = 2, command = function() origpos())
    tkplace(z5.7, relx = 0.04, rely = 0.91, `in` = img5)
    Stat3DPlotMenu <- tkmenu(img3Dstat, tearoff = FALSE)
    tkadd(Stat3DPlotMenu, "command", label = "Pop-Out Static 3D Plot", 
        command = EnlargedStat3DPlot)
    tkadd(Stat3DPlotMenu, "command", label = "Copy Static 3D Plot to Clipboard", 
        command = function() {
            tkrreplot(img3Dstat)
        })
    tkadd(Stat3DPlotMenu, "separator")
    tkadd(Stat3DPlotMenu, "command", label = "Static 3D Plot Options", 
        command = threeDPlotOptions)
    RightClickStat3D <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", img3Dstat))
        rooty <- as.integer(tkwinfo("rooty", img3Dstat))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", Stat3DPlotMenu, xTxt, yTxt)
    }
    tkbind(img3Dstat, "<Button-3>", RightClickStat3D)
    z6.1 <- tkbutton(mytt, text = "<", width = 2, command = function() angle3Ddown())
    tkbind(mytt, "a", angle3Ddown)
    tkplace(z6.1, relx = 0.881, rely = 0.91, `in` = img3Dstat)
    z6.2 <- tkbutton(mytt, text = ">", width = 2, command = function() angle3Dup())
    tkbind(mytt, "d", angle3Dup)
    tkplace(z6.2, relx = 0.96, rely = 0.91, `in` = img3Dstat)
    z6.3 <- tkbutton(mytt, text = "^", width = 2, command = function() scale3Dup())
    tkbind(mytt, "w", scale3Dup)
    tkplace(z6.3, relx = 0.92, rely = 0.865, `in` = img3Dstat)
    z6.4 <- tkbutton(mytt, text = "v", width = 2, command = function() scale3Ddown())
    tkbind(mytt, "s", scale3Ddown)
    tkplace(z6.4, relx = 0.92, rely = 0.955, `in` = img3Dstat)
    z6.5 <- tkbutton(mytt, text = "C", width = 2, command = function() orig3D())
    tkbind(mytt, "C", orig3D)
    tkplace(z6.5, relx = 0.92, rely = 0.91, `in` = img3Dstat)
    Myhscale2 <- 1.2
    Myvscale2 <- 0.95
    mySecondaryNB <- tk2notebook(mytt, tabs = NULL)
    tb1.1 <- tk2frame(mySecondaryNB)
    tkadd(mySecondaryNB, tb1.1, text = "Shepard")
    frameSecTab1 <- tkframe(tb1.1, relief = "groove", borderwidth = 2, 
        background = "black")
    imgshep <- tkrplot(frameSecTab1, function() plotShepard(MGvar$distmat, 
        MGvar$MDSmat), hscale = Myhscale2, vscale = Myvscale2)
    tkplace(imgshep, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = frameSecTab1)
    tkplace(frameSecTab1, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb1.1)
    tkbind(imgshep, "<Button-1>", Shep.OnLeftClick)
    tkconfigure(imgshep, cursor = "hand2")
    Shep.WhichShow.var <- tclVar("Active Plot")
    Shep.WhichShow.ComboBox <- tkwidget(tb1.1, "ComboBox", editable = FALSE, 
        values = c("Active Plot", "Plot1", "Plot2", "Plot3", 
            "Plot4", "Plot5", "Stat3D", "RGL3D"), width = 10, 
        textvariable = Shep.WhichShow.var)
    tkplace(Shep.WhichShow.ComboBox, relx = 0.005, rely = 0.005, 
        `in` = frameSecTab1)
    AShepdistmat <- MGvar$distmat
    AShepdistmat.T1 <- MGvar$distmat
    AShepdistmat.T2 <- MGvar$distmat
    AShepdistmat.T3 <- MGvar$distmat
    AShepdistmat.T4 <- MGvar$distmat
    AShepdistmat.T5 <- MGvar$distmat
    AShepMDSmat <- MGvar$MDSmat
    AShepMDSmat.T1 <- MGvar$MDSmat
    AShepMDSmat.T2 <- MGvar$MDSmat
    AShepMDSmat.T3 <- MGvar$MDSmat
    AShepMDSmat.T4 <- MGvar$MDSmat
    AShepMDSmat.T5 <- MGvar$MDSmat
    OnShep.Go <- function() {
        if (tclvalue(Shep.WhichShow.var) == "Active Plot") {
            tkrreplot(imgshep, function() plotShepard(AShepdistmat, 
                AShepMDSmat, Tab = "Tab1"))
        }
        if (tclvalue(Shep.WhichShow.var) == "Plot1") {
            tkrreplot(imgshep, function() plotShepard(AShepdistmat.T1, 
                AShepMDSmat.T1, Tab = "Tab1"))
        }
        if (tclvalue(Shep.WhichShow.var) == "Plot2") {
            tkrreplot(imgshep, function() plotShepard(AShepdistmat.T2, 
                AShepMDSmat.T2, Tab = "Tab2"))
        }
        if (tclvalue(Shep.WhichShow.var) == "Plot3") {
            tkrreplot(imgshep, function() plotShepard(AShepdistmat.T3, 
                AShepMDSmat.T3, Tab = "Tab3"))
        }
        if (tclvalue(Shep.WhichShow.var) == "Plot4") {
            tkrreplot(imgshep, function() plotShepard(AShepdistmat.T4, 
                AShepMDSmat.T4, Tab = "Tab4"))
        }
        if (tclvalue(Shep.WhichShow.var) == "Plot5") {
            tkrreplot(imgshep, function() plotShepard(AShepdistmat.T5, 
                AShepMDSmat.T5, Tab = "Tab5"))
        }
    }
    tkplace(tkbutton(tb1.1, text = "Go", command = function() OnShep.Go()), 
        relx = 0.19, rely = 0.005, relheight = 0.06, `in` = frameSecTab1)
    tb1.5 <- tkframe(mySecondaryNB)
    tkadd(mySecondaryNB, tb1.5, text = "Stress Plot")
    frameSecTab5 <- tkframe(tb1.5, relief = "groove", borderwidth = 2, 
        background = "black")
    imgstress <- tkrplot(frameSecTab5, function() StressPlot(), 
        hscale = Myhscale2, vscale = Myvscale2)
    tkplace(imgstress, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = frameSecTab5)
    tkplace(frameSecTab5, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb1.5)
    tb1.6 <- tkframe(mySecondaryNB)
    tkadd(mySecondaryNB, tb1.6, text = "Log Stress")
    frameSecTab6 <- tkframe(tb1.6, relief = "groove", borderwidth = 2, 
        background = "black")
    imgstress2 <- tkrplot(frameSecTab6, function() StressPlot2(), 
        hscale = Myhscale2, vscale = Myvscale2)
    tkplace(imgstress2, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = frameSecTab6)
    tkplace(frameSecTab6, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb1.6)
    tb1.2 <- tk2frame(mySecondaryNB)
    tkadd(mySecondaryNB, tb1.2, text = "Scree")
    frameSecTab2 <- tkframe(tb1.2, relief = "groove", borderwidth = 2, 
        background = "black")
    imgscree <- tkrplot(frameSecTab2, function() plotScree(MGvar$scree.stress, 
        MGvar$screepoints.current, MGvar$screepoints.best, MGvar$MDS.dimensions, 
        MGvar$Opt.dim), hscale = Myhscale2, vscale = Myvscale2)
    tkplace(imgscree, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = frameSecTab2)
    tkplace(frameSecTab2, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb1.2)
    Scree.WhichShow.var <- tclVar("Active Plot")
    Scree.WhichShow.ComboBox <- tkwidget(tb1.2, "ComboBox", editable = FALSE, 
        values = c("Active Plot", "Plot1", "Plot2", "Plot3", 
            "Plot4", "Plot5", "Stat3D", "RGL3D"), width = 10, 
        textvariable = Scree.WhichShow.var)
    tkplace(Scree.WhichShow.ComboBox, relx = 0.005, rely = 0.005, 
        `in` = frameSecTab2)
    OnScree.Go <- function() {
        if (tclvalue(Scree.WhichShow.var) == "Active Plot") {
            tkrreplot(imgscree)
        }
        if (tclvalue(Scree.WhichShow.var) == "Plot1") {
            tkrreplot(imgscree, function() plotScree(MGvar$scree.stress.T1, 
                MGvar$screepoints.current.T1, MGvar$screepoints.best.T1, 
                MGvar$MDS.dimensions, MGvar$Opt.dim.T1, Tab = "Tab1"))
        }
        if (tclvalue(Scree.WhichShow.var) == "Plot2") {
            tkrreplot(imgscree, function() plotScree(MGvar$scree.stress.T2, 
                MGvar$screepoints.current.T2, MGvar$screepoints.best.T2, 
                MGvar$MDS.dimensions, MGvar$Opt.dim.T2, Tab = "Tab2"))
        }
        if (tclvalue(Scree.WhichShow.var) == "Plot3") {
            tkrreplot(imgscree, function() plotScree(MGvar$scree.stress.T3, 
                MGvar$screepoints.current.T3, MGvar$screepoints.best.T3, 
                MGvar$MDS.dimensions, MGvar$Opt.dim.T3, Tab = "Tab3"))
        }
        if (tclvalue(Scree.WhichShow.var) == "Plot4") {
            tkrreplot(imgscree, function() plotScree(MGvar$scree.stress.T4, 
                MGvar$screepoints.current.T4, MGvar$screepoints.best.T4, 
                MGvar$MDS.dimensions, MGvar$Opt.dim.T4, Tab = "Tab4"))
        }
        if (tclvalue(Scree.WhichShow.var) == "Plot5") {
            tkrreplot(imgscree, function() plotScree(MGvar$scree.stress.T5, 
                MGvar$screepoints.current.T5, MGvar$screepoints.best.T5, 
                MGvar$MDS.dimensions, MGvar$Opt.dim.T5, Tab = "Tab5"))
        }
    }
    tkplace(tkbutton(tb1.2, text = "Go", command = function() OnScree.Go()), 
        relx = 0.19, rely = 0.005, relheight = 0.06, `in` = frameSecTab2)
    tb1.3 <- tk2frame(mySecondaryNB)
    tkadd(mySecondaryNB, tb1.3, text = "Zoomed Plot")
    frameSecTab3 <- tkframe(tb1.3, relief = "groove", borderwidth = 2, 
        background = "black")
    tkplace(frameSecTab3, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb1.3)
    tb1.4 <- tk2frame(mySecondaryNB)
    tkadd(mySecondaryNB, tb1.4, text = "Notes/Script")
    frameSecTab4 <- tkframe(tb1.4, relief = "groove", borderwidth = 2, 
        background = "black")
    frameSecTabtxt <- tkframe(frameSecTab4, relief = "groove", 
        borderwidth = 2, background = "white")
    Notetxtscr <- tkscrollbar(frameSecTabtxt, repeatinterval = 5, 
        command = function(...) tkyview(Notetxt, ...))
    Notetxt <- tk2text(frameSecTabtxt, width = 65, height = 20, 
        yscrollcommand = function(...) tkset(Notetxtscr, ...))
    tkgrid(Notetxt, Notetxtscr)
    tkgrid.configure(Notetxtscr, sticky = "ns")
    wfile <- ""
    tkplace(tkbutton(frameSecTab4, text = "Run as Code", width = 15, 
        command = function() runtextascode()), relx = 0.04, rely = 0.91, 
        `in` = frameSecTab4)
    tkplace(tkbutton(frameSecTab4, text = "Save Notes", width = 15, 
        command = function() savetext()), relx = 0.37, rely = 0.91, 
        `in` = frameSecTab4)
    tkplace(tkbutton(frameSecTab4, text = "Load Notes", width = 15, 
        command = function() loadtext()), relx = 0.7, rely = 0.91, 
        `in` = frameSecTab4)
    tkplace(frameSecTab4, relx = 0.01, relwidth = 0.98, rely = 0.01, 
        relheight = 0.98, `in` = tb1.4)
    tkplace(frameSecTabtxt, relx = 0.01, relwidth = 0.98, rely = 0.01, 
        relheight = 0.88, `in` = frameSecTab4)
    copyText <- function() tcl("event", "generate", .Tk.ID(Notetxt), 
        "<<Copy>>")
    PasteText <- function() tcl("event", "generate", .Tk.ID(Notetxt), 
        "<<Paste>>")
    SelectAllText <- function() tcl("event", "generate", .Tk.ID(Notetxt), 
        "<<Select All>>")
    TexteditPopupMenu <- tkmenu(Notetxt, tearoff = FALSE)
    tkadd(TexteditPopupMenu, "command", label = "Paste <Ctrl-V>", 
        command = PasteText)
    tkadd(TexteditPopupMenu, "command", label = "Copy <Ctrl-C>", 
        command = copyText)
    RightClickNoteText <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", Notetxt))
        rooty <- as.integer(tkwinfo("rooty", Notetxt))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", TexteditPopupMenu, xTxt, yTxt)
    }
    tkbind(Notetxt, "<Button-3>", RightClickNoteText)
    tkbind(Notetxt, SelectAllText)
    myTableNB <- tk2notebook(mytt, tabs = NULL)
    tb3.1 <- tk2frame(myTableNB)
    tkadd(myTableNB, tb3.1, text = "MDS Configurations")
    frameConfTab <- tkframe(tb3.1, relief = "groove", borderwidth = 2, 
        background = "black")
    tkplace(frameConfTab, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb3.1)
    tableplace()
    tb3.2 <- tk2frame(myTableNB)
    tkadd(myTableNB, tb3.2, text = "Removed Points")
    frameRemPoints <- tkframe(tb3.2, relief = "groove", borderwidth = 2, 
        background = "black")
    tkplace(frameRemPoints, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb3.2)
    tableplaceRemP()
    RemPointsPopupMenu <- tkmenu(MGcomp$table.Remp, tearoff = FALSE)
    tkadd(RemPointsPopupMenu, "command", label = "Replace Point in Active Cell", 
        command = PointReplace)
    tkadd(RemPointsPopupMenu, "separator")
    tkadd(RemPointsPopupMenu, "command", label = "Pop-Out Table")
    RightClickRemPoints <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", MGcomp$table.Remp))
        rooty <- as.integer(tkwinfo("rooty", MGcomp$table.Remp))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", RemPointsPopupMenu, xTxt, yTxt)
    }
    tkbind(MGcomp$table.Remp, "<Button-3>", RightClickRemPoints)
    tb3.3 <- tk2frame(myTableNB)
    tkadd(myTableNB, tb3.3, text = "Removed Axes")
    frameRemAx <- tkframe(tb3.3, relief = "groove", borderwidth = 2, 
        background = "black")
    tkplace(frameRemAx, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = tb3.3)
    tableplaceRemAx()
    RemAxesPopupMenu <- tkmenu(MGcomp$table.RemAx, tearoff = FALSE)
    tkadd(RemAxesPopupMenu, "command", label = "Replace Point in Active Cell", 
        command = AxesReplace)
    tkadd(RemAxesPopupMenu, "separator")
    tkadd(RemAxesPopupMenu, "command", label = "Pop-Out Table")
    RightClickRemAxes <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", MGcomp$table.RemAx))
        rooty <- as.integer(tkwinfo("rooty", MGcomp$table.RemAx))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", RemAxesPopupMenu, xTxt, yTxt)
    }
    tkbind(MGcomp$table.RemAx, "<Button-3>", RightClickRemAxes)
    tkplace(myPlottingNB, relx = 0.43, rely = 0.08, relwidth = 0.56, 
        relheight = 0.86, `in` = mytt)
    tkbind(myPlottingNB, "<ButtonRelease-1>", plottingTabLeftClick)
    tkplace(myTableNB, relx = 0.01, relwidth = 0.41, rely = 0.05, 
        relheight = 0.285, `in` = mytt)
    tkplace(mySecondaryNB, relx = 0.01, rely = 0.35, relwidth = 0.41, 
        relheight = 0.59, `in` = mytt)
    GoToButton <- tkbutton(myPlottingNB, text = "Go to Active Plot", 
        width = 15, command = function() GOtoActiveTab())
    tkplace(GoToButton, relx = 0.81, rely = 0, relheight = 0.038, 
        `in` = myPlottingNB)
    bottomfont <- tkfont.create(family = "times", size = 11, 
        weight = "bold")
    frameBottomStrip = tkframe(mytt, relief = "groove", borderwidth = 5, 
        background = "white")
    labelText <- tclVar("No Active DataSet")
    PlotText <- tclVar("Active Plot is Plot1")
    dimText <- tclVar("p = 2")
    tkplace(tklabel(frameBottomStrip, text = tclvalue(PlotText), 
        textvariable = PlotText, background = "white", font = bottomfont), 
        relx = 0.5, y = 0.01, `in` = frameBottomStrip)
    tkplace(tklabel(frameBottomStrip, text = tclvalue(labelText), 
        textvariable = labelText, background = "white", font = bottomfont), 
        relx = 0.33, y = 0.01, `in` = frameBottomStrip)
    tkplace(tklabel(frameBottomStrip, text = tclvalue(dimText), 
        textvariable = dimText, background = "white", font = bottomfont), 
        relx = 0.7, y = 0.01, `in` = frameBottomStrip)
    developerfont <- tkfont.create(family = "times", size = 6)
    tkplace(tklabel(frameBottomStrip, text = "MDS-GUI Version 1.0. built in 2012\nAndrew Timm & Sugnet Lubbe, University of Cape Town", 
        background = "white", font = developerfont), relx = 0.8, 
        y = 0.001, `in` = frameBottomStrip)
    tkplace(frameBottomStrip, relx = 0.01, rely = 0.95, relwidth = 0.98, 
        relheight = 0.051, `in` = mytt)
    EndProcess <- function() {
        MGvar$endprocess <<- "yes"
        ActiveArrowCursor()
        ActivateAll()
    }
    EndButton <- tkbutton(mytt, text = "End Process", command = function() EndProcess())
    tkplace(EndButton, relx = 0.001, rely = 0.01, `in` = frameBottomStrip)
    MDSGUI.ProgressBar <- ttkprogressbar(mytt, mode = "determinate")
    MDSGUI.ProgressBar.create <- function() {
        tkplace(MDSGUI.ProgressBar, relx = 0.07, rely = 0.01, 
            relwidth = 0.15, relheight = 0.98, `in` = frameBottomStrip)
    }
    MDSGUI.ProgressBar.create()
    iterrat <- tclVar(paste(MGvar$MDS.iter.max, "iters"))
    iterdispfont <- tkfont.create(family = "times", size = 6, 
        weight = "bold")
    tkplace(tklabel(frameBottomStrip, text = tclvalue(iterrat), 
        textvariable = iterrat, background = "white", font = iterdispfont), 
        relx = 0.22, rely = 0.3, `in` = frameBottomStrip)
    Undolast <- function() tcl("event", "generate", .Tk.ID(mytt), 
        "<<Undo>>")
    ShepPopupMenu <- tkmenu(imgshep, tearoff = FALSE)
    tkadd(ShepPopupMenu, "command", label = "         Shepard Plot Menu", 
        font = fontMenuHeading)
    tkadd(ShepPopupMenu, "separator")
    tkadd(ShepPopupMenu, "command", label = "Clear Added Labels", 
        command = ClearShepLabels)
    tkadd(ShepPopupMenu, "command", label = "Label Specific Point", 
        command = Shep.LabelSpecificPoint)
    tkadd(ShepPopupMenu, "separator")
    tkadd(ShepPopupMenu, "command", label = "Pop-Out Enlarged Plot", 
        command = EnlargedShep)
    tkadd(ShepPopupMenu, "command", label = "Copy Shepard Plot to Clipboard", 
        command = function() {
            tkrreplot(imgshep)
        })
    tkadd(ShepPopupMenu, "separator")
    tkadd(ShepPopupMenu, "command", label = "Shepard Plot Options", 
        command = ShepPlotOptions)
    RightClickShep <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", imgshep))
        rooty <- as.integer(tkwinfo("rooty", imgshep))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", ShepPopupMenu, xTxt, yTxt)
    }
    tkbind(imgshep, "<Button-3>", RightClickShep)
    tkbind(imgshep, "<B1-Motion>", BrushingShep)
    tkbind(imgshep, "<ButtonRelease-1>", OnShepRelease)
    EnScree.switch <- "off"
    ScreePopupMenu <- tkmenu(imgscree, tearoff = FALSE)
    tkadd(ScreePopupMenu, "command", label = "           Scree Plot Menu", 
        font = fontMenuHeading)
    tkadd(ScreePopupMenu, "separator")
    tkadd(ScreePopupMenu, "command", label = "Pop-Out Enlarged Plot", 
        command = EnlargedScree)
    tkadd(ScreePopupMenu, "command", label = "Copy Scree/Eigen Plot to Clipboard", 
        command = function() {
            tkrreplot(imgscree)
        })
    tkadd(ScreePopupMenu, "separator")
    tkadd(ScreePopupMenu, "command", label = "Scree Plot Options", 
        command = ScreePlotOptions)
    RightClickScree <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", imgscree))
        rooty <- as.integer(tkwinfo("rooty", imgscree))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", ScreePopupMenu, xTxt, yTxt)
    }
    tkbind(imgscree, "<Button-3>", RightClickScree)
    MGvar$EnStress.switch <- "off"
    StressPopupMenu <- tkmenu(imgstress, tearoff = FALSE)
    tkadd(StressPopupMenu, "command", label = "           Stress Plot Menu", 
        font = fontMenuHeading)
    tkadd(StressPopupMenu, "separator")
    tkadd(StressPopupMenu, "command", label = "Pop-Out Enlarged Plot", 
        command = EnlargedStress)
    tkadd(StressPopupMenu, "command", label = "Copy Stress Plot to Clipboard", 
        command = function() {
            tkrreplot(imgstress)
        })
    tkadd(StressPopupMenu, "separator")
    tkadd(StressPopupMenu, "command", label = "Stress Plot Options", 
        command = StressPlotOps)
    RightClickStress <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", imgstress))
        rooty <- as.integer(tkwinfo("rooty", imgstress))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", StressPopupMenu, xTxt, yTxt)
    }
    tkbind(imgstress, "<Button-3>", RightClickStress)
    MGvar$EnStress2.switch <- "off"
    Stress2PopupMenu <- tkmenu(imgstress2, tearoff = FALSE)
    tkadd(Stress2PopupMenu, "command", label = "           Stress Plot Menu", 
        font = fontMenuHeading)
    tkadd(Stress2PopupMenu, "separator")
    tkadd(Stress2PopupMenu, "command", label = "Pop-Out Enlarged Plot", 
        command = EnlargedStress2)
    tkadd(Stress2PopupMenu, "command", label = "Copy Stress Plot to Clipboard", 
        command = function() {
            tkrreplot(imgstress2)
        })
    tkadd(Stress2PopupMenu, "separator")
    tkadd(Stress2PopupMenu, "command", label = "Stress Plot Options", 
        command = StressPlotOps)
    RightClickStress2 <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", imgstress2))
        rooty <- as.integer(tkwinfo("rooty", imgstress2))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", Stress2PopupMenu, xTxt, yTxt)
    }
    tkbind(imgstress2, "<Button-3>", RightClickStress2)
    imgseczoom <- tkrplot(frameSecTab3, function() seczoomplot(MGvar$newCoords), 
        hscale = 1.1, vscale = 0.9)
    tkplace(imgseczoom, relx = 0.005, relwidth = 0.99, rely = 0.005, 
        relheight = 0.99, `in` = frameSecTab3)
    secZoomplotMenu <- tkmenu(imgseczoom, tearoff = FALSE)
    tkadd(secZoomplotMenu, "command", label = "Zoom Plot Options", 
        command = ZoomPlotOps)
    tkadd(secZoomplotMenu, "command", label = "Copy To Clipboard", 
        command = function() {
            tkrreplot(imgseczoom)
        })
    RightClickSecZoom <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", imgseczoom))
        rooty <- as.integer(tkwinfo("rooty", imgseczoom))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", secZoomplotMenu, xTxt, yTxt)
    }
    tkbind(imgseczoom, "<Button-3>", RightClickSecZoom)
    ClearZoomPoints <- function() {
        MGvar$indexZLabeled <<- c()
        tkrreplot(MGcomp$zooming)
    }
    RightClickZoom <- function(x, y) {
        rootx <- as.integer(tkwinfo("rootx", MGcomp$zooming))
        rooty <- as.integer(tkwinfo("rooty", MGcomp$zooming))
        xTxt <- as.integer(x) + rootx
        yTxt <- as.integer(y) + rooty
        tcl("tk_popup", MGcomp$ZoomPlotMenu, xTxt, yTxt)
    }
    MGvar$seczoomswitch <- "off"
    popzoomswitch <- "off"
    activeX <- 0
    fileMenu <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(fileMenu, "command", label = "New User", command = newUser)
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Save User WorkSpace", 
        state = "disabled")
    tkadd(fileMenu, "command", label = "Load User WorkSpace", 
        state = "disabled")
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Print", command = MDSGUI.print)
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Clear All", command = Clear)
    tkadd(fileMenu, "separator")
    tkadd(fileMenu, "command", label = "Exit MDS-GUI", command = QuitGUI)
    tkadd(topMenu, "cascade", label = "File", menu = fileMenu)
    GeneralMenu <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(GeneralMenu, "command", label = "Undo", command = UndoMove, 
        state = "disabled")
    tkadd(GeneralMenu, "command", label = "Appearance Settings", 
        command = Appearance.Change)
    tkadd(GeneralMenu, "command", label = "General Settings", 
        command = GeneralSettings)
    PlotExportMenu <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(PlotExportMenu, "command", label = "All", command = function() {
        exportsweave(All = TRUE)
    })
    tkadd(PlotExportMenu, "command", label = "Plot1 info", command = function() {
        exportsweave(All = FALSE, P1 = TRUE)
    })
    tkadd(PlotExportMenu, "command", label = "Plot2 info", command = function() {
        exportsweave(All = FALSE, P2 = TRUE)
    })
    tkadd(PlotExportMenu, "command", label = "Plot3 info", command = function() {
        exportsweave(All = FALSE, P3 = TRUE)
    })
    tkadd(PlotExportMenu, "command", label = "Plot4 info", command = function() {
        exportsweave(All = FALSE, P4 = TRUE)
    })
    tkadd(PlotExportMenu, "command", label = "Plot5 info", command = function() {
        exportsweave(All = FALSE, P5 = TRUE)
    })
    tkadd(GeneralMenu, "cascade", label = "Export", menu = PlotExportMenu)
    tkadd(topMenu, "cascade", label = "General", menu = GeneralMenu)
    dataMenu <- tkmenu(topMenu, tearoff = FALSE)
    LoadMenu <- tkmenu(topMenu, tearoff = FALSE)
    SaveMenu <- tkmenu(topMenu, tearoff = FALSE)
    EditDataMenu <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(dataMenu, "cascade", label = "Load Dataset (txt)", command = LoadDataSettxt)
    tkadd(dataMenu, "cascade", label = "Load Dataset (csv)", command = LoadDataSetcsv)
    tkbind(mytt, "<Control-L>", LoadDataSettxt)
    tkadd(dataMenu, "command", label = "Load Dissimilarity Matrix", 
        command = LoadDistSet)
    tkbind(mytt, "<Control-D>", LoadDistSet)
    tkadd(dataMenu, "command", label = "Load Similarity Matrix", 
        command = LoadSimSet)
    tkbind(mytt, "<Control-S>", LoadSimSet)
    tkadd(dataMenu, "command", label = "Load Correlation Matrix", 
        command = LoadSimSet)
    tkadd(dataMenu, "separator")
    tkadd(dataMenu, "command", label = "Data Colour Index", command = ColoursDataOption)
    tkadd(dataMenu, "command", label = "Catagory Colours", command = CatCols, 
        state = "disabled")
    tkadd(dataMenu, "separator")
    tkadd(EditDataMenu, "command", label = "Edit Active Data", 
        command = ShowActiveData, state = "disabled")
    tkadd(EditDataMenu, "command", label = "Edit Active Dissimilarity Matrix", 
        command = ShowActiveDistmat, state = "disabled")
    tkadd(EditDataMenu, "command", label = "Edit Active Similarity Matrix", 
        command = ShowActiveSimmat, state = "disabled")
    tkadd(EditDataMenu, "command", label = "Edit Active MDS Correlation Matrix", 
        command = ShowActiveMDSmat, state = "disabled")
    tkadd(EditDataMenu, "command", label = "Edit Active Coordinate Vectors", 
        command = ShowActiveMDSmat, state = "disabled")
    tkadd(dataMenu, "cascade", label = "Edit", menu = EditDataMenu)
    tkadd(dataMenu, "separator")
    tkadd(SaveMenu, "command", label = "DataSet", command = SaveDataSet)
    tkadd(SaveMenu, "command", label = "Dissimilarity Matrix", 
        command = SaveDistmat)
    tkadd(SaveMenu, "command", label = "Similarity Matrix", command = SaveSimmat)
    tkadd(SaveMenu, "command", label = "Correlation Matrix", 
        command = SaveCormat)
    tkadd(SaveMenu, "command", label = "MDS Coordinate Matrix", 
        command = SaveMDScoords)
    tkadd(dataMenu, "cascade", label = "Save", menu = SaveMenu)
    tkadd(dataMenu, "separator")
    tkadd(dataMenu, "command", label = "Data Options", command = DataOptions)
    tkadd(dataMenu, "separator")
    tkadd(topMenu, "cascade", label = "Data", menu = dataMenu)
    functionMenu <- tkmenu(topMenu, tearoff = FALSE)
    MDSMenu <- tkmenu(topMenu, tearoff = FALSE)
    MDSMenuData <- tkmenu(topMenu, tearoff = FALSE)
    MDSMenuDistmat <- tkmenu(topMenu, tearoff = FALSE)
    MDSMenuDistCal <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(MDSMenu, "separator")
    tkadd(MDSMenu, "command", label = "                 Metric MDS", 
        font = fontMenuHeading)
    tkadd(MDSMenu, "separator")
    tkadd(MDSMenu, "command", label = "Classical Scaling", command = function() {
        freshstart()
        ClasScal()
    })
    tkadd(MDSMenu, "command", label = "Metric Symmetric SMACOF", 
        command = function() {
            freshstart()
            MyMetricSmacofSym()
        })
    tkadd(MDSMenu, "command", label = "Least Squares Scaling", 
        command = function() {
            freshstart()
            MetricLeastSquaresScaling()
        })
    tkadd(MDSMenu, "separator")
    tkadd(MDSMenu, "command", label = "            Non-Metric MDS", 
        font = fontMenuHeading)
    tkadd(MDSMenu, "separator")
    tkadd(MDSMenu, "command", label = "Non-Metric Symmetric SMACOF", 
        command = function() {
            freshstart()
            MyNonMetricSmacofSym()
        })
    tkadd(MDSMenu, "command", label = "Kruskals Non-Metric", 
        command = function() {
            freshstart()
            MyKruskalsMDS()
        })
 	tkadd(MDSMenu, "command", label = "Sammon Mapping", command = function() {
        freshstart()
        SammonMapping()
    })   

    tkadd(MDSMenu, "separator")
    tkadd(MDSMenu, "separator")
    tkadd(functionMenu, "cascade", label = "MDS", menu = MDSMenu)
    tkadd(functionMenu, "cascade", label = "MDS Options", command = function() MDSOptions())
    tkadd(functionMenu, "separator")
    tkadd(functionMenu, "cascade", label = "Procrusts Analysis", 
        command = function() PerformProcrustes())
    tkadd(functionMenu, "separator")
    tkadd(functionMenu, "cascade", label = "Dissimilarity Matrix Calculation", 
        menu = MDSMenuDistCal)
    tkadd(MDSMenuDistCal, "radiobutton", label = "Euclidean", 
        variable = MGvar$dMeasVar, value = "Euc", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Weigthed Euclidean", 
        variable = MGvar$dMeasVar, value = "WEuc", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Mahalanobis Distance", 
        variable = MGvar$dMeasVar, value = "Mah", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "City Block Metric", 
        variable = MGvar$dMeasVar, value = "CBM", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Minkowski Metric", 
        variable = MGvar$dMeasVar, value = "Mink", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Canberra Metric", 
        variable = MGvar$dMeasVar, value = "Can", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Divergence", 
        variable = MGvar$dMeasVar, value = "Div", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Bray Curtis", 
        variable = MGvar$dMeasVar, value = "BC", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Soergel", variable = MGvar$dMeasVar, 
        value = "Soe", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Bhattacharyya Distance", 
        variable = MGvar$dMeasVar, value = "Bhat", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Wave Hedges", 
        variable = MGvar$dMeasVar, value = "WH", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Angular Seperation", 
        variable = MGvar$dMeasVar, value = "Ang", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(MDSMenuDistCal, "radiobutton", label = "Correlation", 
        variable = MGvar$dMeasVar, value = "Corr", command = function() {
            new.dissim.meas(MGvar$activedata)
        })
    tkadd(topMenu, "cascade", label = "Multivariate Tools", menu = functionMenu)
    Aboutfunc <- function() {
        tkmessageBox(title = "MDS-GUI Information", message = "This is the MDS-GUI Version 0.2\n\nDeveloped by Andrew Timm and Sugnet Lubbe, \nDepartment of Statistical Sciences, University of Cape Town.\n\nContact email: timmand@gmail.com\n\nLicense Info: GPL (>= 3)", 
            icon = "info", type = "ok")
    }
    helpMenu <- tkmenu(topMenu, tearoff = FALSE)
    tkadd(helpMenu, "command", label = "Function Code", command = FuncCode.Display)
    tkadd(helpMenu, "separator")
    tkadd(helpMenu, "checkbutton", variable = MGvar$DispHelp.Var, 
        label = "Display Pop-Out Help", command = DispPopUpHelp)
    tkadd(helpMenu, "separator")
    tkadd(helpMenu, "command", label = "Vignette", command = Documents.Help)
    tkadd(helpMenu, "command", label = "User Manual", command = Documents.Help)
    tkadd(helpMenu, "command", label = "HomePage", command = HomePage.Help)
    tkadd(helpMenu, "separator")
    tkadd(helpMenu, "command", label = "About", command = Aboutfunc)
    tkadd(topMenu, "cascade", label = "Help", menu = helpMenu)
    tkbind(mytt, "+", simplezoomin)
    tkbind(mytt, "-", simplezoomout)
    tkbind(mytt, "<Left>", moveleftfunc)
    tkbind(mytt, "<Right>", moverightfunc)
    tkbind(mytt, "<Up>", moveupfunc)
    tkbind(mytt, "<Down>", movedownfunc)
    tkbind(mytt, "c", origpos)
    tkbind(mytt, "1", function() {
        ChangeTabSC(1)
    })
    tkbind(mytt, "2", function() {
        ChangeTabSC(2)
    })
    tkbind(mytt, "3", function() {
        ChangeTabSC(3)
    })
    tkbind(mytt, "4", function() {
        ChangeTabSC(4)
    })
    tkbind(mytt, "5", function() {
        ChangeTabSC(5)
    })
    tkbind(mytt, "<Control-z>", UndoMove)
    tkfocus(mytt)
    tkwm.deiconify(mytt)
}
