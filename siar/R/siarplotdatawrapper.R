siarplotdatawrapper <-
function(siardata, siarversion = 0, grp = NULL, panel = NULL,
        isos = c(1, 2),leg2 = NULL,legloc='topleft') {
        if (!is.null(panel) & is.null(grp)) {
            warning(cat("WARNING. grp set to ALL and panel set to a value.\n Overriding your panel selection and setting to panel=NULL.\n In order to plot all groups on seperate panels please call\n grp=1:siardata$numgroups and panel=1 or panel=c(r,c)\n to specify number of rows and columns"))
            panel <- NULL
        }
        if (all(isos == 0)) {
            isox <- 1
            isoy <- 2
        }
        else {
            isox <- isos[1]
            isoy <- isos[2]
        }
        a <- 1
        if (siardata$numgroups == 1) {
            a <- 0
        }
        if (!is.null(panel)) {
            if (prod(panel) < length(grp)) {
                panel <- c(ceiling(sqrt(length(grp))))
                panel <- c(max(panel, 1), max(ceiling(length(grp)/panel),
                  1))
            }
            split.screen(panel)
        }
        else {
            #newgraphwindow()
        }
        er <- (siardata$sources[, (2 * isox) + 1]^2 + siardata$corrections[,
            (2 * isox) + 1]^2)^0.5
        xmins <- min(c(siardata$sources[, 2 * isox] + siardata$corrections[,
            2 * isox] - 3 * er, siardata$targets[, isox + a]))
        xmaxs <- max(c(siardata$sources[, 2 * isox] + siardata$corrections[,
            2 * isox] + 3 * er, siardata$targets[, isox + a]))
        er <- (siardata$sources[, (2 * isoy) + 1]^2 + siardata$corrections[,
            (2 * isoy) + 1]^2)^0.5
        ymins <- min(c(siardata$sources[, 2 * isoy] + siardata$corrections[,
            2 * isoy] - 3 * er, siardata$targets[, isoy + a]))
        ymaxs <- max(c(siardata$sources[, 2 * isoy] + siardata$corrections[,
            2 * isoy] + 3 * er, siardata$targets[, isoy + a]))
        if (is.null(panel)) {
            plot(1, 1, type = "n", xlim = c(xmins, xmaxs), ylim = c(ymins,
                ymaxs), main = siardata$TITLE, xlab = colnames(siardata$targets)[isox +
                a], ylab = colnames(siardata$targets)[isoy +
                a])
        }
        for (k in 1:length(grp)) {
            if (!is.null(panel)) {
                screen(k)
                plot(1, 1, type = "n", xlim = c(xmins, xmaxs),
                  ylim = c(ymins, ymaxs), main = paste("Group",
                    grp[k]), xlab = colnames(siardata$targets)[isox +
                    a], ylab = colnames(siardata$targets)[isoy +
                    a])
            }
            if (!is.null(grp)) {
                siarplottarget(siardata, isox, isoy, a, grps = grp[k])
            }
            else {
                siarplottarget(siardata, isox, isoy, a, grps = grp)
            }
            for (i in 1:nrow(siardata$sources)) {
                dx <- siardata$sources[i, 2 * isox] + siardata$corrections[i,
                  2 * isox]
                dex <- 2 * (siardata$sources[i, (2 * isox) +
                  1]^2 + siardata$corrections[i, (2 * isox) +
                  1]^2)^0.5
                dy <- siardata$sources[i, 2 * isoy] + siardata$corrections[i,
                  2 * isoy]
                dey <- 2 * (siardata$sources[i, (2 * isoy) +
                  1]^2 + siardata$corrections[i, (2 * isoy) +
                  1]^2)^0.5
                siaraddcross(x = dx, ex = dex, y = dy, ey = dey,
                  upch = 15, clr = i)
            }
        }
        if (!is.null(panel)) {
            close.screen(all.screens = TRUE)
        }
        if (siarversion > 0) {
            mtext(paste("siar v", siarversion), side = 1, line = 4,
                adj = 1, cex = 0.6)
        }
        if (siardata$numgroups == 0) {
            grp <- 1
        }
        if (is.null(grp)) {
            grp <- 1
        }
        pchseq <- c(1:2, 4:20)
        if (leg2 == 1) {
            datalabs <- NULL
            if (siardata$numgroups == 1) {
                datalabs <- "data"
            }
            else {
                for (k in 1:length(grp)) {
                  datalabs <- c(datalabs, as.character(paste("Group",
                    grp[k])))
                }
            }
            legend(legloc, legend = c(as.character(siardata$sources[,
                1]), datalabs), lty = c(rep(1, nrow(siardata$sources)),
                rep(-1, length(grp))), pch = c(rep(15, nrow(siardata$sources)),
                pchseq[grp]), col = c(seq(1, nrow(siardata$sources)),
                rep("grey50", length(grp))), bty = "n")
        }
        if (leg2 == 2) {
            datalabs <- NULL
            if (siardata$numgroups == 1) {
                datalabs <- "data"
            }
            else {
                for (k in 1:length(grp)) {
                  datalabs <- c(datalabs, as.character(paste("Group",
                    grp[k])))
                }
            }
            #newgraphwindow()
            plot(0, 0, "n", xaxt = "n", yaxt = "n", bty = "n")
            legend(0, 0, legend = c(as.character(siardata$sources[,
                1]), datalabs), lty = c(rep(1, nrow(siardata$sources)),
                rep(-1, length(grp))), pch = c(rep(15, nrow(siardata$sources)),
                pchseq[grp]), col = c(seq(1, nrow(siardata$sources)),
                rep("grey50", length(grp))), bty = "n")
        }
    }
