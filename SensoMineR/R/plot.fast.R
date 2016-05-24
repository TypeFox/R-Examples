plot.fast <- function (x, choix = "ind", axes = c(1, 2), xlim = NULL, ylim = NULL,
    invisible = NULL, col.ind = "blue", col.var = "red", col.quali.sup = "darkred",
    col.ind.sup = "darkblue", col.quanti.sup = "black", label = "all",
    cex = 1, lab.grpe = TRUE, title = NULL, habillage = "none",
    palette = NULL, new.plot = TRUE, ...)
{
    res.mca <- x
    if (!inherits(res.mca, "fast"))
        stop("non convenient data")
    if (is.null(palette))
        palette(c("black", "red", "green3", "blue", "cyan", "magenta",
            "darkgray", "darkgoldenrod", "darkgreen", "violet",
            "turquoise", "orange", "lightpink", "lavender", "yellow",
            "lightgreen", "lightgrey", "lightblue", "darkkhaki",
            "darkmagenta", "darkolivegreen", "lightcyan", "darkorange",
            "darkorchid", "darkred", "darksalmon", "darkseagreen",
            "darkslateblue", "darkslategray", "darkslategrey",
            "darkturquoise", "darkviolet", "lightgray", "lightsalmon",
            "lightyellow", "maroon"))
    lab.ind <- lab.var <- lab.quali.sup <- lab.ind.sup <- FALSE
    if (length(label) == 1 && label == "all")
        lab.ind <- lab.var <- lab.quali.sup <- lab.ind.sup <- TRUE
    if ("ind" %in% label)
        lab.ind <- TRUE
    if ("var" %in% label)
        lab.var <- TRUE
    test.invisible <- vector(length = 5)
    if (!is.null(invisible)) {
        test.invisible[1] <- match("ind", invisible)
        test.invisible[2] <- match("var", invisible)
    }
    else test.invisible <- rep(NA, 5)
    coord.var <- res.mca$var$coord[, axes]
    coord.ind <- res.mca$ind$coord[, axes]
    if (is.null(xlim)) {
        xmin <- xmax <- 0
        if (is.na(test.invisible[1]))
            xmin <- min(xmin, coord.ind[, 1])
        if (is.na(test.invisible[1]))
            xmax <- max(xmax, coord.ind[, 1])
        if (is.na(test.invisible[2]))
            xmin <- min(xmin, coord.var[, 1])
        if (is.na(test.invisible[2]))
            xmax <- max(xmax, coord.var[, 1])
        xlim <- c(xmin, xmax) * 1.2
    }
    else {
        xmin = xlim[1]
        xmax = xlim[2]
    }
    if (is.null(ylim)) {
        ymin <- ymax <- 0
        if (is.na(test.invisible[1]))
            ymin <- min(ymin, coord.ind[, 2])
        if (is.na(test.invisible[1]))
            ymax <- max(ymax, coord.ind[, 2])
        if (is.na(test.invisible[2]))
            ymin <- min(ymin, coord.var[, 2])
        if (is.na(test.invisible[2]))
            ymax <- max(ymax, coord.var[, 2])
        ylim <- c(ymin, ymax) * 1.2
    }
    else {
        ymin = ylim[1]
        ymax = ylim[2]
    }
    if (habillage == "quali") {
        aux = 1
        col.var = NULL
        for (j in res.mca$call$quali) {
            col.var <- c(col.var, rep(aux, nlevels(res.mca$call$X[,
                j])))
            aux = aux + 1
        }
    }
    if ((habillage != "none") & (habillage != "quali")) {
        if (!is.factor(res.mca$call$X[, habillage]))
            stop("The variable ", habillage, " is not qualitative")
        col.ind <- as.numeric(as.factor(res.mca$call$X[, habillage]))
        n.mod <- nlevels(as.factor(res.mca$call$X[, habillage]))
    }
    if (choix == "ind" | choix == "var") {
        sub.titre <- NULL
        titre = title
        if (is.null(title))
            titre <- "SortingTask factor map"
        else sub.titre <- "SortingTask factor map"
        if (is.na(test.invisible[1]) | is.na(test.invisible[2]) |
            is.na(test.invisible[4]) | is.na(test.invisible[5])) {
            if (new.plot)
                dev.new()
            plot(0, 0, main = titre, xlab = paste("Dim ", axes[1],
                " (", signif(res.mca$eig[axes[1], 2], 4), "%)",
                sep = ""), ylab = paste("Dim ", axes[2], " (",
                signif(res.mca$eig[axes[2], 2], 4), "%)", sep = ""),
                xlim = xlim, ylim = ylim, col = "white", asp = 1,
                cex = cex)
            if (!is.null(sub.titre))
                title(sub = sub.titre, cex.sub = cex, font.sub = 2,
                  col.sub = "steelblue4", adj = 0, line = 3.8)
            abline(v = 0, lty = 2, cex = cex)
            abline(h = 0, lty = 2, cex = cex)
            if (is.na(test.invisible[1])) {
                points(coord.ind, pch = 16, col = col.ind)
                if (lab.ind)
                  text(coord.ind[, 1], y = coord.ind[, 2], labels = rownames(coord.ind),
                    pos = 3, col = col.ind, cex = cex)
            }
            if (is.na(test.invisible[2])) {
                points(coord.var[, 1], y = coord.var[, 2], pch = 17,
                  col = col.var, cex = cex)
                if (lab.var)
                  text(coord.var[, 1], y = coord.var[, 2], labels = rownames(coord.var),
                    pos = 3, col = col.var, cex = cex)
                if ((habillage != "none") & (habillage != "quali") &
                  (is.na(test.invisible[1]) | is.na(test.invisible[2])))
                  legend("topleft", legend = levels(res.mca$call$X[,
                    habillage]), text.col = 1:n.mod, cex = 0.8)
            }
        }
    }
    group = nrow(res.mca$group$coord)
    if (choix == "group") {
        if (new.plot)
            dev.new()
        if (is.null(title))
            title <- "Consumers representation"
        else sub.title <- "Consumers representation"
        coord.actif <- res.mca$group$coord[, axes]
        plot(coord.actif, xlab = paste("Dim ", axes[1], " (",
            signif(res.mca$eig[axes[1], 2], 4), "%)", sep = ""),
            ylab = paste("Dim ", axes[2], " (", signif(res.mca$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = c(0, 1), ylim = c(0,
                1), pch = 17, cex = cex, main = title, cex.main = cex *
                1.2, asp = 1)
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
        if (lab.grpe)
            text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif),
                pos = 3)
    }
}
