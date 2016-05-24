plot.fasnt=function (x,choix="ind", axes = c(1, 2), xlim = NULL, ylim = NULL, invisible = NULL,
    col.ind = "blue", col.var = "red", lab.ind=TRUE,lab.var=TRUE, lab.coord=TRUE, lab.partial=TRUE,
    cex = 1, lab.grpe = TRUE, title = NULL, habillage = "none", palette = NULL,
    new.plot = TRUE, ...)
{
    res.fasnt <- x
    if (!inherits(res.fasnt, "fasnt"))
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
#    lab.ind <- lab.var <- FALSE
#    if (length(label) == 1 && label == "all")
#        lab.ind <- lab.var <- TRUE
#    if ("ind" %in% label)
#        lab.ind <- TRUE
#    if ("var" %in% label)
#        lab.var <- TRUE
if (choix=="ind"){
    test.invisible <- vector(length = 5)
    if (!is.null(invisible)) {
        test.invisible[1] <- match("ind", invisible)
        test.invisible[2] <- match("var", invisible)
    }
    else test.invisible <- rep(NA, 5)
    coord.var <- res.fasnt$quali.var$coord[, axes]
    coord.ind <- res.fasnt$ind$coord[, axes]
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
    if (habillage != "none") {
        if (!is.factor(res.fasnt$call$X[, habillage]))
            stop("The variable ", habillage, " is not qualitative")
        col.ind <- as.numeric(as.factor(res.fasnt$call$X[, habillage]))
        n.mod <- nlevels(as.factor(res.fasnt$call$X[, habillage]))
    }

    sub.titre <- NULL
    titre = title
    if (is.null(title))
        titre <- "Sorted Napping Task factor map"
    else sub.titre <- "Sorted Napping Task factor map"
    if (is.na(test.invisible[1]) | is.na(test.invisible[2]) |
        is.na(test.invisible[4]) | is.na(test.invisible[5])) {
        if (new.plot)
            dev.new(width = 8, height = 8)
        plot(0, 0, main = titre, xlab = paste("Dim ", axes[1],
            " (", signif(res.fasnt$eig[axes[1], 2], 4), "%)", sep = ""),
            ylab = paste("Dim ", axes[2], " (", signif(res.fasnt$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = xlim, ylim = ylim,
            col = "white", asp = 1, cex = cex)
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
            legend("topleft", legend = levels(res.fasnt$call$X[,
                habillage]), text.col = 1:n.mod, cex = 0.8)
    }
    }
}

group=nrow(res.fasnt$group$coord[[2]])

if (choix == "group") {
        if (new.plot)
            dev.new(width = 8, height = 8)
        if (is.null(title))
            title <- "Subjects representation"
        else sub.title <- "Subjects representation"
        coord.actif <- res.fasnt$group$coord[[2]][, axes]

        plot(coord.actif, xlab = paste("Dim ", axes[1], " (", signif(res.fasnt$eig[axes[1],
                2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(res.fasnt$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = c(0,
            1), ylim = c(0, 1), pch = 17,
            cex = cex, main = title, cex.main = cex * 1.2, asp = 1)
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
         if (lab.grpe)
            text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif),
                pos = 3)
}

if (choix == "coord") {
        if (new.plot)
        dev.new(width = 8, height = 8)
        if (is.null(title))
            title <- "Napping data representation"
        else sub.title <- "Napping representation"
        lab.x <- paste("Dim ", axes[1], " (", signif(res.fasnt$eig[axes[1],2], 4), " %)", sep = "")
        lab.y <- paste("Dim ", axes[2], " (", signif(res.fasnt$eig[axes[2],2], 4), " %)", sep = "")
        plot(0, 0, xlab = lab.x, ylab = lab.y,xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1),main = title, col = "white",asp = 1)
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
        x.cercle <- seq(-1, 1, by = 0.01)
        y.cercle <- sqrt(1 - x.cercle^2)
        lines(x.cercle, y = y.cercle)
        lines(x.cercle, y = -y.cercle)
        abline(v = 0, lty = 2)
        abline(h = 0, lty = 2)

        coord.var <- res.fasnt$quanti.var$cor[, axes]
        for (v in 1:nrow(coord.var)) {
        arrows(0, 0, coord.var[v, 1], coord.var[v,2], length = 0.1, angle = 15, code = 2)
        if (abs(coord.var[v, 1]) > abs(coord.var[v,2])) {
        if (coord.var[v, 1] >= 0)
        pos <- 4
        else pos <- 2
          }
         else {
        if (coord.var[v, 2] >= 0)
          pos <- 3
        else pos <- 1
        }
        if (lab.coord){
        text(coord.var[v, 1], y = coord.var[v, 2],labels = rownames(coord.var)[v], pos = pos)}

}
}


if (choix=="partial"){
        if (new.plot)
        dev.new(width = 8, height = 8)
        if (is.null(title))
            title <- "Partial representations"
        else sub.title <- "Partial representations"

        coord.ind <- res.fasnt$ind$coord[, axes]

        if (is.null(xlim)) {
        inter <- res.fasnt$ind$partial[[2]][,, 1][,axes]
        for (i in 2:nrow(res.fasnt$group$coord[[2]])) inter <- rbind(inter, res.fasnt$ind$partial[[2]][,, i][,axes])
        xmin <- min(res.fasnt$ind$coord[, axes[1]],inter[, axes[1]])
        xmax <- max(res.fasnt$ind$coord[, axes[1]],inter[, axes[1]])
        ymin <- min(res.fasnt$ind$coord[, axes[2]],inter[, axes[2]])
        ymax <- max(res.fasnt$ind$coord[, axes[2]],inter[, axes[2]])
        xlim <- c(xmin, xmax) * 1.1
        ylim <- c(ymin, ymax) * 1.1 }

        plot(0, 0, main = title, xlab = paste("Dim ", axes[1],
            " (", signif(res.fasnt$eig[axes[1], 2], 4), "%)", sep = ""),
            ylab = paste("Dim ", axes[2], " (", signif(res.fasnt$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = xlim, ylim = ylim,
            col = "white", asp = 1, cex = cex)
        abline(v = 0, lty = 2)
        abline(h = 0, lty = 2)

        points(res.fasnt$ind$coord[, axes], pch = 15,col=1:nrow(res.fasnt$ind$coord))
        if (lab.ind)
        text(res.fasnt$ind$coord[, axes[1]],res.fasnt$ind$coord[, axes[2]],col=1:nrow(res.fasnt$ind$coord), rownames(res.fasnt$ind$coord),pos = 3, offset = 0.2, cex = 0.8*cex)

        for (j in 1:nrow(res.fasnt$group$coord[[2]])) {
        points(res.fasnt$ind$partial[[2]][,, j][,axes], col=rep(1:nrow(res.fasnt$ind$coord),times=nrow(res.fasnt$group$coord[[2]])),pch = 20, cex = 0.8*cex)
        if (lab.partial)
        text(res.fasnt$ind$partial[[2]][,, j][,axes], col=rep(1:nrow(res.fasnt$ind$coord),times=nrow(res.fasnt$group$coord[[2]])),labels=rownames(res.fasnt$group$coord[[2]])[j],pos=3,cex=0.5*cex)
        for (i in 1:nrow(res.fasnt$ind$partial[[2]]))
        lines(c(res.fasnt$ind$coord[i,axes[1]], res.fasnt$ind$partial[[2]][,,j][i, axes[1]]), c(res.fasnt$ind$coord[i, axes[2]], res.fasnt$ind$partial[[2]][,, j][i,axes[2]]),col=i)
        }


}

}
