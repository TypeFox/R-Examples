plot.fahst=function (x,choix="ind", axes = c(1, 2), xlim = NULL, ylim = NULL, invisible = NULL,
    col.ind = "blue", col.var = "red", lab.ind=TRUE,lab.var=TRUE,
    cex = 1,lab.lev=TRUE, lab.grpe = TRUE, title = NULL, habillage = "none", habillage.lev= "none", traj=FALSE, palette = NULL,
    new.plot = TRUE, ...)
{
    res.fahst <- x
    if (!inherits(res.fahst, "fahst"))
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
    test.invisible <- vector(length = 5)
    if (!is.null(invisible)) {
        test.invisible[1] <- match("ind", invisible)
        test.invisible[2] <- match("var", invisible)
    }
    else test.invisible <- rep(NA, 5)
    coord.var <- res.fahst$var$coord[, axes]
    coord.ind <- res.fahst$ind$coord[, axes]
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
        for (j in res.fahst$call$quali) {
            col.var <- c(col.var, rep(aux, nlevels(res.fahst$call$X[,
                j])))
            aux = aux + 1
        }
    }
    if ((habillage != "none") & (habillage != "quali")) {
        if (!is.factor(res.fahst$call$X[, habillage]))
            stop("The variable ", habillage, " is not qualitative")
        col.ind <- as.numeric(as.factor(res.fahst$call$X[, habillage]))
        n.mod <- nlevels(as.factor(res.fahst$call$X[, habillage]))
    }
    if (choix=="ind" | choix=="var"){
    sub.titre <- NULL
    titre = title
    if (is.null(title))
        titre <- "Hierarchical Sorting Task factor map"
    else sub.titre <- "Hierarchical Sorting Task factor map"
    if (is.na(test.invisible[1]) | is.na(test.invisible[2]) |
        is.na(test.invisible[4]) | is.na(test.invisible[5])) {
        if (new.plot)
            dev.new(width = 8, height = 8)
        plot(0, 0, main = titre, xlab = paste("Dim ", axes[1],
            " (", signif(res.fahst$eig[axes[1], 2], 4), "%)", sep = ""),
            ylab = paste("Dim ", axes[2], " (", signif(res.fahst$eig[axes[2],
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
            legend("topleft", legend = levels(res.fahst$call$X[,
                habillage]), text.col = 1:n.mod, cex = 0.8)
    }
    }
    }

group=nrow(res.fahst$group$coord)
    if (choix == "group") {
        if (new.plot)
            dev.new(width = 8, height = 8)
        if (is.null(title))
            title <- "Subjects representation"
        else sub.title <- "Subjects representation"
        coord.actif <- res.fahst$group$coord[, axes]

        plot(coord.actif, xlab = paste("Dim ", axes[1], " (", signif(res.fahst$eig[axes[1],
                2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(res.fahst$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = c(0,
            1), ylim = c(0, 1), pch = 17,
            cex = cex, main = title, cex.main = cex * 1.2, asp = 1)
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
         if (lab.grpe)
            text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif),
                pos = 3)
}


if (choix == "level") {
        if (new.plot)
            dev.new(width = 8, height = 8)
        if (is.null(title))
            title <- "Levels representation"
        else sub.title <- "Levels representation"
        coord.actif <- res.fahst$var$coord.lev[, axes]
        
        if (habillage.lev == "none"){
        plot(coord.actif, xlab = paste("Dim ", axes[1], " (", signif(res.fahst$eig[axes[1],
                2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(res.fahst$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = c(0,
            1), ylim = c(0, 1), pch = 17,
            cex = cex, main = title, cex.main = cex * 1.2, asp = 1)
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
         if (lab.lev)
            text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif),
                pos = 3)
        
          if (traj == TRUE){
            subj=0 
            for (j in 1:group){
            if (res.fahst$call$group[j]!=1){
            for (i in 1:(res.fahst$call$group[j]-1)){
            lines(x=res.fahst$var$coord.lev[(subj+i):(subj+i+1),axes[1]],y=res.fahst$var$coord.lev[(subj+i):(subj+i+1),axes[2]])}}
            subj=subj+res.fahst$call$group[j]
            }
          }
        }
        
        if (habillage.lev == "subject"){
        plot(coord.actif, xlab = paste("Dim ", axes[1], " (", signif(res.fahst$eig[axes[1],
                2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(res.fahst$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = c(0,
            1), ylim = c(0, 1), pch = 17,
            cex = cex, main = title, cex.main = cex * 1.2, asp = 1,col=rep(1:group,times=res.fahst$call$group))
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
         if (lab.lev)
            text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif),
                pos = 3,col=rep(1:group,times=res.fahst$call$group))
         if (traj == TRUE){
            subj=0 
            for (j in 1:group){
            if (res.fahst$call$group[j]!=1){
            for (i in 1:(res.fahst$call$group[j]-1)){
            lines(x=res.fahst$var$coord.lev[(subj+i):(subj+i+1),axes[1]],y=res.fahst$var$coord.lev[(subj+i):(subj+i+1),axes[2]],col=j)}}
            subj=subj+res.fahst$call$group[j]
            }
          }
         legend("topleft", legend = rownames(res.fahst$group$coord), text.col = 1:group, cex = 0.8)
        }

lev=c()
for (i in 1:group){
lev=c(lev,1:res.fahst$call$group[i])}
        
        if (habillage.lev == "level"){
        plot(coord.actif, xlab = paste("Dim ", axes[1], " (", signif(res.fahst$eig[axes[1],
                2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(res.fahst$eig[axes[2],
                2], 4), "%)", sep = ""), xlim = c(0,
            1), ylim = c(0, 1), pch = 17,
            cex = cex, main = title, cex.main = cex * 1.2, asp = 1,col=lev)
        title(cex.sub = cex, font.sub = 2, col.sub = "steelblue4",
            adj = 0, line = 3.8)
            
         if (lab.lev)
            text(coord.actif[, 1], y = coord.actif[, 2], labels = rownames(coord.actif),
                pos = 3,col=lev)
                
         if (traj == TRUE){
            subj=0 
            for (j in 1:group){
            if (res.fahst$call$group[j]!=1){
            for (i in 1:(res.fahst$call$group[j]-1)){
            lines(x=res.fahst$var$coord.lev[(subj+i):(subj+i+1),axes[1]],y=res.fahst$var$coord.lev[(subj+i):(subj+i+1),axes[2]])}}
            subj=subj+res.fahst$call$group[j]
            }
          }
         legend("topleft", legend = paste("Level",levels(as.factor(lev))), text.col = 1:max(lev), cex = 0.8)
        }
}

}
