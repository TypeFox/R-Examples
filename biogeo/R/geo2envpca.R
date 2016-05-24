geo2envpca <-
function (edat, g1, group1 = "Species", group2 = "", world, scaling = 1, 
    vars = c("AMT", "AP", "MTCM", "MTWM", "PWQ", "PCQ"), showrecord = "", 
    ext = c(-180, 180, -60, 90)) 
{
    fg <- dev.list()
    fig1 <- fg[1]
    fig2 <- fg[2]
    fig3 <- fg[3]
    cn <- names(edat)
    f3 <- match(group1, cn)
    ff <- which(edat[, f3] == g1 & edat$Exclude == 0)
    nv <- length(vars)
    vrs <- rep(0, nv)
    for (i in 1:nv) {
        vrs[i] <- match(vars[i], cn)
    }
    env <- edat[ff, vrs]
    envna <- env * 0
    for (i in 1:ncol(env)) {
        envna[, i] <- is.na(env[i]) * 1
    }
    x <- coord2numeric(edat$x[ff])
    y <- coord2numeric(edat$y[ff])
    nas <- cbind(is.na(x) * 1, is.na(y) * 1, envna)
    nas1 <- rowSums(nas)
    f <- which(nas1 >= 1)
    ff2 <- ff[f]
    edat$Exclude[ff2] <- 2
    ff <- which(edat[, f3] == g1 & edat$Exclude == 0)
    env <- edat[ff, vrs]
    x <- coord2numeric(edat$x[ff])
    y <- coord2numeric(edat$y[ff])
    sc <- scaling
    pca <- rda(env, scale = TRUE, scaling = sc)
    ei <- pca$CA$eig
    va <- round((ei/sum(ei)) * 100, 1)
    sites <- scores(pca, choices = 1:2, scaling = sc, display = "sites")
    vars <- scores(pca, choices = 1:2, scaling = sc, display = "species")
    if (length(ext) == 4) {
        xlm <- ext[1:2]
        ylm <- ext[3:4]
    }
    else {
        mnx <- min(x)
        mxx <- max(x)
        mny <- min(y)
        mxy <- max(y)
        rgx <- (mxx - mnx) * 0.1
        rgy <- (mxy - mny) * 0.1
        xlm <- c(mnx - rgx, mxx + rgx)
        ylm <- c(mny - rgy, mxy + rgy)
    }
    #bringToTop(dev.set(fig1), stay = T)
    dev.set(fig1)
    par(mai = c(1.36, 0.2, 1.093333, 0.2))
    plot(world, border = "gray", xlim = xlm, ylim = ylm)
    box(which = "plot", lty = "solid", col = "brown")
    points(x, y, xlab = "x-coordinate", ylab = "y-coordinate")
    if (nchar(showrecord) > 0) {
        fid <- which(edat$ID == showrecord)
        if (length(fid) > 0) {
            text(edat$x[fid], edat$y[fid], labels = showrecord, 
                adj = 1, pos = 4, cex = 0.7)
            points(edat$x[fid], edat$y[fid], col = "green")
        }
    }
    if (nchar(group2) > 0) {
        f4 <- match(group2, cn)
        g2 <- unique(edat[, f4])
        ng2 <- length(g2)
        if (ng2 > 0) {
            cols <- rep(1, ng2)
            for (j in 1:ng2) {
                fg2 <- which(edat[ff, f4] == g2[j])
                points(x[fg2], y[fg2], col = j + 2)
                cols[j] <- j + 2
            }
            legend("bottomleft", legend = as.character(g2), col = cols, 
                pch = 18)
        }
    }
    mtext(g1, side = 3, line = 1, adj = 0.5)
    #bringToTop(dev.set(fig2), stay = T)
    dev.set(fig2)
    par(mai = c(1.36, 1.093333, 1.093333, 0.5))
    xmin <- min(min(sites[, 1]), min(vars[, 1]))
    xmax <- max(max(sites[, 1]), max(vars[, 1]))
    ymin <- min(min(sites[, 2]), min(vars[, 2]))
    ymax <- max(max(sites[, 2]), max(vars[, 2]))
    xlabel <- paste("PC1", " (", va[1], " %)", sep = "")
    ylabel <- paste("PC2", " (", va[2], " %)", sep = "")
    plot(vars[, 1], vars[, 2], xlim = c(xmin, xmax), ylim = c(ymin, 
        ymax), type = "n", xlab = xlabel, ylab = ylabel)
    arrows(0, 0, vars[, 1], vars[, 2], length = 0, lty = 3, col = "gray")
    points(sites[, 1], sites[, 2])
    txt <- row.names(vars)
    text(vars[, 1], vars[, 2], labels = txt, adj = 0.5, offset = 1)
    box(which = "plot", lty = "solid", col = "blue")
    if (nchar(group2) > 0) {
        if (ng2 > 0) {
            for (j in 1:ng2) {
                fg2 <- which(edat[ff, f4] == g2[j])
                points(sites[fg2, 1], sites[fg2, 2], col = j + 
                  2)
            }
        }
    }
    f <- 0
    while (f != 4) {
        #bringToTop(dev.set(fig3), stay = T)
      dev.set(fig3)
        par(mai = c(0, 0, 0, 0))
        plot(1, 1, type = "n", xlab = "", ylab = "", xaxt = "n", 
            yaxt = "n", xlim = c(0, 5), ylim = c(0, 4))
        x1 <- rep(0.5, 4)
        y1 <- (4:1) - 0.5
        points(x1, y1, pch = 18, col = c("brown", "blue", "blue", 
            "black"))
        text(x1, y1, labels = c("Query Geographical", "Query Environment", 
            "Exclude Environment", "Exit"), adj = 0, pos = 4, 
            cex = 0.7)
        f <- identify(x1, y1, n = 1, plot = FALSE)
        if (f == 1) {
            #bringToTop(dev.set(fig1), stay = T)
            dev.set(fig1)
            fi <- identify(x, y, n = 1, plot = FALSE)
            points(x[fi], y[fi], pch = 18, col = "red")
            fdat <- ff[fi]
            fid <- edat$ID[fdat]
            text(x[fi], y[fi], labels = fid, adj = 1, pos = 4, 
                cex = 0.7)
            #bringToTop(dev.set(fig2))
            dev.set(fig2)
            points(sites[fi, 1], sites[fi, 2], pch = 18, col = "red")
            text(sites[fi, 1], sites[fi, 2], labels = fid, adj = 1, 
                pos = 4, cex = 0.7)
        }
        if (f == 2) {
            #bringToTop(dev.set(fig2), stay = T)
            dev.set(fig2)
            fi <- identify(sites[, 1], sites[, 2], n = 1, plot = FALSE)
            points(sites[fi, 1], sites[fi, 2], pch = 18, col = "red")
            fdat <- ff[fi]
            fid <- edat$ID[fdat]
            text(sites[fi, 1], sites[fi, 2], labels = fid, adj = 1, 
                pos = 4, cex = 0.7)
            #bringToTop(dev.set(fig1))
            dev.set(fig1)
            points(x[fi], y[fi], pch = 18, col = "red")
            text(x[fi], y[fi], labels = fid, adj = 1, pos = 4, 
                cex = 0.7)
        }
        if (f == 3) {
            #bringToTop(dev.set(fig2))
          dev.set(fig2)
            fi <- identify(sites[, 1], sites[, 2], n = 1, plot = FALSE)
            points(sites[fi, 1], sites[fi, 2], pch = 18, col = "black")
            fdat <- ff[fi]
            fid <- edat$ID[fdat]
            text(sites[fi, 1], sites[fi, 2], labels = fid, adj = 1, 
                pos = 4, cex = 0.7)
            edat$Exclude[fdat] <- 5
        }
        if (f == 4) {
            graphics.off()
        }
    }
    if (any(edat$Exclude == 2)) {
        f2 <- which(edat$Exclude == 2)
        edat$Exclude[f2] <- 0
    }
    return(edat)
}
