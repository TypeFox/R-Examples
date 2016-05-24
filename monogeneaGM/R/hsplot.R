#' Scatter plot with marginal histograms
#'
#' This function creates a scatter plot that is augmented with marginal histograms (with smoothed splines) on the x and y-axes.
#' @param dat a matrix with rows representing species samples and two columns representing the variables of
#' interest (typically a selected principal component of size and shape variables)
#' @param seqx range of values on the x-axis, for histogram plotting
#' @param seqy range of values on the y-axis, for histogram plotting
#' @param ybound range of values on the y-axis
#' @param xbound range of values on the x-axis
#' @param y.label title for y-axis
#' @param x.label title for x-axis
#' @param labcol a character vector giving the color annotation for species names in alphabetical order
#' @param type a character that controls where the marginal histograms are located: right-up ("ru"), left-up ("lu"), left-down ("ld"), right-down ("rd")
#' @param sep1 tick mark separation for x-axis
#' @param sep2 tick mark separation for y-axis
#' @param phylo if TRUE, the centroids of species are joined into a tree-like object corresponding to their phylogeny given in \code{phy}
#' @param phy a tree object of class \code{phylo}
#' @param pointscale a constant specifying the symbol size of the centroids
#' @param supp.hist.y if TRUE, suppresses the plotting of marginal histogram at the y-axis
#' @param supp.hist.x if TRUE, suppresses the plotting of marginal histogram at the x-axis
#' @param y.axis.label if FALSE, suppresses the y-axis title
#' @param x.axis.label if FALSE, suppresses the x-axis title
#' @param addline an optional numeric vector of length two specifying the intercept and slope of linear line 
#" to be superimposed on the main plot
#' @details This function is similar to \code{tpColorPlot2d}, and is mainly useful for its additional
#' option for visualising the marginal distributions. The \code{ape} and \code{phytools} packages are required.
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Paradis E, Claude J & Strimmer K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
#'
#' Revell LJ. (2012). phytools: An R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217-223.
#' @examples
#' library(phytools)
#'
#' data(pwed_pd)
#' data(ligotree)
#' data(spcolmap)
#'
#' #principal component analysis of size variables from ventral and dorsal anchors
#' pca2d_pwedv <- pca2d(pwed_pd[,1:55], sgn=-1, labcol=spcolmap$color, 
#' phylo=TRUE, phy=ligotree, genus="L. ", bound.y=c(-0.2, 0.2),
#' bound.x1=c(-0.2,0.2), bound.x2 = c(-0.2,0.2))
#"
#' pca2d_pwedd <- pca2d(pwed_pd[,56:110], sgn=-1, labcol=spcolmap$color,
#' phylo=TRUE, phy=ligotree, genus="L. ", bound.y=c(-0.2, 0.2),
#' bound.x1=c(-0.2,0.2), bound.x2 = c(-0.2,0.2))
#'
#' #comparing size PC1 between ventral and dorsal anchors
#' hsplot(cbind(pca2d_pwedv$scores[,1], pca2d_pwedd$scores[,1]), 
#' seqx=seq(-0.2,0.15,0.01), seqy=seq(-0.2,0.15,0.01),
#' xbound=c(-0.2,0.15), ybound=c(-0.2,0.15), y.label="Dorsal size PC1", 
#' x.label="Ventral size PC1", labcol=spcolmap$color, type="ru", 
#' sep1=0.1, sep2=0.1, phylo=FALSE, addline=c(0,1))
#' 

hsplot <- function (dat, seqx, seqy, ybound, xbound, y.label = "", x.label = "", 
    labcol, type = "ru", sep1, sep2, phylo = FALSE, phy, pointscale = 1.5, 
    supp.hist.y = FALSE, supp.hist.x = FALSE,  y.axis.label = TRUE,
    x.axis.label = TRUE, addline=NULL) {

    left.flip <- 1
    down.flip <- 1
    if (type == "ru") {
        nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), 
            c(3, 1), c(1, 3), TRUE)
        layout.show(nf)
        par(mar = c(5, 4, 0.5, 0.2))
        dir.x <- 1
        dir.y <- 2
        prettyset.x <- c(0, 4, 1, 0)
        prettyset.y <- c(5, 1, 1, 1)
    }
    else if (type == "lu") {
        nf <- layout(matrix(c(0, 2, 3, 1), 2, 2, byrow = TRUE), 
            c(1, 3), c(1, 3), TRUE)
        layout.show(nf)
        left.flip <- -1
        par(mar = c(5, 0.2, 0.5, 4))
        dir.x <- 1
        dir.y <- 4
        prettyset.x <- c(0, 0, 1, 4)
        prettyset.y <- c(5, 1, 1, 1)
    }
    else if (type == "ld") {
        nf <- layout(matrix(c(3, 1, 0, 2), 2, 2, byrow = TRUE), 
            c(1, 3), c(3, 1), TRUE)
        layout.show(nf)
        down.flip <- -1
        left.flip <- -1
        par(mar = c(0.5, 0.2, 5, 4))
        dir.x <- 3
        dir.y <- 4
        prettyset.x <- c(3, 0, 0, 4)
        prettyset.y <- c(1, 1, 5, 1)
    }
    else if (type == "rd") {
        nf <- layout(matrix(c(1, 3, 2, 0), 2, 2, byrow = TRUE), 
            c(3, 1), c(3, 1), TRUE)
        layout.show(nf)
        down.flip <- -1
        par(mar = c(0.5, 4, 5, 0.2))
        dir.x <- 3
        dir.y <- 2
        prettyset.x <- c(3, 4, 0, 0.5)
        prettyset.y <- c(0, 1, 5, 1)
    }
    plot(dat[, 1], dat[, 2], cex = 1, col = "white", ylim = ybound, 
        xlim = xbound, ylab = "", xlab = "", xaxt = "n", yaxt = "n")
    if(is.null(addline)==FALSE){
	abline(a=addline[1], b=addline[2])
	}
    if (y.axis.label == TRUE) {
        axis(dir.y, round(seq(seqy[1], seqy[length(seqy)], sep2), 
            2), round(seq(seqy[1], seqy[length(seqy)], sep2), 
            2))
        mtext(y.label, side = dir.y, line = 2.5)
    }
    if (x.axis.label == TRUE) {
        axis(dir.x, round(seq(seqx[1], seqx[length(seqx)], sep1), 
            2), round(seq(seqx[1], seqx[length(seqx)], sep1), 
            2))
        mtext(x.label, side = dir.x, line = 3)
    }
    species <- levels(as.factor(rownames(dat)))
    nsp <- length(species)
    if (phylo == TRUE) {
        mean.mat <- matrix(0, nsp, ncol(dat))
        for (i in 1:ncol(dat)) {
            temp <- stack(dat[, i])
            mean.mat[, i] <- tapply(temp[, 1], temp[, 2], mean)
        }
        rownames(mean.mat) <- species
        anc.coord <- matrix(0, phy$Nnode, 2)
        for (i in 1:2) {
            anc.coord[, i] <- fastAnc(phy, mean.mat[, i])
        }
        all.nodes <- mean.mat
        all.nodes <- all.nodes[phy$tip.label, ]
        all.nodes <- rbind(all.nodes, anc.coord)
        for (i in 1:nrow(phy$edge)) {
            lines(all.nodes[phy$edge[i, ], 1], all.nodes[phy$edge[i, 
                ], 2])
        }
        for (i in 1:nrow(anc.coord)) {
            points(anc.coord[i, 1], anc.coord[i, 2], pch = 16, 
                cex = 0.5 * pointscale, col = "white")
            points(anc.coord[i, 1], anc.coord[i, 2], cex = 0.5 * 
                pointscale)
        }
    }
    for (i in 1:nrow(dat)) {
        sp <- which(rownames(dat) == species[i])
        red.green.blue <- as.numeric(col2rgb(labcol[i]))
        points(dat[sp, 1], dat[sp, 2], col = rgb(red = red.green.blue[1], 
            green = red.green.blue[2], blue = red.green.blue[3], 
            alpha = 60, maxColorValue = 255), pch = 16, cex = 1)
        points(mean(dat[sp, 1]), mean(dat[sp, 2]), col = rgb(red = red.green.blue[1], 
            green = red.green.blue[2], blue = red.green.blue[3], 
            alpha = 255, maxColorValue = 255), pch = 16, cex = pointscale)
        points(mean(dat[sp, 1]), mean(dat[sp, 2]), pch = 1, cex = pointscale)
    }
    xhist <- hist(dat[, 1], breaks = seqx, plot = FALSE)
    yhist <- hist(dat[, 2], breaks = seqy, plot = FALSE)
    if (supp.hist.x == TRUE) {
        plot(1:1, 1:1, pch = "", bty = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n")
    }
    else if (supp.hist.x == FALSE) {
        par(mar = prettyset.x)
        cd <- down.flip * xhist$counts
        xdat <- 1:length(cd) - 0.5
        barplot(cd, axes = TRUE, space = 0, main = "", col = "gray90", 
            plot = TRUE, xaxt = "n", yaxt = "n")
        smoothingSpline <- smooth.spline(xdat, cd, spar = 0.3)
        xl <- seq(min(xdat), max(xdat), (max(xdat) - min(xdat))/1000)
        lines(xl, predict(smoothingSpline, xl)$y, col = "red", 
            lwd = 1)
    }
    if (supp.hist.y == TRUE) {
        plot(1:1, 1:1, pch = "", bty = "n", xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n")
    }
    else if (supp.hist.y == FALSE) {
        par(mar = prettyset.y)
        cd <- left.flip * yhist$counts
        xdat <- 1:length(cd) - 0.5
        barplot(cd, axes = TRUE, space = 0, main = "", horiz = TRUE, 
            col = "gray90", plot = TRUE, xaxt = "n", yaxt = "n")
        smoothingSpline <- smooth.spline(xdat, cd, spar = 0.3)
        xl <- seq(min(xdat), max(xdat), (max(xdat) - min(xdat))/1000)
        lines(predict(smoothingSpline, xl)$y, xl, , col = "red", 
            lwd = 1)
    }
}


