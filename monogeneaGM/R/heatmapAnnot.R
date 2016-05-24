#' Annotated heat map
#'
#' This function produces a heat map with hierarchical clustering of samples using multivariate size or shape data.
#' @param x a matrix with rows representing species samples and columns representing morphometrical variables of
#' interest
#' @param labcol a character vector giving the color annotation for the species 
#' @param xlab title for x-axis
#' @param ylab title for y-axis
#' @param genus single character abbreviation for genus
#' @param rowlab if FALSE, suppresses row labels on the heat map
#' @param pt if the rows have been ranked, two lines are drawn - one at the \code{pt}x100th percentile (bottom), 
#' the other at the (1-\code{pt})x100th percentile (top); set to 0 to disable
#' @param dist.method a distance metric implemented in the \code{dist} function; defaults to \code{manhattan}
#' @param clust.method a clustering algorithm implemented in the \code{hclust} function; defaults to \code{ward.D}
#' @details This function customizes the \code{heatmap.2} function in the \code{gplots} package (Version 2.17.0). Row standardization is switched on, and a color legend for species is given in the left panel.  
#' @seealso \code{\link{dist}}, \code{\link{hclust}}
#' @author Tsung Fei Khang \email{tfkhang@@um.edu.my}
#' @references Khang TF, Soo OYM, Tan WB, Lim LHS. (2016). Monogenean anchor morphometry: systematic value, phylogenetic signal, and evolution. PeerJ 4:e1668.
#'
#' Warnes GR, Bolker B, Bonebakker L, Gentleman R, Huber W, Liaw A, Lumley T,Maechler M, Magnusson M, Moeller S, 
#' Schwartz M, Venables B. (2015). gplots:Various R programming tools for plotting data. R package version 2.17.0.
#' Available at: http://CRAN.R-project.org/package=gplots.
#' 
#' @examples
#' library(gplots)
#'
#' data(ligophorus_shape)
#' data(spcolmap)
#'
#' dendrogram <- hclust(dist(ligophorus_shape))
#'
#' #check dendrogram and note cut-off for the two main clades
#' plot(dendrogram)
#'
#' clade_id <- cutree(dendrogram, h=0.5)
#'
#' f_s <- numeric(ncol(ligophorus_shape))
#' 
#' for(i in 1:ncol(ligophorus_shape)){
#' dat <- stack(ligophorus_shape[,i])
#' #replace species label with clade label
#' dat[,2] <- clade_id
#' f_s[i] <- t.test(values~ind, data=dat)$statistic
#' }
#'
#' rank_s <- order(f_s, decreasing=TRUE)
#'
#' heatmapAnnot(ligophorus_shape[,rank_s],labcol=spcolmap$color,
#' xlab="Specimens", genus="L. ")
#'

heatmapAnnot <- function (x, labcol, xlab = "", ylab = "", genus = "", rowlab=TRUE, pt = 0, dist.method="manhattan", clust.method="ward.D")
{

    splabel <- as.character(levels(as.factor(rownames(x))))
    colorlabel <- unlist(mapply(function(k) rep(labcol[k], sum(rownames(x) == 
        splabel[k])), k = 1:length(splabel)))
    species_full <- mapply(function(k) paste(c(genus, splabel[k]), 
        collapse = ""), k = 1:length(splabel))
    f <- vector("expression", length(species_full))
    for (s in 1:length(species_full)) {
        f[[s]] <- substitute(italic(nn), list(nn = species_full[s]))
    }
    lmat <- matrix(c(5, 3, 4, 0, 3, 1, 6, 3, 2), byrow = TRUE, 
        ncol = 3)
    lhei <- c(1.5, 0.2, 4)
    lwid <- c(2.5, 0.2, 4)
    myplot <- function() {
        oldpar <- par("mar")
        par(mar = c(5.1, 4.1, 0.5, 0.5))
        plot(1:10, 1:10, pch = "", bty = "n", xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "")
        text(5, 10, "Species", cex = 1.5)
        legend(1, 9, pch = rep(16, length(splabel)), pt.cex = 2, 
            col = labcol, f)
        legend(1, 9, pch = rep(1, length(splabel)), pt.cex = 2, 
            f)
    }
    rwb <- c("#99000D", "#FB6A4A", "white", "#6BAED6", "#084594")
    rwbtones <- colorRampPalette(rwb, space = "Lab")(100)

    if(rowlab == FALSE){
    heatmap.2(t(x), col = rwbtones, scale = "row", Colv = TRUE, 
        hclustfun = function(k) hclust(k, method = clust.method), dendrogram = "column", 
        distfun = function(k) dist(k, method = dist.method), 
        symkey = FALSE, Rowv = FALSE, keysize = 0.25, ColSideColors = colorlabel, 
        margins = c(3, 3), density.info = "none", trace = "none", 
        labRow = FALSE, labCol = FALSE, xlab = xlab, ylab = ylab, 
        lmat = lmat, lhei = lhei, lwid = lwid, key = TRUE, extrafun = myplot, 
        colsep = c(0, nrow(x)), rowsep = c(0, round(ncol(x) * 
            pt), round(ncol(x) * (1 - pt)), ncol(x)), sepcolor = "black", 
        sepwidth = c(0.001, 0.001))
	}

	else{
    heatmap.2(t(x), col = rwbtones, scale = "row", Colv = TRUE, 
        hclustfun = function(k) hclust(k, method = clust.method), dendrogram = "column", 
        distfun = function(k) dist(k, method = dist.method), 
        symkey = FALSE, Rowv = FALSE, keysize = 0.25, ColSideColors = colorlabel, 
        margins = c(3, 3), density.info = "none", trace = "none", 
        labRow = colnames(x), labCol = FALSE, xlab = xlab, ylab = ylab, 
        lmat = lmat, lhei = lhei, lwid = lwid, key = TRUE, extrafun = myplot, 
        colsep = c(0, nrow(x)), rowsep = c(0, round(ncol(x) *  pt), 
        round(ncol(x) * (1 - pt)), ncol(x)), sepcolor = "black", 
        sepwidth = c(0.001, 0.001))
	}

}


