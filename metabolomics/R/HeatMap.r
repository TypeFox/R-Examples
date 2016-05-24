HeatMap <- function(inputdata, 
    colramp=redgreen(75),    
    scale=c("row", "column", "none"), 
    dendrogram=c("column", "row", "both", "none"), distmethod="euclidean",
    aggmethod="complete", margins=c(5, 5), key=TRUE,
    keysize=1.5, cexRow=0.5, ColSideColors=NULL, ...)
{
    # Create dendogram
    scale <- match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    
    # Create groups information, colour scale
    groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
    unique.groups <- levels(groups)
    if(is.null(ColSideColors))  {
        cols <- ColList(length(unique.groups))
        ColSideColors <- c(rep(NA, length(rownames(inputdata))))
        for (ii in 1:length(inputdata[, 1])) {
            selected <- which(unique.groups == inputdata[, 1][ii])
            ColSideColors[ii] <- cols[selected]
        }
    }
    
    # Create heatmap
    inputdata <- t(editcolnames(inputdata)[, -1])
    heatmap_fn(inputdata, col=colramp, scale=scale,
        dendrogram=dendrogram, margins=margins,
        key=key, symkey=FALSE, density.info="none", lvtrace="none",
        ColSideColors=ColSideColors, cexRow=cexRow, keysize =keysize, 
        distfun=function(x) dist(x, method=distmethod),
        hclustfun=function(x) hclust(x, method=aggmethod),
        ...
    )
}
