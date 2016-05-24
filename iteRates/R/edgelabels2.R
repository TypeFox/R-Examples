edgelabels2 <- function (text, edge, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
    thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
    bg = "lightgreen", horiz = FALSE, width = NULL, height = NULL, 
    date = NULL, ...) 
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(edge)) {
        sel <- 1:dim(lastPP$edge)[1]
        subedge <- lastPP$edge
    }
    else {
        sel <- edge
        subedge <- lastPP$edge[sel, , drop = FALSE]
    }
    if (lastPP$type == "phylogram") {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[, 
                2]])/2
            YY <- lastPP$yy[subedge[, 2]]
        }
        else {
            XX <- lastPP$xx[subedge[, 2]]
            YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[, 
                2]])/2
        }
    }
    else {
        XX <- (lastPP$xx[subedge[, 1]] + lastPP$xx[subedge[, 
            2]])/2
        YY <- (lastPP$yy[subedge[, 1]] + lastPP$yy[subedge[, 
            2]])/2
    }
    if (!is.null(date)) 
        XX[] <- max(lastPP$xx) - date
    BOTHlabels2(text, sel, XX, YY, adj, frame, pch, thermo, pie, 
        piecol, col, bg, horiz, width, height, ...)
}
