#' Plot a \code{bicomp} Object
#' 
#' Plot method for an object of class \code{bicomp} (see \code{\link{bicomp}}).
#' 
#' @param x An object of class \code{bicomp}.
#' @param which A numeric vector indicating which matrices to plot, with 0 = original data, 
#' 1 = overall means, 2 = row means, 3 = column means and 4 = interactions.
#' @param arrange Logical indicating whether the arrange  the plots side-by-side
#' via \code{\link{grid.arrange}} or not.
#' @param col A character vector of length three giving the parameters 
#' \code{low}, \code{mid} and \code{high} for \code{\link{scale_fill_gradient2}}.
#' @param strip.legend Logical indicating whether to strip the legend off the plot or not.
#' @param add.titles Logical indicating whether to add titles to the plots or not.
#' @param \dots Additional arguments to \code{\link{theme}}.
#' @importFrom reshape2 melt
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid textGrob gpar
#' @export
plot.bicomp <- function(x, which = 0L:4L, arrange = TRUE, 
                        col = c("red4", "beige", "blue4"), strip.legend = TRUE,
                        add.titles = FALSE, ...) {
  ## Avoid global variable issues
  Row <- Column <- Rating <- Value <- NULL
  
  ## Check x and create list if needed
  if (length(attributes(x)$return_type) == 1) {
    x <- list(x)
    names(x) <- attributes(x[[1]])$return_type
  }
  
  ## Set up show and override arrange and which if necessary
  show <- rep(FALSE, 5)
  which <- which[match(names(x), c("original", "overall", "rows", "columns", "interactions"))]
  show[which + 1] <- TRUE
  if (sum(show) != 5) arrange <- FALSE
  plots <- vector(mode = "list", length = 5L)
#   lims <- range(do.call(cbind, x))
  lims <- c(-1, 1) * max(abs(do.call(cbind, x)))

  ## Only arrange if which is of length 5
  if (length(which) != 5) arrange <- FALSE

  ## Original matrix
  if (show[1L]) {
   dfOrig <- reshape2::melt(x$original, as.is = TRUE)
   colnames(dfOrig) <- c("Row", "Column", "Value")
   dfOrig$Column <- factor(dfOrig$Column, levels = colnames(x$original))
   dfOrig$Row <- factor(dfOrig$Row, levels = rev(rownames(x$original)))
   plots[[1L]] <- ggplot(data = dfOrig, aes(x = Column, y = Row, fill = Value)) + 
      geom_tile(colour = "white") + 
     scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], limits = lims) +
     theme_bw() + theme(...)
   if (add.titles) plots[[1L]] <- plots[[1L]] + ggtitle("Data")
#      ggtitle("Data") + theme(...)
#       ggtitle("Data") + theme(axis.text = element_blank(), axis.ticks = element_blank(), ...)
    
  }

  if (show[2L]) {
    dfOvl <- reshape2::melt(x$overall, as.is = TRUE)
    colnames(dfOvl) <- c("Row", "Column", "Value")
    dfOvl$Column <- factor(dfOvl$Column, levels = colnames(x$overall))
    dfOvl$Row <- factor(dfOvl$Row, levels = rev(rownames(x$overall)))
    plots[[2L]] <- ggplot(data = dfOvl, aes(x = Column, y = Row, fill = Value)) + 
      geom_tile(colour = "white") + 
      scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], limits = lims) +
      theme_bw() +theme(legend.position = "none", ...) 
    if (add.titles) plots[[2L]] <- plots[[2L]] + ggtitle("Overall mean")
#       ggtitle("Overall mean") + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), ...)
#       ggtitle("Overall mean") + theme(...)
  }

  if (show[3L]) {
    dfRows <- reshape2::melt(x$rows, as.is = TRUE)
    colnames(dfRows) <- c("Row", "Column", "Value")
    dfRows$Column <- factor(dfOrig$Column, levels = colnames(x$rows))
    dfRows$Row <- factor(dfOrig$Row, levels = rev(rownames(x$rows)))
    plots[[3L]] <- ggplot(data = dfRows, aes(x = Column, y = Row, fill = Value)) + 
      geom_tile(colour = "white") + 
      scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], limits = lims) +
      theme_bw() +theme(legend.position = "none", ...) 
    if (add.titles) plots[[3L]] <- plots[[3L]] + ggtitle("Row means")
#       ggtitle("Row means") + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), ...)
#       ggtitle("Row means") + theme(...)
  }

  if (show[4L]) {
    dfCols <- reshape2::melt(x$columns, as.is = TRUE)
    colnames(dfCols) <- c("Row", "Column", "Value")
    dfCols$Column <- factor(dfCols$Column, levels = colnames(x$columns))
    dfCols$Row <- factor(dfCols$Row, levels = rev(rownames(x$columns)))
    plots[[4L]] <- ggplot(data = dfCols, aes(x = Column, y = Row, fill = Value)) + 
      geom_tile(colour = "white") + 
      scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], limits = lims) +
      theme_bw() + theme(legend.position = "none", ...) 
    if (add.titles) plots[[4L]] <- plots[[4L]] + ggtitle("Column means")
#       ggtitle("Column means") + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), ...)
#       ggtitle("Column means") + theme(...)
  }

  if (show[5L]) {
    dfInt <- reshape2::melt(x$interactions, as.is = TRUE)
    colnames(dfInt) <- c("Row", "Column", "Value")
    dfInt$Column <- factor(dfInt$Column, levels = colnames(x$interactions))
    dfInt$Row <- factor(dfInt$Row, levels = rev(rownames(x$interactions)))
    plots[[5L]] <- ggplot(data = dfInt, aes(x = Column, y = Row, fill = Value)) + 
      geom_tile(colour = "white") + 
#       scale_fill_gradient(low = "red", high = "blue", limits = lims)
      scale_fill_gradient2(low = col[1], mid = col[2], high = col[3], limits = lims) +
      theme_bw() + theme(legend.position = "none", ...) 
  if (add.titles) plots[[5L]] <- plots[[5L]] + ggtitle("Interactions")
#       ggtitle("Interactions") + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), ...)
#       ggtitle("Interactions") + theme(...)
  }

  if (arrange) {
#     do.call(gridExtra::grid.arrange, c(plots[which], nrow = 1))
    grid.arrange(gridExtra::arrangeGrob(plots[[1]], grid::textGrob("=", gp = grid::gpar(fontsize=50)), plots[[2]], 
                                        grid::textGrob("+", gp = grid::gpar(fontsize=50)), 
                                        nrow = 1, widths = c(0.365, 0.08, 0.28, 0.08), default.units = "npc"), 
                 gridExtra::arrangeGrob(plots[[3]], 
                                        grid::textGrob("+", gp = grid::gpar(fontsize=50)), plots[[4]],
                                        grid::textGrob("+", gp = grid::gpar(fontsize=50)), plots[[5]], 
                                        widths = c(0.28, 0.08, 0.28, 0.08, 0.28), nrow = 1, default.units = "npc"), ncol = 1)
  }
  else return(plots[which + 1])
}
