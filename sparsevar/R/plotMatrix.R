#' @title Matrix plot
#' 
#' @description Plot a sparse matrix
#' 
#' @param M the matrix to plot
#' @return An \code{image} plot with a particular color palette (black zero entries, red 
#' for the negative ones and green for the positive)
#' @usage plotMatrix(M)
#' 
#' @export
plotMatrix <- function(M) {
  
  nr <- nrow(M)
  nc <- ncol(M)
  M <- t(M)[, nc:1]
  ggplot2::ggplot(reshape2::melt(M), ggplot2::aes_string(x='Var1', y='Var2', fill='value')) + ggplot2::geom_raster() + 
           ggplot2::scale_fill_gradient2(low='red', high='green', mid='black') + ggplot2::xlab("Row") + ggplot2::ylab("Col")
  
}

#' @title VAR plot
#' 
#' @description Plot all the matrices of a VAR model
#' 
#' @param A the list containing the VAR matrices to be plotted
#' @return An \code{image} plot with a particular color palette (black zero entries, red 
#' for the negative ones and green for the positive)
#' @usage plotVAR(A)
#' 
#' @export
plotVAR <- function(A) {
  
  p <- length(A)
  
  pl <- list()
  # par(mfrow = c(1,p))
  for (i in 1:p) {
    pl[[i]] <- plotMatrix(A[[i]])
  }
  
  multiplot(plotlist = pl, cols = p)
}

#' @title Plot VAR models for comparison
#' 
#' @description Plot all the matrices of a two VAR models
#' 
#' @param var1 the list containing the first VAR model matrices to be plotted
#' @param var2 the list containing the second VAR model matrices to be plotted
#' @return An \code{image} plot with a particular color palette (black zero entries, red 
#' for the negative ones and green for the positive)
#' @usage plotComparisonVAR(var1, var2)
#' 
#' @export
plotComparisonVAR <- function(var1, var2) {
  
  p <- length(var2$A)
  
  pl <- list()
  # par(mfrow = c(2,p))
  for (i in 1:p) {
    pl[[i]] <- plotMatrix(var1$A[[i]])
  }
  #title("Estimate")
  for (i in 1:p) {
    pl[[i + p]] <- plotMatrix(var2$A[[i]])
  }
  #title("Simulation")
  multiplot(plotlist = pl, cols = p, layout = matrix(1:(2*p), nrow = 2, byrow = TRUE))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# From R Cookbook
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}