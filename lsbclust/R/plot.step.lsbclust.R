#' Plot method for class 'step.lsbclust'
#' 
#' Plot 'step.lsbclust' objects.
#' 
#' @param x An object of class \code{step.lsbclust}
#' @param which Which type of plot to produce. 
#' @param col.all A character vector of length one indicating which of \code{"overall"}, 
#' \code{"rows"}, \code{"columns"} or \code{"interactions"} should be mapped to colour in the
#' plot for all possible models. Care needs to be taken that the stated component is included in
#' the fit.
#' @param arrange Logical indicating whether the arrange  the plots side-by-side
#' via \code{\link{grid.arrange}} or not.
#' @param chull Logical indicating whether to plot the estimated convex hull or not.
#' @param \dots additional arguments passed to \code{\link{theme}}.
#' @keywords hplot
#' @author Pieter C. Schoonees
#' @method plot step.lsbclust
#' @export
plot.step.lsbclust <- function(x, which = 1L:5L, col.all = NULL, arrange = FALSE, 
                               chull = FALSE, ...){
  
  ## Avoid global variable issues
  DF <- Loss <- Nclust <- Overall <- Rows <- Columns <- Interactions <- NULL
  
  ## Determine components
  nms <- names(x$ind.loss)
  cind <- match(nms, c("overall", "rows", "columns", "interactions"))
  
  ## Process which
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  
  ## Do not plot missing components
  show[-1][!(1:4 %in% cind)] <- FALSE
  
  plots <- vector(mode = "list", length = 5)
  
  ## Plot all losses
  if (show[1]) {
    
    dfAll <- data.frame(Loss = x$loss$total, DF = x$df$total, Nclust = x$nclust)
    
    plots[[1]] <- ggplot(data = dfAll, mapping = aes(x = DF, y = Loss))
          
    ## Possibly map to colour
    plots[[1]] <- plots[[1]] + geom_point(aes_string(colour = col.all))
    
    ## Better labels
    plots[[1]] <- plots[[1]] + ggtitle("All Models") + xlab("Degrees-of-Freedom")
    
    ## Possibly add convex hull
    if (chull) {
      ch <- chull(dfAll$DF, dfAll$Loss)
      plots[[1]] <- plots[[1]] + geom_line(data = dfAll[ch, ], colour = "red")
    }
  }
  
  ## Scree plot for overall means
  if (show[2]) {
    
    dfOvl <- data.frame(Loss = x$ind.loss$overall, DF = x$ind.df$overall, Nclust = x$ind.nclust$overall)
    plots[[2]] <- ggplot(data = dfOvl, mapping = aes(x = Nclust, y = Loss)) + geom_line() + geom_point() +
      ggtitle("Scree plot: Overall means") + xlab("Number of Clusters")
  }

  ## Scree plot for row margins
  if (show[3]) {
    
    dfRows <- data.frame(Loss = x$ind.loss$rows, DF = x$ind.df$rows, Nclust = x$ind.nclust$rows)
    plots[[3]] <- ggplot(data = dfRows, mapping = aes(x = Nclust, y = Loss)) + geom_line() + geom_point() +
      ggtitle("Scree plot: Row margins") + xlab("Number of Clusters")
  }

  ## Scree plot for column margins
  if (show[4]) {
    
    dfColumns <- data.frame(Loss = x$ind.loss$columns, DF = x$ind.df$columns, Nclust = x$ind.nclust$columns)
    plots[[4]] <- ggplot(data = dfColumns, mapping = aes(x = Nclust, y = Loss)) + geom_line() + geom_point() +
      ggtitle("Scree plot: Column margins") + xlab("Number of Clusters")
  }

  ## Scree plot for overall means
  if (show[5]) {
    
    dfInteractions <- data.frame(Loss = x$ind.loss$interactions, DF = x$ind.df$interactions, 
                                 Nclust = x$ind.nclust$interactions)
    plots[[5]] <- ggplot(data = dfInteractions, mapping = aes(x = Nclust, y = Loss)) + geom_line() + geom_point() +
      ggtitle("Scree plot: Interactions") + xlab("Number of Clusters")
  }

  ## Return plots
  if(arrange) do.call(gridExtra::grid.arrange, plots[show])
  else return(plots[show])
}