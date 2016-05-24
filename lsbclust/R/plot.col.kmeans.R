#' Plot method for class 'col.kmeans'
#' 
#' Simple plot method for object of class 'col.kmeans' as output by \code{\link{orc.lsbclust}}.
#' 
#' @param x An object of class \code{col.kmeans}
#' @param which Which type of plot to produce (only 3 types are implemented).
#' @param \dots additional arguments passed to \code{\link{theme}}.
#' @keywords hplot
#' @author Pieter C. Schoonees
#' @method plot col.kmeans
#' @export
#' @examples
#' data("dcars")
#' m <- orc.lsbclust(data = dcars, margin = 3, delta = c(1,1,1,1), nclust = 5, type = "columns")
#' plot(m)
plot.col.kmeans <- function(x, which = 1L, ...){
  
  ## Avoid global variable issues
  Cluster <- Column <- Value <- NULL
  
  ## Set up show
  show <- rep(FALSE, 3)
  show[which] <- TRUE
  
  ## Determine lengths etc.
  k <- length(x$size)
  N <- sum(x$size)
  K <- ncol(x$centers)
  levs <- paste0("C", 1:k, " (", round(100 * x$size / N, 1), "%)")
#   levs <- paste0("Segment ", 1:k, " (n = ", x$size, ")")
  df <- data.frame(Value = as.numeric(t(x$centers)), 
                   Cluster = factor(rep(levs, each = K), levels = levs), 
                   Column = factor(rep(colnames(x$centers), k), levels = rev(colnames(x$centers))))
  if (show[1L]){
    plot1 <- ggplot(data = df, mapping = aes(x = Value, y = Column, group = Cluster, fill = Cluster)) + 
      geom_vline(colour = "grey", linetype = 2, xintercept = 0) + 
      geom_segment(mapping = aes(yend = Column), xend = 0, size = 0.5) + 
      geom_point(size = 3, shape = 21) + facet_wrap(~ Cluster) + theme_bw(base_size = 16) + 
      theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
            panel.grid.major.y = element_line(linetype = 2), legend.position = "none") + 
      theme(...)
  }
  if (show[2L]) {
    plot2 <- ggplot(data = df) + 
      geom_line(mapping = aes(x = Column, y = Value, colour = Cluster, group = Cluster), size = 1.5) + 
      theme(legend.title = element_blank())  + ylim(c(min(0, min(df$Value)), max(df$Value)))  + 
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(...)
  }
  if (show[3L]) {
    plot3 <- ggplot(data = df, mapping = aes(x = Cluster, y = Column, fill = Value)) + 
      geom_tile() + scale_fill_gradient2(low = "red4", mid = "beige", high = "blue4") +
      theme_bw(base_size = 16)
  }
  mget(paste0("plot", which))
}
