plot.R2s <- function(x, ...) {
  if ((length(x$Comp) == 1)) {
    stop("You don't need a graph for a single component model; just use 'print'")
  }
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
  R2.values <- data.frame(Comp = x[[1]], CVR2 = x[[2]],  R2X = x[[3]],  R2Y = x[[4]])
  R2.values.m <- melt(R2.values, id = 1, value.name = "R2")
  R2.values.m$variable <- factor(R2.values.m$variable, levels = levels(R2.values.m$variable), 
                               labels = c("CVR2", "R2Y", "R2X"))
  X_R2.values <- data.frame(Comp = R2.values$Comp, R2.values$R2X)
  Y_R2.values <- data.frame(Comp = R2.values$Comp, CVR2 = R2.values$CVR2, 
                          R2.values$R2Y)
  Y_R2.values.m <- melt(Y_R2.values, id = 1, value.name = "R2")
  Y_R2.values.m$variable <- factor(Y_R2.values.m$variable, levels = levels(Y_R2.values.m$variable), 
                                 labels = c("CVR2", "R2Y"))
  XR2.values.m <- melt(X_R2.values, id = 1, value.name = "R2")
  XR2.values.m$variable <- factor(XR2.values.m$variable, levels = levels(XR2.values.m$variable), 
                                labels = c("R2X"))
  p1 <- with(Y_R2.values.m, ggplot(Y_R2.values.m, aes_string(x = Comp, y = R2, colour = variable, group = variable)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    geom_point(aes_string(col = variable, pch = variable), size = 5) + 
    geom_line(aes_string(col = variable)) + 
    xlab("Latent Variable Number") + 
    ylab("R2") + 
    ggtitle("R2 Diagnostics") + 
    scale_colour_manual(name = "Metric", values = c("red", "blue")) + 
    scale_shape_manual(name = "Metric", values = c(15, 17)))
  p2 <- with(XR2.values.m, ggplot(XR2.values.m, aes_string(x=Comp, y=R2, colour=variable, group=variable)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
    geom_point(aes_string(col = variable), size = 5,colour="black") + 
    geom_line(aes_string(col = variable),colour="black") + 
    xlab("Latent Variable Number") + 
    ylab("R2X") + 
    ggtitle("R2X Diagnostics"))
multiplot(p1, p2, cols=2)
}








