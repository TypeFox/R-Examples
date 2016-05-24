### This files contains functions for ploting multivariate data for clusters.
### Written: Wei-Chen Chen on 2009/02/02.

plotmd <- function(x, class = NULL, xlab = "Variables", ylab = "Data", ...){
  x.a <- 1:ncol(x)
  xlim <- range(x.a)
  ylim <- range(x)

  plot(NULL, NULL, type = "n", axes = FALSE,
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
  if(is.null(class)){
    tmp <- apply(x, 1, function(y) lines(list(x = x.a, y = y)))
  } else{
#    color.my <- c("#00EEEEFF", "#A0EEEEFF", "#AAC6EEFF", "#EEC6BBFF",
#                  "#EE66EEFF", "#DEAA00FF", "#EECF5EFF", "grey65", "grey80",
#                  "#60FF00FF", "#C9FF00FF")
    # color.my <- c("#00EEEEFF", "#AAC6EEFF", "#EE66EEFF", "#EEC6BBFF",
    #               "grey65", "#60FF00FF", "#A0EEEEFF", "#EECF5EFF", "grey80",
    #               "#DEAA00FF", "#C9FF00FF")
    color.my <- color.class
    color <- color.my[class %% length(color.my) + 1]

    tmp <- apply(cbind(color, x), 1,
                 function(y) lines(list(x = x.a, y = y[-1]), col = y[1]))
  }
  box()
  axis(1, at = x.a, labels = x.a)
  axis(2)
}
