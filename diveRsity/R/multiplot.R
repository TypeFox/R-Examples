multiplot <- function(..., plotlist=NULL, cols) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  plotCols = cols
  plotRows = ceiling(numPlots/plotCols)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
}