# (c) Kevin Dunn, 2015.

# require(png)

tradeOffTable <- function(){
  plot.new()
  img <- readPNG(system.file("trade-off-table.png", package="pid"))
  grid::grid.raster(img)
}