# define a baseline etopo colour palette
etopo <- read.csv(textConnection(
"altitudes,colours
10000,#FBFBFB
4000,#864747
3900,#7E4B11
2000,#9B8411
1900,#BD8D15
300,#F0CF9F
0,#307424
-1,#AFDCF4
-12000,#090B6A
"
), stringsAsFactors=FALSE)
etopo$altitudes01 <- scales::rescale(etopo$altitudes)


etopo.colors <- function(n) {
  colorRampPalette(etopo$colours)(n)
}

scale_fill_etopo <- function(...) {
  ggplot2::scale_fill_gradientn(colours=etopo$colours, values=etopo$altitudes01, limits=range(etopo$altitudes), ...)
}

scale_colour_etopo <- function(...) {
  ggplot2::scale_colour_gradientn(colours=etopo$colours, values=etopo$altitudes01, limits=range(etopo$altitudes), ...)
}
