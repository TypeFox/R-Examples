FlinnPlot <-
function(oss = 1, lp = 0, out.file, max.k = 5, plot.title = "Flinn diagram", labs){
  pdf(file = out.file, width = 3, height = 3.5, useDingbats = FALSE, family = 'serif')
  plot(1, 1, type = 'n', asp = 1, ylim = c(1, max.k), xlim = c(1, max.k), axes = FALSE, xlab = "Y/Z", ylab = "X/Y", main = plot.title)
  lines(c(1, max.k), c(1, max.k), col = "gray")
  par(xaxs = 'i', yaxs = 'i', ps = 10, mar = c(4, 3, 3, 1))
  for(j in 1:length(oss)){
    temp <- EllipAxes(es = oss[j], nu = lp[j])
x <- temp[1]
y <- temp[2]
z <- temp[3]

y.coord <- x / y
    x.coord <- y / z

points(x.coord, y.coord, pch = 19, cex = .5)

if(!missing(labs)){
  text(x.coord, y.coord, labs[j], cex = .75, col = "gray", adj = c(0,0))
}
  }
  axis(1)
  axis(2)
  dev.off()
}
