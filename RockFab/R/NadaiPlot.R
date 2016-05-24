NadaiPlot <-
function(oss = 1, lp = 0, out.file, oss.int = 1, max.oss = 3, plot.title = "Nadai plot", labs){
  #Coerce max.oss to integer
  max.oss <- as.integer(max.oss)
  #Error for a value of zero
  if(max.oss <= 0){stop("Variable: 'max.oss' too small!")}
  
  #Set up PDF device
  pdf(file = out.file, width = 3, height = 3.5, useDingbats = FALSE, family = 'serif')
  par(mai = c(0, 0, .25, 0), omi = c(0, 0, .25, 0))
  plot(0, 0, type = 'n', asp = 1, ylim = c(0, max.oss + (max.oss / 32)), axes = FALSE, ann = FALSE)
  
  #Loop through integers from 0 to max.oss and add curves
  for(j in seq(from = 0, to = max.oss, by = oss.int)){
    theta <- seq(from = (pi / 3), to = (2 * pi / 3), length = 50)
x <- j * cos(theta)
y <- j * sin(theta)
text(x[1], y[1], j, adj = c(0, 1))
if(j < max.oss){
  lines(x, y, col = "gray", lwd = .5)
}
  }
  
  #Add iso-lode lines
  lines(c(0, max.oss * cos((5 * pi) / 12)), c(0, max.oss * sin((5 * pi) / 12)), col = "gray", lwd = .5)
  lines(c(0, max.oss * cos((7 * pi) / 12)), c(0, max.oss * sin((7 * pi) / 12)), col = "gray", lwd = .5)
  lines(c(0,0), c(0, max.oss))
  lines(c(0, max.oss * cos(pi / 3)), c(0, max.oss * sin(pi / 3)), lwd = 1)
  lines(c(0, max.oss * cos(2 * pi / 3)), c(0, max.oss * sin((4 * pi) / 6)), lwd = 1)
  lines(x, y, lwd = 1)
  
  #Add standard text
  text(0, max.oss + (max.oss / 32), expression(nu), adj = c(.5, 0), cex = 1.5)
  text(max.oss * cos(2 * pi / 3) / 2, max.oss * sin( 2 * pi / 3) / 2, expression(bar(epsilon)[s]), adj = c(1, 1), cex = 1.5)
  text((max.oss + (max.oss / 32)) * cos((5 * pi) / 12), (max.oss + (max.oss / 32)) * sin((5 * pi) / 12), "Oblate", adj = c(.5, 0), srt = 345)
  text((max.oss + (max.oss / 32)) * cos((7 * pi) / 12), (max.oss + (max.oss / 32)) * sin((7 * pi) / 12), "Prolate", adj = c(.5, 0), srt = 15)  
  text((max.oss + (max.oss / 32)) * cos(pi / 3), (max.oss + (max.oss / 32)) * sin(pi / 3), "+1", adj = c(1, 0), srt = 330)
  text((max.oss + (max.oss / 32)) * cos(2 * pi / 3), (max.oss + (max.oss / 32)) * sin( 2 * pi / 3), "-1", adj = c(0, 0), , srt = 30)
  mtext(text = plot.title, side = 3, outer = TRUE, cex = 1.25)
  for(j in 1:length(oss)){
    x.coord <- oss[j] * sin((pi / 6) * lp[j])
    y.coord <- oss[j] * cos((pi / 6) * lp[j])
    points(x.coord, y.coord, pch = 19, cex = .5)
if(!missing(labs)){
  text(x.coord, y.coord, labs[j], cex = .75, col = "gray", adj = c(0,0))
}
  }
  dev.off()
}
