bubbleLegend <-
function (shape="t",
                        attribute,
                        colPalette=NULL,
                        legendName="Legend",
                        bgc='white',
                        scale.level=1,
                        strokeColor="#FFAA00",
                        temp=FALSE,
                        dirname="" ) {
  if (strokeColor == "") {
    strokeColor = "black"
  }
  png(filename = ifelse(temp, paste(tempdir(), "/", legendName, 
                                    ".png", sep = ""), paste(dirname,"/", legendName, ".png", sep = "")), 
      width = 150, height = nlevels(factor(attribute)) * 40, 
      units = "px", bg = "white", pointsize = 10)
  
  par(mar = c(0, 0, 0, 0), bg = bgc)
  plot(0, xlab = "", ylab = "", type = "n", axes = F, asp = 1, 
       xlim = c(0, 15), ylim = c(0, 3 * nlevels(factor(attribute))))
  pal <- colorRampPalette(c("green", "orange", "brown"), space = "Lab")
  niv <- levels(factor(attribute))
  if (is.null(colPalette)) {
    cols <- pal(length(niv))
  }
  else {
    cols <- colPalette
  }
  k = 1
  if (shape == "t") {
    for (i in nlevels(factor(attribute)):1) {
      x <- 1.5 + cos(seq(-pi/6, 2 * pi, by = pi * 2/3)) * 
        1.5 * scale.level[k]
      y <- i * 3 - 1.5 + sin(seq(-pi/6, 2 * pi, by = pi * 
                                   2/3)) * 1.5 * scale.level[k]
      polygon(x, y, col = cols[k], border = strokeColor)
      text(10, i * 3 - 1.5, niv[k], cex = 1)
      k = k + 1
    }
    graph1 <- dev.cur()
    dev.off(graph1)
  }
  else if (shape == "q") {
    for (i in nlevels(factor(attribute)):1) {
      x <- 1.5 + cos(seq(-pi/4, 2 * pi, by = pi/2)) * 1.5 * 
        scale.level[k]
      y <- i * 3 - 1.5 + sin(seq(-pi/4, 2 * pi, by = pi/2)) * 
        1.5 * scale.level[k]
      polygon(x, y, col = cols[k], border = strokeColor)
      text(10, i * 3 - 1.5, niv[k], cex = 1)
      k = k + 1
    }
    graph1 <- dev.cur()
    dev.off(graph1)
  }
  else {
    for (i in nlevels(factor(attribute)):1) {
      x <- 1.5 + cos(c(seq(0, 2 * pi, 0.2), 0)) * 1.5 * 
        scale.level[k]
      y <- i * 3 - 1.5 + sin(c(seq(0, 2 * pi, 0.2), 0)) * 
        1.5 * scale.level[k]
      polygon(x, y, col = cols[k], border = strokeColor)
      text(10, i * 3 - 1.5, niv[k], cex = 1)
      k = k + 1
    }
    graph1 <- dev.cur()
    dev.off(graph1)
  }
}