grandistrib <-
function (x, main="", scale = "fine", xlab = "Stations", ylab = "Percentage") 
{
  if (scale == "fine") 
    Descript <- granstat(x, aggr = F)$sedim[-c(1:5), ]
  if (scale == "large") 
    Descript <- as.data.frame(granstat(x, aggr = F)$sedim[c(2:5), ])
  Descript = Descript[nrow(Descript):1, ]
  layout(matrix(c(1, 1, 1, 1, 1, 2), nrow = 1), widths = c(0.7, 
                                                           0.3))
  par(mar = c(8, 4, 4, 2) + 0.5)
  barplot(as.matrix(Descript), col = rainbow(nrow(Descript)), 
          font.lab = 2, xlab = xlab, ylab = ylab, main=main, cex.lab = 1.9,cex.names=1.7,cex.axis=1.5,las=2)
  par(mar = c(1, 1, 1, 1) + 0.1)
  plot(1, 1, type = "n", axes = F, xlab = "", ylab = "")
  legend("left", legend = rev(row.names(Descript)), fill = rev(rainbow(nrow(Descript))),
         text.width=6,cex=2,bty="n")
}
