draw.oriloc <- function(ori, main = "Title",
  xlab = "Map position in Kb",
  ylab = "Cumulated combined skew in Kb", las = 1, las.right = 3,
  ta.mtext = "Cumul. T-A skew", ta.col = "pink", ta.lwd = 1,
  cg.mtext = "Cumul. C-G skew", cg.col = "lightblue", cg.lwd = 1,
  cds.mtext = "Cumul. CDS skew", cds.col = "lightgreen", cds.lwd = 1,
  sk.col = "black", sk.lwd = 2,
  add.grid = TRUE, ...){

#
# Get data and use Kb units for skews and map positions:
#
  ta <- ori$x/1000          # Cumulated T-A skew in Kb
  cg <- ori$y/1000          # Cumulated C-G skew in Kb
  skew <- ori$skew/1000     # Cumulated combined skew in Kb
  cdsskew <- ori$CDS.excess # Cumulated CDS orientation bias
  start.kb <- ori$start.kb
  end.kb <- ori$end.kb
#
# This is to deal with CDS that wrapp around the genome:
#
  wrapped <- which(abs(ori$start.kb - ori$end.kb) >= 50)

  if(length(wrapped)!=0){
    if(wrapped == length(cdsskew)){
      start.kb[wrapped] <- max(start.kb[wrapped], end.kb[wrapped])
      end.kb[wrapped] <- max(start.kb[wrapped], end.kb[wrapped])
    }
  if(wrapped == 1){
    end.kb[wrapped] <- min(start.kb[wrapped], end.kb[wrapped])
    start.kb[wrapped] <- min(start.kb[wrapped], end.kb[wrapped])
  }
  }
#
# Use CDS midpoints as x-coordinates:
#
  meancoord <- (start.kb + end.kb)/2

  ymin <- min(ta, cg, skew)
  ymax <- max(ta, cg, skew)
  xmin <- min(meancoord)
  xmax <- max(meancoord)
  
  ticks <- pretty(cdsskew)
  
  ticks.y <- (ymax-ymin)/(max(cdsskew) - min(cdsskew))*(ticks - min(cdsskew)) + ymin
  cds.y   <- (ymax-ymin)/(max(cdsskew) - min(cdsskew))*(cdsskew - min(cdsskew)) + ymin

  plot(meancoord, cg, type="l", xlab = xlab,
      ylab = ylab, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
      cex.lab = 1.35, col = cg.col, main = main, lwd = cg.lwd, las = las, ...)

  axis(side = 4, at = ticks.y, labels = ticks, col = cds.col, las = las.right)

#
# Add vertical grid when required:
#
  if(add.grid){
    tmp <- pretty(meancoord)
    abline(v = tmp, col = "grey", lty=3)
    tmp <- tmp[-length(tmp)] + diff(tmp)/2
    abline(v = tmp, col = "grey", lty=3)
  }
  
  lines(meancoord, ta, col = ta.col, lwd = ta.lwd)
  lines(meancoord, skew, col = sk.col, lwd = sk.lwd)
  lines(meancoord, cds.y, col = cds.col, lwd = cds.lwd)
  
  mtext(ta.mtext, col = ta.col, adj = 0)
  mtext(cg.mtext, col = cg.col)
  mtext(cds.mtext, col = cds.col, adj=1)

}
