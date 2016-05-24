overallComponentsPlot <- function(comp, ref){
  diffMR <- differenceMR(comp, ref, eval="original")
  old.par <- par(no.readonly = TRUE)
  par(oma = c(0, 0, 0, 14))
  
  barplot(t(diffMR[2:4]), ylab="Difference (percent of domain)")
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(x=0.1, y=0.3, c(colnames(diffMR)[4:2]), bty="n", 
         fill = c(rgb(230,230,230, maxColorValue=255), rgb(174,174,174, maxColorValue=255), 
                  rgb(77,77,77, maxColorValue=255)))
  
  par(old.par)
}