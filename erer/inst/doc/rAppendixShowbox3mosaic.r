# A. Understand data; library(vcd) has more options for mosaic plots
setwd("C:/aErer")
Titanic; as.data.frame(Titanic)
str(Titanic)
ftable(Titanic)

# B. Draw a mosaic plot for two categorical variables
windows(width = 5.3, height = 2.5, pointsize = 9)
bringToTop(stay = TRUE)
par(mai = c(0.4, 0.4, 0.1, 0))
mosaicplot(formula = ~ Class + Survived, data = Titanic, 
  color = c("red", "green"), main = "", cex.axis = 1)
showMosaic <- recordPlot()

# C. Save the graph on a file device
pdf(file = "fig_showMosaic.pdf", width = 5.3, height = 2.5)
replayPlot(showMosaic)
dev.off()