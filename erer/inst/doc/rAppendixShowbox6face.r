# Load library and data
setwd("C:/aErer"); library(aplpack)
data(longley); str(longley)
longley[as.character(c(1947, 1952, 1957, 1962)), 1:4]

# Some practices
windows(); bringToTop(stay = TRUE)
faces()
faces(face.type = 0)
faces(xy = rbind(1:4, 5:3, 3:5, 5:7), face.type = 2)
faces(xy = longley[c(1, 6, 11, 16), ], face.type = 1)
faces(xy = longley[c(1, 6, 11, 16), ], face.type = 0)

# Display on the screen device
windows(width = 5.3, height = 2.5, pointsize = 9); bringToTop(stay = TRUE)
par(mar = c(0, 0, 0, 0), family = "serif")
aa <- faces(xy = longley[c(1, 6, 11, 16), ], plot = FALSE)
class(aa)  # "faces"
plot.faces(x = aa, face.type = 2, width = 1.1, height = 0.9)
showFace <- recordPlot()

# Save the graph on a file device
pdf(file = "fig_showFace.pdf", width = 5.3, height = 2.5, pointsize = 9)
replayPlot(showFace); dev.off()