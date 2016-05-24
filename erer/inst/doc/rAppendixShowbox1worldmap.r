# A. Load the package and understand region names
setwd("C:/aErer"); library(maps)
dat <- map(database = "world", plot = FALSE); str(dat)

# B. Display the map on the screen device
windows(width = 5.3, height = 2.5); bringToTop(stay = TRUE)
map(database = "world", fill = FALSE, col = "green", mar = c(0, 0, 0, 0))
map(database = "world", regions = c("Brazil", "China"), fill = TRUE,
  col = c("yellow", "red"), add = TRUE)
showWorld <- recordPlot()

# C. Save the map on a file device
pdf(file = "fig_showWorld.pdf", width = 5.3, height = 2.5)
replayPlot(showWorld); dev.off()