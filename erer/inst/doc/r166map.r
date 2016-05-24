# A. Data preparation
library(maps); options(stringsAsFactors = FALSE)
info <- map(database = 'state', plot = FALSE)
str(info); info[['names']][c(1:3, 38:40, 53:55)]
library(mapdata); map(database = 'china', col = "red")

mbr2 <- data.frame(region = c('florida', 'texas', 'alabama', 
    'north carolina:knotts', 'north carolina:main', 'north carolina:spit',
    'virginia:chesapeake', 'virginia:chincoteague', 'virginia:main', 
    'mississippi', 'georgia', 'louisiana', 'south carolina', 
    'arkansas', 'tennessee', 'oklahoma', 'kentucky'), 
  color = c(rep('grey50', times = 10), rep('grey75', times = 3), 
    rep('grey98', times = 4)))
mbr <- mbr2[order(mbr2$region), ]; rownames(mbr) <- 1:nrow(mbr); tail(mbr)

center <- data.frame(
  abb = c('AL', 'AR', 'FL', 'GA', 'KY', 'LA', 'MS', 'NC', 'OK', 'SC', 
    'TN', 'TX', 'VA', 'Gulf of Mexico', 'Atlantic Ocean'),
  long = c(-86.7, -92.1, -81.5, -83.2, -85.7, -92.0, -89.9, -78.5,
           -97.0, -81.0, -86.0, -99.0, -78.5, -89.0, -78.5),
  lati = c(32.6, 34.8, 28.1, 32.7, 37.3, 30.5, 32.6, 35.2,
           35.3, 33.6, 35.8, 31.2, 37.2, 27.5, 31.0))
tail(center)

# B1. Screen device - 13 states within USA
win.graph(width = 5.5, height = 3.4); bringToTop(stay = TRUE)
map(database = 'state', mar = c(0, 0, 0, 0))
map(database = 'state', region = mbr$region, fill = TRUE, col = 'black',
  add = TRUE)
map(database = 'state', region = mbr$region, fill = FALSE, col = 'white',
  add = TRUE)  
out1 <- recordPlot()

# B2. Screen device - regulations in 13 states
win.graph(width = 5.5, height = 3.4, pointsize = 8); bringToTop(stay=TRUE)
par(mgp = c(1, 0.5, 0), tcl = -0.3, family = "serif")
map(database = 'state', region = mbr$region, fill = TRUE, col = mbr$color,
  xlim = c(-107, -75), ylim = c(24.5, 40), mar = c(4, 4, 0.2, 0.2))
axis(side = 1, at = c(-104, -98, -92, -86, -80))
axis(side = 2, at = c(27, 30, 33, 36, 39), las = 1); box()
# title(xlab = 'Longitude', ylab = 'Latitude'); locator(n = 1)
text(x = center$long, y = center$lati, labels = center$abb)
legend(x = -107, y = 39.8, legend = c('I', 'II', 'III'), 
   box.col = 'white', fill = sort(unique(mbr$color)))
out2 <- recordPlot()

# C. Save a pdf copy
pdf(file = 'fig_map1.pdf', width = 5.5, height = 3.4)  # ratio w/h = 1.6
replayPlot(out1); dev.off()
pdf(file = 'fig_map2.pdf', width = 5.5, height = 3.4)
replayPlot(out2); dev.off()