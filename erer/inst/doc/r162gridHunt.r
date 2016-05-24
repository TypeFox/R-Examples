# A. Run the program and generate the data
setwd("C:/aErer"); source("r072sunSJAF.r", echo = FALSE); library(grid)
pr <- p1$trend; class(pr); tail(pr); plot(p1)
xa <- seq(from = 0, to = max(pr[, 1]), by = 10)
ya <- seq(from = 0.10, to = max(pr[, 2:4]), by = 0.05)
yb <- substr(sprintf("%.2f", ya), start = 2, stop = 4)
mv <- colMeans(p1$q$w$x)[p1$nam.c]; u <- 'native'

# B. Draw the graph by grid
# B1. Screen device
windows(width = 4, height = 3); bringToTop(stay = TRUE)

# B2. Create plot region and figure margins
vp.plot <- plotViewport(margins = c(2.6, 3, 1, 1), name = "reg.plot",
  gp = gpar(fontfamily = 'serif', fontsize = 9))
vp.data <- dataViewport(xscale = c(0, max(pr$HuntYrs)),
  yscale = c(0.10, max(pr[, 2:4])), name = "reg.data")
vp.lege <- viewport(x = unit(5, u), y = unit(0.4, u), 
  width = stringWidth("Nonresident"), height = unit(3, "lines"), 
  just = c("left", "top"), name = "reg.lege") 
pushViewport(vpStack(vp.plot, vp.data, vp.lege))
current.viewport(); current.vpTree()

# B3. Data viewport: four lines, axes, texts
seekViewport(name = "reg.data")
grid.xaxis(at = xa, label = as.character(xa));  # grid.rect() 
grid.yaxis(at = ya, label = yb)
grid.lines(x = unit(pr$HuntYrs, u), y = unit(pr[, 2], u), gp = gpar(lty=1))
grid.lines(x = unit(pr$HuntYrs, u), y = unit(pr[, 3], u), gp = gpar(lty=2))
grid.lines(x = unit(pr$HuntYrs, u), y = unit(pr[, 4], u), gp = gpar(lty=3))
grid.lines(x = unit(mv, u), y = unit(c(0, 1), "npc"), gp = gpar(lty = 4))

grid.text("Hunting experience (Year)", y = unit(-2.5, "lines"))
grid.text("Prob(Insurance purchase = yes)", x = unit(-3, "lines"), rot=90)
grid.text(label = c("Nonresident", "All", "Resident"), 
  x = unit(c(39, 50, 60), u), y = unit(c(0.35, 0.22, 0.18), u))

# B4. Legend viewport: add legend
downViewport(name = "reg.lege")
grid.legend(labels = c("Nonresident", "All", "Resident"), ncol = 1,
  gp = gpar(lty = c(2, 1, 3)), vgap = unit(0, "lines"))
out <- recordPlot()  # save the display list

# C. Save PDF on a file device
pdf(file = "fig_gridHunt.pdf", width = 4, height = 3, family = "serif")
replayPlot(out); dev.off() 