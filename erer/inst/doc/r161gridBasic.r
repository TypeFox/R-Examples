# A. Three scenarios for "vp" in the plotting function
library(grid)
vp1.descrip <- viewport(width = 0.9, height = 0.9, name = "vp1.pushed")
vp2.descrip <- viewport(width = 0.5, height = 0.5, name = "vp2.pushed")

# A1. vp = NULL (default)
dev.new(); current.viewport()  # ROOT
pushViewport(vp1.descrip); grid.rect(vp = NULL)  # vp1.pushed

# A2. vp = viewport object
dev.new()
grid.text(label = "Output B", gp = gpar(col= 'red'), vp = vp1.descrip)
current.viewport()  # ROOT

dev.new()
pushViewport(vp1.descrip)  # vp1.pushed
grid.text("Output B", gp = gpar(col= 'green'), vp = NULL)
popViewport()       # remove "vp1.pushed"
current.viewport()  # ROOT

# A3. vp = pushed viewport name
dev.new()
pushViewport(vp1.descrip, vp2.descrip); current.vpTree()  # vp2.pushed
upViewport(n = 0)                                         # ROOT
grid.rect(gp = gpar(lty = 'solid'), vp = "vp1.pushed")
downViewport(name = "vp1.pushed")                         # vp1.pushed
grid.rect(gp = gpar(lty = 'dashed'), vp = "vp2.pushed")
upViewport(n = 1)                                         # Root

dev.new()
pushViewport(vp1.descrip, vp2.descrip)  # vp2.pushed
upViewport(0)                           # ROOT
downViewport(name = "vp1.pushed")       # vp1.pushed
grid.rect(gp = gpar(lty = 'solid'))
downViewport(name = "vp2.pushed")       # vp2.pushed
grid.rect(gp = gpar(col = 'purple'))
upViewport(n = 1)                       # vp1.pushed

# B. Plotting functions, parameters, and units
convertX(x = unit(2.54, "cm"), unitTo = "inches")  # 1 inch
unit.c(unit(1:3, "inches"), unit(2:4, "cm"))

# C. Draw a graph with several viewports
windows(width = 5.4, height = 3, pointsize = 9, family = "serif")
vp3.descrip <- viewport(x = 0.25, y = 0.4, width = 0.4, height = 0.4,
  angle = 10, name = "vp3.pushed")
vp4.descrip <- viewport(x = 0.75, y = 0.4, width = 0.4, height = 0.4,
  angle = -10, name = "vp4.pushed")
pushViewport(vpList(vp3.descrip, vp4.descrip))
current.viewport(); current.vpTree()  # a viewport list

upViewport(n = 0)
grid.rect(width = unit(1, "npc") - unit(3, "mm"), 
  height = unit(1, "npc") - unit(3, "mm"))
grid.text(label = "Learning viewports and primitive functions in grid",
  y = unit(1, "npc") - unit(2, "lines"), gp = gpar(fontsize = 9))

seekViewport(name = "vp3.pushed")
grid.rect(gp = gpar(lty = 2))
grid.curve(x1 = 0.1, y1 = 0.1, x2 = 0.9, y2 = 0.9, curvature = 0.4,
  arrow = arrow(), gp = gpar(lwd = 3, col = 'gray'))

seekViewport(name = "vp4.pushed")
grid.roundrect(gp = gpar(lty = 3), r = unit(0.3, "snpc"))
grid.points(x = 0.5, y = 0.5, pch = 11, size = unit(0.8, 'inches'))
out <- recordPlot()

pdf(file = "C:/aErer/fig_gridBasic.pdf", width = 5.4, height = 3)
replayPlot(out); dev.off() 