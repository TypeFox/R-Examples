# A. Scree device; default settings on the root viewport
library(grid)
win.graph(width = 5.5, height = 5.5); bringToTop(stay = TRUE)
pushViewport(viewport(gp = gpar(fontfamily = 'serif', fontsize = 9)))

# B1. Outside box; two connection curves
grid.roundrect(width = unit(1, "npc") - unit(0.2, "inches"),
  height = unit(1, "npc") - unit(0.2, "inches"), r = unit(0.05, "snpc"))

grid.curve(x1 = 0.06, y1 = 0.84, x2 = 0.88, y2 = 0.52, curvature = -0.95,
  angle = 10, shape = 0.9)

grid.curve(x1 = 0.06, y1 = 0.16, x2 = 0.88, y2 = 0.48, curvature = 0.05,
  angle = 10, shape = 0.9)

# B2. Middle part: author, reader, two arrows, paper
grid.roundrect(x = 0.05, width = unit(0.8, "npc"), just = "left",
  height = unit(0.17, "inches"), r = unit(0.5, "snpc"),
  gp = gpar(col = "white", fill = "gray80"))
grid.roundrect(x = 0.05, y = 0.56, width = unit(0.43, "npc"),
  height = unit(2, "lines"), r = unit(0.4, "snpc"), just = "left",
  gp = gpar(fill = 'white'))
grid.roundrect(x = 0.05, y = 0.44, width = unit(0.43, "npc"),
  height = unit(2, "lines"), r = unit(0.4, "snpc"), just = "left",
  gp = gpar(fill = 'white'))

grid.lines(x = 0.12, y = c(0.61, 0.71), arrow = arrow(),  # two arrows
  gp = gpar(lwd = 10, col = 'gray80', lineend = "round"))
grid.lines(x = 0.12, y = c(0.39, 0.29), arrow = arrow(),
  gp = gpar(lwd = 10, col = 'gray80', lineend = "round"))

grid.roundrect(x = 0.9, y = 0.51, width = 0.11, height = 0.14,  # book
  gp = gpar(fill = 'gray80'), r = unit(0.1, "snpc"))
grid.roundrect(x = 0.89, y = 0.5, width = 0.11, height = 0.14,
  gp = gpar(fill = 'white'), r = unit(0.1, "snpc"))
grid.text(label = c("An author's work", "A reader's memory", "Paper"), 
  x = c(0.15, 0.15, 0.86), y = c(0.56, 0.44, 0.5), just = "left")

# B3. Six boxes and texts
x6 <- c(0.05, 0.30, 0.55, 0.05, 0.30, 0.55)
y6 <- c(0.90, 0.78, 0.66, 0.10, 0.22, 0.34)
for (i in 1:length(x6)) {
  grid.roundrect(x = x6[i], y = y6[i], width = unit(0.32, "npc"),
    height = unit(4, "lines"), r = unit(0.4, "snpc"), just = "left",
    gp = gpar(fill = 'white'))
}

grid.text(x = x6 + 0.05, y = y6, just = "left", label = c(
  "1. Research Idea\n- One sentence\n- A few days or months",
  "2. Outline\n- 2 pages\n- Several months",
  "3. Details\n- 30 pages\n- Several weeks",
  "3. Research Idea\n- Topic only\n- After several years",
  "2. Outline\n- Structure only\n- In a year",
  "1. Details\n- Most facts\n- In several days"))
memo.grid <- recordPlot()  # record the graph

# C. Save a PDF copy
pdf(file = "C:/aErer/fig_gridMemo.pdf", width = 5.5, height = 5.5)
replayPlot(memo.grid)
dev.off()