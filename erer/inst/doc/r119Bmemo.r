# Some constants for curves
a <- 10; b <- 5; k <- 0.1; h <- 0.01; ww <- 5.5; hh <- 5.5

# Graph version A          
win.graph(width = ww, height = hh, pointsize = 9)
par(mai = c(0.4, 0.4, 0.1, 0.1), family = "serif")
curve(expr = k + b * sqrt(1 - ((x - h) / a) ^ 2), 
  from = 0, to = 10, n = 500, xlim = c(0, 10.5), ylim = c(-6.5, 6.5))

# Graph version B
win.graph(width = ww, height = hh, pointsize = 9)
bringToTop(stay = TRUE)
par(mai = c(0.1, 0.1, 0.1, 0.1), family = "serif")
curve(expr = k + b * sqrt(1 - ((x - h) / a)^2), axes = FALSE, ann = FALSE,
  from = 0, to = 10, n = 500, xlim = c(0, 10.5), ylim = c(-6.5, 6.5))

# Shared commands for versions A and B
curve(expr = -k - b * sqrt(1 - ((x - h) / a)^2), 
  from = 0, to = 10, n = 500, add = TRUE) 
box()
segments(x0 = 0, y0 = 0, x1 = 10, y1 = 0, col = "gray80", lwd = 10, lend=2)
arrows(x0 = 1, y0 = c(-1.5, 1.5), x1 = 1, y1 = c(-3,3), col = "gray80",
  lwd = 10, code = 2, length = 0.15, angle = 20, lend = 1, ljoin = 1)
rect(xleft = 0.2, ybottom = c(-0.3, 0.3), xright = 5, ytop = c(-1.2, 1.2),
  col = "white")
  
text(x = 1.1, y = 0.7, labels = "An author's work", pos = 4, adj = c(1, 1))
text(x = 0.9, y = -0.8, labels = "A reader's memory", pos = 4, adj=c(1, 1))
rect(xleft = c(0, 3, 6, 0, 3, 6), 
     ybottom = c(5.2, 3.4, 1.6, -5.2, -3.4, -1.6),  
     xright = c(3.5, 6.5, 9.5, 3.5, 6.5, 9.5), 
     ytop = c(6.6, 4.8, 3.0, -6.6, -4.8, -3.0), col = "white")
rect(xleft = 9.0, ybottom = -1.0, xright = 10.50, ytop = 1.20,
  col = "gray90", border = NULL)
rect(xleft = 8.9, ybottom = -1.1, xright = 10.35, ytop = 1.05, col="white")

text(x = c(0.1, 3.1, 6.1), y = c(5.8, 3.95, 2.15, -6.05, -4.25, -2.4), 
  pos = 4, adj = c(1, 1), labels = c(
  "1. Research Idea\n- One sentence\n- A few days or months",
  "2. Outline\n- 2 pages\n- Several months",
  "3. Details\n- 30 pages\n- Several weeks",
  "3. Research Idea\n- Topic only\n- After several years",
  "2. Outline\n- Structure only\n- In a year",
  "1. Details\n- Most facts\n- In several days"))
text(x = 9.1, y = 0.4, labels = "Paper", pos = 4, adj = c(1, 1))
memo.graph <- recordPlot()  # record the graph for version C

# Graph version C
pdf(file = "C:/aErer/fig_memo.pdf", width = ww, height = hh, pointsize = 9,
  family = "serif")
replayPlot(memo.graph)
dev.off()