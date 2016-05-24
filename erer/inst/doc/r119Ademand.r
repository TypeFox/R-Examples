# 1. Create a screen device
win.graph(width = 3, height = 4, pointsize = 9); bringToTop(stay = TRUE)

# 2a. Parameter choices for version A
# in.mai <- c(0.5, 0.5, 0.1, 0.1)
# in.axes <- in.ann <- in.draw <- TRUE; in.x <- 1.5; in.y <- 1.5

# 2b. Parameter choices for version B
in.mai <- c(0.3, 0.3, 0.1, 0.1)
in.axes <- in.ann <- in.draw <- FALSE; in.x <- 0.2; in.y <- 0.3

# 3. Shared commands
par(mai = in.mai, mgp = c(2, 0.6, 0), las = 1, family = "serif")
plot(0, xlim = c(0, 10), ylim = c(1, 9.5), type = "n", axes = in.axes,
  ann = in.ann, xaxs = "i", yaxs = "i", xlab = "", ylab = "")

xa <- 4.5; xb <- 5.25; xc <- 4.5; ya <- 5.0; yb <- 5.50; yc <- 6.0
polygon(x = c(0, xa, xb, 0), y = c(ya, ya, yb, yb), border = NA,
  density = 30, lty = "dotted")  # producer surplus
polygon(x = c(0, xb, xc, 0), y = c(yb, yb, yc, yc), border = NA,
  col = "grey90")  # consumer surplus

if (in.draw) {  # adjustment on the axis lines
  box()
} else {
  axis(side = 1, labels = FALSE, lwd.ticks = -1, at = c(0, 10))
  axis(side = 2, labels = FALSE, lwd.ticks = -1, at = c(1, 9.5))
}
lines(x = c(0, 9), y = c(2, 8))  # Supply curve
lines(x = c(0, 9), y = c(8, 2))  # Initial demand curve
lines(x = c(0, 9), y = c(9, 3))  # New demand curve

lines(x = c(0,  xa), y = c(ya, ya), lty = 2)  # five dashed lines
lines(x = c(0,  xb), y = c(yb, yb), lty = 2)
lines(x = c(0,  xc), y = c(yc, yc), lty = 2)
lines(x = c(xc, xc), y = c(0,  yc), lty = 2)
lines(x = c(xb, xb), y = c(0,  yb), lty = 2)

points(x = c(xa, xb), y = c(ya, yb), pch = 16, cex = 1.5)  # 3 points
points(x = xc, y = yc, pch = 21, bg = "grey70")
text(x = c(xa - 0.6, xb, xc), y = c(ya + 0.1, yb, yc),
  labels = c("a", "b", "c"), pos = c(1,3,3))  # labels at plot region
text(x = c(9, 9, 9), y = c(2, 3, 8) + 0.2,
  labels = c(expression(D[a]), expression(D[b]), "S"), pos = 3 )

mtext(text = c(expression(Q[a]), expression(Q[b])),
  side = 1, line = in.x, at = c(xa, xb))  # labels in figure margin
mtext(text = c("d", "e", "f", expression(P[b]), expression(P[a]), "g"),
  side = 2, line = in.y,  at = c(9, 8, yc, yb, ya, 2))
out.demand <- recordPlot()

# Graph version C
pdf(file = "C:/aErer/fig_demandB.pdf", width = 3, height = 4, 
  pointsize = 9, family = "serif")
replayPlot(out.demand)
dev.off()