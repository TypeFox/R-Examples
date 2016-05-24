# A. Window device for a chess board
win.graph(width = 3, height = 3)
bringToTop(stay = TRUE)
par(mai = c(0.1, 0.1, 0.1, 0.1))

# B. Draw cells with different colors
plot(x = 1:9, y = 1:9, type = "n", xaxs = "i", yaxs = "i",
  axes = FALSE, ann = FALSE)
for (a in 1:8) {
  for (b in 1:8) {
    colo <- ifelse(test = (a + b) %% 2 == 0, yes = "gray50", no = "gray98")
    rect(xleft = a, xright = a + 1, ybottom = b, ytop = b + 1,
      col = colo, border = "white")
  }
}
box()

# C. Add chess pieces
points(x = c(2.5, 3.5), y = c(3.5, 6.5), pch = 16, cex = 3)
points(x = c(5.5, 7.5), y = c(4.5, 3.5), pch = 21, cex = 3, bg = "white")
out <- recordPlot()

# D. Save a pdf copy
pdf(file = "C:/aErer/fig_chess.pdf", width = 3, height = 3, 
  useDingbats = FALSE)
replayPlot(out); dev.off()