# A1. Run the program and generate the default graph
setwd("C:/aErer"); ls()
source("r072sunSJAF.r", echo = FALSE); lss()

windows(width = 4, height = 3, pointsize = 9); bringToTop(stay = TRUE)
par(mai = c(0.7, 0.7, 0.1, 0.1), family = "serif"); plot(p1)

pdf(file = "fig_huntingold.pdf", width = 4, height = 3, pointsize = 9)
par(mai = c(0.7, 0.7, 0.1, 0.1), family = "serif"); plot(p1); dev.off()

# A2. Understand the plotting method and data 
names(p1); class(p1); plot.maTrend

# B1. Data preparation 
pr <- p1$trend; class(pr); head(pr)
ya <- seq(from = 0.10, to = 0.45, by = 0.05)
yb <- substr(sprintf("%.2f", ya), start = 2, stop = 4)
mv <- colMeans(p1$q$w$x)[p1$nam.c]

# B2. Graphics device
windows(width = 4, height = 3, pointsize = 9); bringToTop(stay = TRUE)
par(mai = c(0.5, 0.6, 0.1, 0.1), family = "serif", mgp = c(2, 0.6, 0))

# B3. High-level plotting function
plot(pr[, 2] ~ pr[, 1], type = "l", lty = 1, xaxs = "i", yaxs = "i", 
  ylim = range(ya), las = 1, ann = FALSE, axes = FALSE)

# B4. Low-level plotting functions
axis(side = 1); axis(side = 2, las = 1, at = ya, labels = yb)
lines(pr[, 3] ~ pr[, 1], type = "l", lty = 2)
lines(pr[, 4] ~ pr[, 1], type = "l", lty = 3)
abline(v = mv, lty = 4)
title(xlab = "Hunting experience (Year)", 
      ylab = "Prob(Insurance purchase = yes)")  
text(x = c(39, 50, 60), y = c(0.35, 0.22, 0.18), 
  labels = c("Nonresident", "All", "Resident"))
fig.new <- recordPlot()  # save the graph on the screen device

# B5. Save a PDF version
pdf(file = "fig_huntingnew.pdf", width = 4, height = 3, pointsize = 9)
replayPlot(fig.new)
dev.off()