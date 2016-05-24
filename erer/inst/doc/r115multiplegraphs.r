# A. Multiple series or columns by plot()
library(erer); data(daBed, daBedRaw)
toto <- window(daBedRaw[, "vWD"], start = c(2001, 1), end = c(2008, 12))
tot <- toto / 10 ^ 6
sha <- daBed[, c('sCN', 'sVN', 'sID', 'sMY', 'sCA', 'sBR', 'sIT')] * 100
bed <- ts.union(tot, sha)
colnames(bed) <- c('All value ($ mil.)', 'China (%)', 'Vietnam (%)',
  'Indonesia (%)', 'Malaysia (%)', 'Canada (%)', 'Brazil (%)', 'Italy (%)')
head(bed)

windows(width = 5.5, height = 5, family = 'serif', pointsize = 9)
plot(x = bed, xlab = "", main = "", oma.multi = c(2.5, 0, 0.2, 0))

pdf(file = "C:/aErer/fig_mtimeseries.pdf", width = 5.5, height = 5,
  family = "serif", pointsize = 9)
plot(x = bed, xlab = "", main = "", oma.multi = c(2.5, 0, 0.2, 0))
dev.off()

# pair of data frame objects
dnor <- data.frame(aa = rnorm(30), bb = rnorm(30), cc = rnorm(30))
windows(); bringToTop(stay = TRUE); plot(dnor) 

# B. split.screen()
windows(); bringToTop(stay = TRUE)
split.screen(figs = c(3, 2))  # split display into 6 screens
for (i in 1:6) {
  screen(n = i)  # open screen
  plot(rnorm(10), pch = as.character(i))
}
close.screen(all.screens = TRUE)  # close the split status

# C. Parameter setting by par(): mfcol, mfrow
windows(); bringToTop(stay = TRUE); par(mfrow = c(3, 2)) 
for (i in 1:6) {
  plot(rnorm(10), pch = as.character(i))
}

# D. layout()
windows(); bringToTop(stay = TRUE)
layout(mat = matrix(1:6, byrow = FALSE, ncol = 2)); layout.show(n = 6)
for (i in 1:6) {
  plot(rnorm(10), pch = as.character(i))
}

# E. use "plt" and "fig" parameter for flexible arrangements
# An empty plot  + first graph + overlapping second graph
win.graph(); bringToTop(stay = TRUE)
plot(x = 0:1, y = 0:1, type = "n", xlab = "", ylab = "", axes = FALSE)
box(which = "figure", lty = "dashed", col = "black", lwd = 5)
box(which = "plot", lty = "dashed", col = "blue", lwd = 5)
fp1 <- par(c("fig", "plt"))

par(new = TRUE, plt = c(x1 = 0.15, x2 = 0.6, y1 = 0.4, y2 = 0.85))
plot(x = rnorm(1000), col = "red", xlab = "", ylab = "")
box(which = "figure", lty = "dashed", col = "yellow", lwd = 5)
box(which = "plot", lty = "dashed", col = "orange", lwd = 5)
fp2 <- par(c("fig", "plt"))

par(new = TRUE, fig = c(x1 = 0.4, x2 = 0.98, y1 = 0.02, y2 = 0.6), 
                plt = c(x1 = 0.2, x2 = 0.9, y1 = 0.3,  y2 = 0.8))
plot(x = rnorm(1000), col = "blue", , xlab = "", ylab = "")
box(which = "figure", lty = "dashed", col = "yellow", lwd = 5)
box(which = "plot", lty = "dashed", col = "orange", lwd = 5)
fp3 <- par(c("fig", "plt"))

# show the figure and plot region locations
fp <- round(cbind(unlist(fp1), unlist(fp2), unlist(fp3)), digits = 2)
rownames(fp) <- c("fig.x1", "fig.x2", "fig.y1", "fig.y2", 
                  "plt.x1", "plt.x2", "plt.y1", "plt.y2")
colnames(fp) <- c("plot.A", "plot.B", "plot.C"); fp