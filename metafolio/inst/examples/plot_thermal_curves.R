# plot thermal curves across populations as an example

#par(cex = 0.8, mgp = c(2.4, 0.65, 0), tck = -0.02)
#optim_temps <- seq(13, 19, length.out = 10)
#x <- seq(3, 29, length.out = 200)
#plot(1, 1, xlim = c(4, 28), ylim = c(-0.01, 1.4), ylab = "Ricker productivity parameter (a)", xlab = expression(Temperature~(degree*C)), type = "n", yaxs = "i", las = 1)
#for(i in 1:10) {
  #lines(x, thermal_curve_a(x, optim_temp = optim_temps[i], max_a = 1.3, width_param = 0.02), ylab = "a", xlab = "Temperature", type = "l", col = col_pal[i], lwd = 1.5)
#}

pdf_eps("thermal-curves", width = 6, height = 4, type = TYPE)
par(cex = 0.8, mgp = c(2.4, 0.65, 0), tck = -0.02, mfrow = c(1,1), mar = c(4,4,.5,.5))
optim_temps <- seq(13, 19, length.out = 10)
widths <- c(seq(0.05, 0.02, length.out = 5), rev(seq(0.05, 0.02, length.out = 5)))
heights <- c(seq(2.8, 2.2, length.out = 5), rev(seq(2.8, 2.2, length.out = 5)))

x <- seq(3, 29, length.out = 200)
#plot(1, 1, xlim = c(4, 28), ylim = c(-0.01, 2.9), ylab = "Ricker productivity parameter (a)", xlab = expression(Temperature~(degree*C)), type = "n", yaxs = "i", las = 1)
plot(1, 1, xlim = c(4, 28), ylim = c(-0.01, 2.9), ylab = "Ricker productivity parameter (a)", xlab = "Environmental value", type = "n", yaxs = "i", las = 1)
for(i in 1:10) {
  a <- thermal_curve_a(x, optim_temp = optim_temps[i], max_a = heights[i], width_param = widths[i])
  a[a<0] <- 0
  lines(x, a, col = col_pal[i], lwd = 1.5)
}

abline(v = 16, lty = 2, lwd = 1.5)
abline(v = c(16-3.5, c = 16+3.5), lty = 3, lwd = 1.5)
dev.off()
