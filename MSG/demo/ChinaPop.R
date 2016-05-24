library(KernSmooth)
x = ChinaPop
x[, 1:2] = apply(x[, 1:2], 2, function(z) 20 * (z -
    min(z))/(max(z) - min(z)) + 5)
symbols(x[, 4], x[, 5], thermometers = x[, 1:3], fg = "gray40",
    inches = 0.5, xlab = '\u4EBA\u5747\u9884\u671F\u5BFF\u547D', ylab = '\u9AD8\u5B66\u5386\u8005\u4EBA\u6570')
est = bkde2D(x[, 4:5], apply(x[, 4:5], 2, dpik))
contour(est$x1, est$x2, est$fhat, add = TRUE, lty = "12")
for (i in 1:nrow(x)) text(x[i, 4], x[i, 5], rownames(x)[i],
    cex = 0.75, adj = attr(x, "adj")[i, ])
rug(x[, 4], 0.02, side = 3, col = "gray40")
rug(x[, 5], 0.02, side = 4, col = "gray40")
boxplot(x[, 4], horizontal = TRUE, pars = list(boxwex = 7000,
    staplewex = 0.8, outwex = 0.8), at = -6000, add = TRUE, notch = TRUE, col = "skyblue",
    xaxt = "n")
boxplot(x[, 5], at = 63, pars = list(boxwex = 1.4,
    staplewex = 0.8, outwex = 0.8), add = TRUE, notch = TRUE, col = "skyblue",
    yaxt = "n")
text(67, 60000, "2005", cex = 3.5, col = "gray")
