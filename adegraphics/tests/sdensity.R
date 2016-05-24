library(adegraphics)
pdf("sdensity.pdf")

xx1 <- rnorm(1000, 1, 2)
yy1 <- rnorm(1000, 1, 2)

g1 <- s.density(cbind(xx1, yy1), paxes.draw = T, gridsize = c(40, 40))
g2 <- s.density(cbind(xx1, yy1), paxes.draw = T, gridsize = c(80, 80), col = colorRampPalette(c("red", "blue"))(58), storeData = FALSE, region = TRUE)
g3 <- s.density(cbind(yy1 + 3, xx1 + 3), gridsize = c(400, 400))
g4 <- s.density(cbind(xx1, yy1), paxes.draw = T, gridsize = c(200, 200), add = TRUE)
g5 <- s.density(cbind(c(rnorm(50000, 1, 1), rnorm(50000, -1, 1)), c(rnorm(50000, -1, 0.5), rnorm(50000, 1, 0.5))), paxes.draw = T, gridsize = c(200, 200), region = TRUE, contour = TRUE, plabels.cex = 1, plabels.srt = "vertial")
g6 <- s.density(cbind(rnorm(300, 3, 0.3), rnorm(300, -2, 0.5)), gridsize = c(500, 500), thres = 0.01, nr = 10, regions = list(alpha = 0.5), col = colorRampPalette(c("red", "blue"))(108))
g7 <- s.density(cbind(c(rnorm(50000, 1, 1), rnorm(50000, -1, 1)), c(rnorm(50000, -1, 0.5), rnorm(50000, 1, 0.5))), paxes.draw = T, gridsize = c(100, 100), region = TRUE, contour = TRUE, plabels.cex = 1, nclass = 5)
