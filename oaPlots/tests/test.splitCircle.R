library(oaPlots)
library(oaColors)

pdf("testCircles.pdf", width = 8, height = 8)
blankPlot( xlim = c(0, 4), ylim = c(0,4))



x <- 0.5; y <- 0.5; radius <- 0.48; col = oaColors("blue"); lty <- 1; lwd <- 1; nv <- 1000; r <- 0.5; border = NA
col1 <- oaColors("blue")
col2 <- oaColors("red")
splitAngle <- pi/4


splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 1.5; splitAngle <- pi/2
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 2.5; splitAngle <- 3*pi/4
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 3.5; splitAngle <- pi
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)


col1 <- oaColors("orange"); col2 <- oaColors("pink")
y <- 1.5; splitAngle <- pi/4; x <- 0.5
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 1.5; splitAngle <- pi/2
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 2.5; splitAngle <- 3*pi/4
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 3.5; splitAngle <- pi
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)


col1 <- oaColors("yellow"); col2 <- oaColors("green")
y <- 2.5; splitAngle <- pi/4; x <- 0.5
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 1.5; splitAngle <- pi/2
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 2.5; splitAngle <- 3*pi/4
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 3.5; splitAngle <- pi
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)


col1 <- oaColors("limegreen"); col2 <- oaColors("cyan")
y <- 3.5; splitAngle <- pi/4; x <- 0.5
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 1.5; splitAngle <- pi/2
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 2.5; splitAngle <- 3*pi/4
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)

x <- 3.5; splitAngle <- pi
splitCircle(x = x, y = y, radius = radius, splitAngle = splitAngle, nv = nv, border = border,
		col1 = col1, col2 = col2)
dev.off()
