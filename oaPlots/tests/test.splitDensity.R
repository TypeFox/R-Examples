library(oaPlots)
library(RColorBrewer)

set.seed(1)

# example Code
par(plt = c(0,1,0,1))
data <- rnorm(1000)
x <- density(data)$x
y <- density(data)$y

colVec <- brewer.pal(9, "Blues")[3:8]
outerCol <- brewer.pal(9, "Blues")[9]

oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(x = x, y = y, colVec = colVec, split = c(-1, -2, 0, 2, 1), 
		outerCol = outerCol,
		yScale = 0.95, yshift = 0)

oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(x = x, y = y, colVec = colVec, split = c(0), 
		outerCol = outerCol,
		yScale = 0.95, yshift = 0)

oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(x, y, colVec = colVec, split = c(-100), 
		outerCol = outerCol,
		yScale = 0.95, yshift = 0)

# should give warning
oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(x, y, colVec = colVec, 
		split = c(-1, -2, 0, 1, 1.1, 1.6, 1.8, 2.1, 2.2), 
		outerCol = outerCol,
		yScale = 0.95, yshift = 0)

# should give warning
oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(densityObj = density(data), x = x, y = y, colVec = colVec,
		split = c(-1, -2, 0, 2, 1), outerCol = outerCol,
		yScale = 0.95, yshift = 0)

oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(densityObj = density(data), colVec = colVec,
		split = c(-1, -2, 0, 2, 1), outerCol = outerCol,
		yScale = 0.95, yshift = 0)

oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(densityObj = density(data), colVec = colVec[3],
		split = NULL, outerCol = outerCol,
		yScale = 0.95, yshift = 0)






# Example 2
data <- rnorm(1000)
x <- density(data)$x
y <- density(data)$y
colVec <- brewer.pal(9, "Reds")[3:8]
outerCol <- brewer.pal(9, "Reds")[9]

oaTemplate(xlim = range(x), ylim = c(0, 1), ygrid = 0, cex.axis = 1.2)
drawSplitDensity(x, y, colVec = colVec, split = c(-8), 
		outerCol = outerCol,
		yScale = 0.95, yshift = 0)