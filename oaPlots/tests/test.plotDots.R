library(oaPlots)

x <- sample(1:5, size = 25, replace = TRUE)
plot(x = -1, y = -1, xlim = c(0.5,1.5), ylim = range(x), 
    ylab = "", xlab = "")
colVec <- c(rep("olivedrab", 15), rep("goldenrod", 5), rep("red", 5))

boxplot(x)
plotDots(vec = x, xLeft = 0.9, xRight = 1.1, pch = 19, 
    col = colVec, cex = 2)
