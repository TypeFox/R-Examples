# A. Data: number of plots, points, and colors
setwd("C:/aErer"); library(copula); library(grid)
num.plot <- 500; num.points <- 10000
set.rho <- seq(from = 0, to = 1, length.out = num.plot)
pie(x = rep(x = 1, times = 15), col = rainbow(15))  # understand rainbow()
set.col <- rainbow(num.plot)
str(set.col); head(set.col, n = 4)

# B. Display the video on the screen device; need about 35 seconds
# Two random variables with bivariate normal copula
windows(width = 4, height = 4); bringToTop(stay = TRUE)
par(mar = c(2.5, 2.5, 1, 1))
for (i in 1:num.plot) {
  sam <- rCopula(n = num.points, copula = normalCopula(param = set.rho[i]))
  plot(x = sam, col = set.col[i], xlim = c(0, 1), ylim = c(0, 1), 
    type = "p", pch = ".")
}
str(sam); head(sam, n = 3)
    
# C. Three screenshots for ERER book
windows(width = 5.3, height = 2.5); bringToTop(stay = TRUE)
v1 <- viewport(x = 0.02, y = 0.98, width = 0.55, height = 0.5, 
  just = c("left", "top"))
pushViewport(v1); grid.rect(gp = gpar(lty = "dashed"))
sam <- rCopula(n = 3000, copula = normalCopula(param = 0))
grid.points(x = sam[, 1], y = sam[, 2], pch = 1, 
  size = unit(0.001, "char"), gp = gpar(col = "red"))
          
upViewport(0); current.viewport()
v2 <- viewport(width = 0.55, height = 0.5)
pushViewport(v2); grid.rect(gp = gpar(fill = "white", lty = "dashed"))
sam <- rCopula(n = 3000, copula = normalCopula(param = 0.8))
grid.points(x = sam[, 1], y = sam[, 2], pch = 1, 
  size = unit(0.001, "char"), gp = gpar(col = "darkgreen"))
   
upViewport(0)
v3 <- viewport(x = 0.98, y = 0.02, width = 0.55, height = 0.5, 
  just = c("right", "bottom"))
pushViewport(v3); grid.rect(gp = gpar(fill = "white", lty = "dashed"))
sam <- rCopula(n = 3000, copula = normalCopula(param = 0.99))
grid.points(x = sam[, 1], y = sam[, 2], pch = 1, 
  size = unit(0.001, "char"), gp = gpar(col = "blue"))
showCorrelation <- recordPlot()

# D. Save the three screen shots on a file device
pdf(file = "fig_showCorrelation.pdf", width = 5.3, height = 2.5,
  useDingbats = FALSE)
replayPlot(showCorrelation); dev.off() 