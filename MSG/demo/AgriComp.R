source(system.file("extdata", "AgriComp.R", package = "MSG"))
x = AgriComp
library(maps)
layout(matrix(1:2, 2), height = c(1, 8))
par(mar = c(0, 1.3, 2, 1.3))
barplot(rep(1, 98), col = heat.colors(98), axes = FALSE,
    space = 0)
mtext(side = 4, round(max(x[, 1]), 3), cex = 0.8)
mtext(side = 2, round(min(x[, 1]), 3), cex = 0.8)
mtext(side = 3, "Agricultural products competitiveness",
    line = 0.5)
coeff = round((x[, 1] - min(x[, 1]))/(max(x[, 1]) -
    min(x[, 1])) * 97) + 1
par(mar = c(2, 4, 0, 0))
map(col = "gray80", lty = 2)
for (i in 1:nrow(x)) {
    map(region = x[i, 2], fill = TRUE, col = heat.colors(98)[coeff[i]],
        add = TRUE)
}
abline(h = 0, lty = 2)
rect(72.26818, -55.8565, 168.93766, 54.589, lty = 2,
    col = NULL, lwd = 1)
map.axes()
