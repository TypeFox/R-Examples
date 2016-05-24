library(plotrix)
par(mar = rep(1, 4))
plot(1:5, ylim=c(1,4),type = "n", ann = FALSE, axes = FALSE)
polygon(c(0, 0, 5.5, 5.5), c(0, 3, 3, 0), border = NA,
    col = "#44aaff")
pie.col = rgb(1:0, 1:0, c(0, 0), c(0.7, 1))
floating.pie(1.7, 3, c(0.8, 0.2), radius = 0.5, col = pie.col,
    startpos = -0.3 * pi)
floating.pie(3.1, 3, c(0.7, 0.3), radius = 0.5, col = pie.col,
    startpos = -0.2 * pi)
floating.pie(4, 1.5, c(0.15, 0.85), radius = 0.5,
    col = pie.col, startpos = 0.35 * pi)
points(c(3.9, 4.01, 4.07, 4.03), c(2.1, 2.3, 2.55,
    2.9), pch = 21, cex = 1.5, bg = "white")
text(c(1.7, 3.1, 4), c(3.7, 3.7, 3.7), c("\u5408\u683C", "\u5408\u683C",
    "\u4E0D\u5408\u683C"))
box()
