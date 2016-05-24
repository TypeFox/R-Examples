par(xaxt = "n", yaxt = "n", mar = c(2, 0.1, 3, 0.1))
y = seq(2 * pi, 0, length = 150)
x = sin(y)
plot(x, y, type = "n", xlim = c(-1.3, 1.3), ylim = c(-0.93, 
    7.21), xlab = "", ylab = "", main = "IN THE ORBIT OF STATISTICS", 
    asp = 1)
rx = seq(-2, 2, length = 300)
ry = dnorm(rx, 0, 0.2) + y[75]
rty = ry - y[75]
agl = pi/4 + ifelse(rty/rx > 0, atan(rty/rx), atan(rty/rx) + 
    pi)
r = sqrt(rx^2 + rty^2)
rtx = r * cos(agl)
rty = r * sin(agl) + y[75]
lines(rtx, rty)
segments(-0.45, y[75] - 0.45, 0.45, y[75] + 0.45, 
    col = "white", lwd = 2)
q95 = qnorm(0.025, 0, 0.2)
d95 = dnorm(q95, 0, 0.2)
rg = rx > q95 & rx < -q95
polygon(c(rtx[rg][1], rtx[rg], rtx[rg][length(rg)]), 
    c(rty[rg][1], rty[rg], rty[rg][length(rg)]), col = "lavender")
radius = ifelse(0.5 * abs(cos(y)) + 0.1 < 0.15, 0.15, 
    0.5 * abs(cos(y)) + 0.1)
symbols(x, y, circles = radius, inches = FALSE, fg = heat.colors(150), 
    add = TRUE)
idx = rep(c(TRUE, FALSE), 75)
segments(x[idx], y[idx], x[!idx], y[!idx], col = heat.colors(200)[1:150][idx])
arrows(x[1], y[1], x[1] + 0.75, y[1] + 0.75, col = "red", 
    lwd = 2, angle = 75)
arrows(x[150] - 0.75, y[150] - 0.75, x[150], y[150], 
    col = "yellow", lwd = 2, angle = 75)
mtext("``STAT'' made by Yihui XIE using R, School of Statistics, RUC, 2007", 
    side = 1, cex = 0.8) 
