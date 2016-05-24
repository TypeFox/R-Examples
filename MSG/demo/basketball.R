if (require('animation')) {
    data(CLELAL09)
    draw.court = function() {
        rect(0, 0, 94, 50)
        circle = function(x, y, r, from = 0, to = 2 * pi, lines = FALSE, ...) {
            theta = seq(from, to, length = 100)
            if (lines)
                lines(x + r * cos(theta), y + r * sin(theta), ...)
            else polygon(x + r * cos(theta), y + r * sin(theta), ...)
        }
        points(c(5.25, 94 - 5.25), c(25, 25), cex = 2)
        segments(47, 0, 47, 50)
        circle(47, 25, 8)
        circle(47, 25, 2, col = "lightgray")
        theta1 = acos((25 - 35/12)/23.75)
        circle(5.25, 25, 23.75, -pi/2 + theta1, pi/2 - theta1, TRUE)
        circle(94 - 5.25, 25, 23.75, pi/2 + theta1, 3 * pi/2 - theta1, TRUE)
        segments(0, 35/12, 5.25 + 23.75 * sin(theta1), 35/12)
        segments(0, 50 - 35/12, 5.25 + 23.75 * sin(theta1), 50 - 35/12)
        segments(94, 35/12, 94 - 5.25 - 23.75 * sin(theta1), 35/12)
        segments(94, 50 - 35/12, 94 - 5.25 - 23.75 * sin(theta1), 50 - 35/12)
        circle(19, 25, 6, -pi/2, pi/2, TRUE)
        circle(19, 25, 6, pi/2, 3 * pi/2, TRUE, lty = 2)
        circle(94 - 19, 25, 6, pi/2, 3 * pi/2, TRUE)
        circle(94 - 19, 25, 6, -pi/2, pi/2, TRUE, lty = 2)
        circle(5.25, 25, 4, -pi/2, pi/2, TRUE)
        circle(94 - 5.25, 25, 4, pi/2, 3 * pi/2, TRUE)
        rect(0, 17, 19, 33, border = "gray")
        rect(94, 17, 94 - 19, 33, border = "gray")
    }
    with(CLELAL09, {
        jitterx = jitter(realx, amount = 1)
        jittery = jitter(realy, amount = 1)
        par(mar = c(0.05, 1, 0.05, 1), xpd = TRUE)
        smoothScatter(realx, realy, postPlotHook = NULL, xlim = c(0, 94), ylim = c(0, 50),
                      asp = 1, xaxs = "i", yaxs = "i", axes = F)
        draw.court()
        text(0, 25, "CLE", adj = c(0.5, 0), srt = 90)
        text(94, 25, "LAL", adj = c(0.5, 1), srt = 90)
        points(jitterx, jittery, pch = c(1, 4)[as.integer(result)])
    })
}
