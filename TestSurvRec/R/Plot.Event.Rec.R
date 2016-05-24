Plot.Event.Rec <-
function (yy, xy = 1, xf = 1) {

    XL <- data.frame(yy)
    fit1 <-survfitr(Survr(id, time, event) ~ as.factor(group), 
        data = XL, type = "pe")
    y <- survfitr(Survr(id, time, event) ~ 1, data = XL, type = "pe")
    failed <- matrix(y$failed)
    censored <- matrix(y$censored)
    time <- matrix(y$time)
    n.event <- matrix(y$n.event)
    AtRisk <- matrix(y$AtRisk)
    m <- y$m
    k <- matrix(m)
    mm <- 0
    m <- 0
    nn <- 0
    xf <- xf
    for (z in 1:xf) {
        nom <- as.character(paste("Unit = ", z, "\n"))
        m <- k[z, 1]
        nn <- mm + 1
        timebase <- 0
        tiempocalendario1 = matrix(0, nrow = m + 2, ncol = 1)
        x1 <- matrix(0, nrow = m + 2, ncol = 1)
        x2 <- matrix(0, nrow = m + 2, ncol = 1)
        y1 <- matrix(0, nrow = m + 2, ncol = 1)
        y2 <- matrix(0, nrow = m + 2, ncol = 1)
        tiempocalendario1[1, 1] <- 0
        x1[1, 1] <- 0
        x2[1, 1] <- 0
        y1[1, 1] <- 0
        y2[1, 1] <- 0
        for (i in 2:(m + 2)) {
            if (i < m + 2) {
                timebase <- timebase + failed[i - 1 + mm, 1]
                tiempocalendario1[i, 1] <- timebase
            }
            if (i == (m + 2)) {
                timebase <- timebase + censored[z, 1]
                tiempocalendario1[i, 1] <- timebase
            }
        }
        rx <- max(tiempocalendario1[1:(m + 2), 1])
        ry <- max(failed[nn:(nn + m - 1), 1], censored[z, 1])
        plot(0:1 * rx, 0:1 * ry, xlab = "Calendar time", ylab = "Gap time")
        title(main = nom)
        abline(h = 0, col = gray(0.9))
        abline(v = tiempocalendario1[(m + 2), 1], col = gray(0.9))
        r <- 1
        for (r in (1:(m + 2))) {
            if (r < (m + 2 - 1)) {
                x1[r, 1] <- tiempocalendario1[r, 1]
                y1[r, 1] <- tiempocalendario1[1, 1]
                x2[r, 1] <- tiempocalendario1[(r + 1), 1]
                y2[r, 1] <- failed[(r + mm), 1]
                if (r < (m + 2 - 1)) {
                  segments(x1[r, 1], y1[r, 1], x2[r, 1], y2[r, 
                    1], lty = 1, col = "red")
                  abline(h = y2[r, 1], col = gray(0.9))
                  abline(v = x2[r, 1], col = gray(0.9))
                  text(x2[r, 1], y2[r, 1], labels = y2[r, 1], 
                    adj = c(0, 1), cex = 0.7, pos = "2")
                }
            }
            if (r == (m + 2)) {
                x1[r, 1] <- tiempocalendario1[(r - 1), 1]
                y1[r, 1] <- tiempocalendario1[1, 1]
                x2[r, 1] <- tiempocalendario1[(m + 2), 1]
                y2[r, 1] <- censored[z, 1]
                segments(x1[r, 1], y1[r, 1], x2[r, 1], y2[r, 
                  1], lty = 1, col = "red")
                abline(h = y2[r, 1], col = gray(0.9))
                abline(v = x2[r, 1], col = gray(0.9))
                text(x2[r, 1], y2[r, 1], labels = y2[r, 1], adj = c(0, 
                  1), cex = 0.7, pos = "2")
            }
        }
        mm <- mm + m
    }
}
