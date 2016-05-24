"img.matcor" <-
function (correl, type = 1) 
{
    matcorX = correl$Xcor
    matcorY = correl$Ycor
    matcorXY = correl$XYcor
    lX = ncol(matcorX)
    lY = ncol(matcorY)
    def.par <- par(no.readonly = TRUE)
    if (type == 1) {
        par(mfrow = c(1, 1), pty = "s")
        image(1:(lX + lY), 1:(lX + lY), t(matcorXY[nrow(matcorXY):1, ]), 
            zlim = c(-1, 1), main = "XY correlation", 
            col = tim.colors(64), 
            axes = FALSE, , xlab = "", ylab = "")
        box()
        abline(h = lY + 0.5, v = lX + 0.5, lwd = 2, lty = 2)
        image.plot(legend.only = TRUE, zlim = c(-1, 1), 
            col = tim.colors(64), horizontal = TRUE)
    }
    else {
        layout(matrix(c(1, 2, 3, 3, 0, 0), ncol = 2, nrow = 3, 
             byrow = TRUE), widths = 1, heights = c(0.8, 1, 0.06))

    #layout 1 
        par(pty = "s", mar = c(2, 2, 2, 1))
        image(1:lX, 1:lX, t(matcorX[lX:1, ]), zlim = c(-1, 1), 
            main = "X correlation", axes = FALSE, xlab = "", ylab = "",
            col = tim.colors(64))
        box()

        image(1:lY, 1:lY, t(matcorY[lY:1, ]), zlim = c(-1, 1), 
            main = "Y correlation", 
            col = tim.colors(64), axes = FALSE, 
            xlab = "", ylab = "")
        box()

    #layout 2
        partXY = matcorXY[lX:1, (lX + 1):(lX + lY)]
        if (lX > lY) {
            partXY = matcorXY[(lX + lY):(lX + 1), 1:lX]
            lX = ncol(matcorY)
            lY = ncol(matcorX)
        }

        par(pty = "m", mar = c(5,4,2,3), mai = c(0.8, 0.5, 0.3, 0.4))
        image(1:lY, 1:lX, t(partXY), zlim = c(-1, 1), 
            main = "Cross-correlation", axes = FALSE, xlab = "", 
            ylab = "", col = tim.colors(64))
        box()

        image.plot(legend.only = TRUE, zlim = c(-1, 1), 
            horizontal = TRUE, col = tim.colors(64), 
            legend.width = 2.5)
    }
    par(def.par)
}

