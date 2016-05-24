shapes3dShapes <- function (x,loop=0,type="p",color=2,joinline=c(1:1),axes3=FALSE,rglopen=TRUE,main=main){ 
    if (is.matrix(x)) {
        xt <- array(0, c(dim(x), 1))
        xt[, , 1] <- x
        x <- xt
    }
    if (is.array(x) == FALSE) {
      stop("Data not in right format : require an array \n")
    }
    if (is.array(x) == TRUE) {
        if (rglopen) {
            rgl.open()
            rgl.bg(color = "white")
        }
        if (dim(x)[2] == 2) {
            x <- as.3d(x)
        }
        if (loop == 0) {
            k <- dim(x)[1]
            sz <- centroid.size(x[, , 1])/sqrt(k)/30
            plotshapes3d(x, type = type, color = color, size = sz, 
                joinline = joinline)
            if (axes3) {
                axes3d(color = "black", cex = 0.1, nticks = 4) 
                title3d(main = main, xlab = "x", ylab = "y", zlab = "z", color = "black")
            }
        }
        if (loop > 0) {
            for (i in 1:loop) {
                plotshapestime3d(x, type = type)
            }
        }
    }
}
