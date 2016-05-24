"hexgrid" <-
function(x=0, y=0, r=100, nro=10, nco=20) {
    odd <- function(x) {x != as.integer(x/2) * 2}
    even <- function(x) {x == as.integer(x/2) * 2}
    delta.y <- sqrt(r^2-(r/2)^2)
    delta.x <- r/2
    xm <- matrix(x, nrow=nro, ncol=nco)
    for (i in 2:nco) {
        xm[,i] <- xm[,(i-1)]+delta.x
    }
    xm <- data.frame(xm)
    names(xm) <- c(1:nco)
    rownames(xm) <- c(1:nro)
    xmodd <- xm[odd(as.numeric(rownames(xm))),odd(as.numeric(names(xm)))]
    xmeven <- xm[even(as.numeric(rownames(xm))),even(as.numeric(names(xm)))]
    xmall <- c(as.vector(as.matrix(xmodd)),as.vector(as.matrix(xmeven)))
    ym <- matrix(y, nrow=nro, ncol=nco)
    for (i in 2:nro) {
        ym[i,] <- ym[(i-1),]-delta.y
    }
    ym <- data.frame(ym)
    names(ym) <- c(1:nco)
    rownames(ym) <- c(1:nro)
    ymodd <- ym[odd(as.numeric(rownames(ym))),odd(as.numeric(names(ym)))]
    ymeven <- ym[even(as.numeric(rownames(ym))),even(as.numeric(names(ym)))]
    ymall <- c(as.vector(as.matrix(ymodd)),as.vector(as.matrix(ymeven)))
    am <- matrix("1", nrow=nro, ncol=nco)
    for (i in 2:nro) {
        am[i,] <- i
    }
    am <- data.frame(am)
    names(am) <- c(1:nco)
    rownames(am) <- c(1:nro)
    amodd <- am[odd(as.numeric(rownames(am))),odd(as.numeric(names(am)))]
    ameven <- am[even(as.numeric(rownames(am))),even(as.numeric(names(am)))]
    amall <- c(as.vector(as.matrix(amodd)),as.vector(as.matrix(ameven)))
    bm <- matrix("1", nrow=nro, ncol=nco)
    for (i in 2:nco) {
        bm[,i] <- i
    }
    bm <- data.frame(bm)
    names(bm) <- c(1:nco)
    rownames(bm) <- c(1:nro)
    bmodd <- bm[odd(as.numeric(rownames(bm))),odd(as.numeric(names(bm)))]
    bmeven <- bm[even(as.numeric(rownames(bm))),even(as.numeric(names(bm)))]
    bmall <- c(as.vector(as.matrix(bmodd)),as.vector(as.matrix(bmeven)))
    coords <- data.frame(amall, bmall, xmall, ymall)
    names(coords) <- c("ROW", "COL", "X", "Y")
    coords
}
