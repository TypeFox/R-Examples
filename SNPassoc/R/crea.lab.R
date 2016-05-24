`crea.lab` <-
function (x, pos.ini, cex, dist) 
{
    n <- nrow(x)
    nn <- dimnames(x)[[1]]
    text(pos.ini + 0.1, 1.2, "frequency", cex = cex, adj = 0)
    text(pos.ini + 0.5, 1.2, "percentage", cex = cex, adj = 0)
    for (i in 1:n) {
        control <- (i - 1) * dist
        text(pos.ini, 1 - control, nn[i], cex = cex)
        text(pos.ini + 0.4, 1 - control, x[i, 1], adj = 1, cex = cex)
        text(pos.ini + 0.8, 1 - control, formatC(x[i, 2], 2,2,format="f"), adj = 1, 
            cex = cex)
    }
}

