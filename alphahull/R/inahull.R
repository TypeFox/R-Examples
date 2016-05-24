inahull <-
function (ahull.obj, p) 
{
    p <- matrix(p, ncol = 2)
    compl <- ahull.obj$complement
    halfpl <- which(compl[, "r"] < 0)
    n.halfpl <- length(halfpl)
    ball <- which(compl[, "r"] > 0)
    n.ball <- length(ball)
    in.compl <- logical(length = dim(p)[1])
    if (n.halfpl >= 1) {
        h <- 1
        while ((h <= n.halfpl)) {
            sig = compl[halfpl[h], 3]
            a = compl[halfpl[h], 1]
            b = compl[halfpl[h], 2]
            if (sig <= -3) {
                in.compl[p[, 1] > a & sig == 3] <- TRUE
                in.compl[p[, 1] < a & sig == -4] <- TRUE
            }
            else {
                in.compl[p[, 2] > (a + b * p[, 1]) & sig == -1] <- TRUE
                in.compl[p[, 2] < (a + b * p[, 1]) & sig == -2] <- TRUE
            }
            h <- h + 1
        }
    }
    wp <- which(!in.compl)
    k <- 1
    while ((k <= n.ball)) {
        r = compl[ball[k], 3]
        c1 = compl[ball[k], 1]
        c2 = compl[ball[k], 2]
        d <- sqrt((p[wp, 1] - c1)^2 + (p[wp, 2] - c2)^2)
        in.compl[wp[d < r]] <- TRUE
        k <- k + 1
    }
    return(in.ahull = !in.compl)
}
