gx.add.chisq <-
function (p = c(0.98, 0.95, 0.9), df = NULL, ifflip = FALSE, 
    cex = 0.6) 
{
    if (is.null(df)) 
        stop("df (no. of variables) must be defined\n")
    q <- qchisq(p, df)
    nq <- length(q)
    limits <- par("usr")
    for (i in 1:nq) {
        if (q[i] < limits[2]) {
            abline(v = q[i], lty = 3)
            xpos <- q[i] - (limits[2] - limits[1]) * 0.005
            ypos <- limits[3] + (limits[4] - limits[3]) * 0.03
            adj <- 1
            if (ifflip) {
                xpos <- q[i] + (limits[2] - limits[1]) * 0.005
                adj <- 0
            }
            text(xpos, ypos, paste("p =", p[i]), adj = adj, cex = cex)
        }
    }
    invisible()
}
