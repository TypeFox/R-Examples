arc <-
function (c, r, v, theta, ...) 
{
    angles <- anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- c[1] + r * cos(seqang)
    y <- c[2] + r * sin(seqang)
    lines(x, y, ...)
}
