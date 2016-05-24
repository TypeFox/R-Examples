DrawArc <- 
function (x = 1, y = NULL, radius = 1, angle1 = deg1 * pi/180, 
    angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 35, col = 1, 
    ...) 
{
    angle <- angle1 + seq(0, length = n) * (angle2 - angle1)/n
    p1x <- x + radius * cos(angle)
    p1y <- y + radius * sin(angle)
    angle <- angle1 + seq(length = n) * (angle2 - angle1)/n
    p2x <- x + radius * cos(angle)
    p2y <- y + radius * sin(angle)
    segments(p1x, p1y, p2x, p2y, col = col, ...)
}
