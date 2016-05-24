`draw.cusp.bifset` <-
function (rx = par("usr")[1:2], ry = par("usr")[3:4], xpos = min(rx) + 
    0.01 * diff(rx)[1], ypos = max(ry) - 0.01 * diff(ry)[1], 
    xscale =  0.1 * diff(rx), yscale =  0.1 * diff(ry) / xscale, 
    aspect = 1, mark = 1, col = hsv(0.7, s = 0.8, alpha = 0.5), 
    border = NA, density = NA, bifurcation.set.fill = gray(0.8), 
    background = hsv(0.1, s = 0.1, alpha = 0.5), ..., X) 
{
    sx = xscale
    sy = yscale
    x = xpos + sx
    y = ypos - 1 * sy
    rect(-1 * sx + x, -1 * sy + y, 1 * sx + x, 1 * sy + y, col = background, 
        ...)
    bif <- cusp.bifset(seq_range(-1:1))
    polygon(c(bif[, 2], rev(bif[, 3])) * sx + x, c(bif[, 1], 
        rev(bif[, 1])) * sy + y, col = bifurcation.set.fill)
    segments(0 * sx + x, -1 * sy + y, 0 * sx + x, 1 * sy + y, 
        lty = 3, col = gray(0.3))
    segments(-1 * sx + x, 0 * sy + y, 1 * sx + x, 0 * sy + y, 
        lty = 3, col = gray(0.3))
    tmp <- switch(as.character(mark), "0" = NA, "1" = polygon(c(bif[, 
        2], -1, -1) * sx + x, c(bif[, 1], 1, 0) * sy + y, col = col, 
        border = border, density = density, ...), "2" = polygon(c(bif[, 
        3], 1, 1) * sx + x, c(bif[, 1], 1, 0) * sy + y, col = col, 
        border = border, density = density, ...), "3" = polygon(c(-1, 
        -1, 1, 1) * sx + x, c(-1, 0, 0, -1) * sy + y, col = col, 
        border = border, density = density, ...), "4" = polygon(c(bif[, 
        2], rev(bif[, 3])) * sx + x, c(bif[, 1], rev(bif[, 1])) * 
        sy + y, col = col, border = border, density = density, 
        ...))
}

