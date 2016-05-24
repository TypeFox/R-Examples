shade = function(x, y1, y2, col = hsv(alpha=0.5), border = NA, ...){
    x = c(x, rev(x))
    y = c(y1, rev(y2))
    polygon(x, y, col=col, border=border, ...)
}

