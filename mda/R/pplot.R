pplot <-
function (x, g, colors, pch, add = FALSE, type = "p", ...) 
{
    g <- as.factor(g)
    cc <- as.numeric(g)
    np <- seq(levels(g))
    if (missing(colors)) 
        colors <- np + 1
    else colors <- rep(colors, length = length(np))
    if (missing(pch)) 
        pch <- paste(np)
    else pch <- rep(pch, length = length(np))
    if (!add) 
        plot(x, type = "n", ...)
    for (i in unique(cc)) points(x[cc == i, , drop = FALSE], col = colors[i], 
        pch = pch[i], type = type)
}

