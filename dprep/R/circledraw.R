circledraw <-
function (numpts = 200, radius = 1) 
{
    xy = rep(360/numpts, numpts)
    xy <- c(0, cumsum(xy))
    xy <- xy * pi/180
    t2xy = list(x = radius * cos(xy), y = radius * sin(xy))
    plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1] > pin[2]) 
        xlim <- (pin[1]/pin[2]) * xlim
    else ylim <- (pin[2]/pin[1]) * ylim
    plot.window(xlim, ylim, "", asp = 1)
    polygon(t2xy$x, t2xy$y, border = "black", col = "white")
    invisible(NULL)
}
