draw.radial.line <- function(start, end, center=c(0, 0), angle=0, deg=NA, expand=FALSE,
    col=NA, lwd=NA, ...) 
    {
    if (is.na(col))
        col <- par("col")
    if (is.na(lwd))
        lwd <- par("lwd")
    if (is.na(deg))
        deg <- angle*180/pi
    angle <- deg*pi/180
    cosang <- cos(angle)
    sinang <- sin(angle)
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
    if (!expand || end <= 0)
        {
        px <- c(start, end)
        lines(center[1]+px*cosang, center[2]+px*sinang*ymult, col=col, lwd=lwd, ...)
        }
    else
        {
        devunits <- dev.size("px")
        lwdend <- (lwd/2)*(xylim[2]-xylim[1])/devunits[1]
        px <- c(start, start, end, end)
        lwdstart <- lwdend * start/end
        py <- c(lwdstart, -lwdstart, -lwdend, lwdend)
        pxt <- center[1] + px*cosang - py*sinang*ymult
        pyt <- center[2] + px*sinang + py*cosang*ymult
        polygon(pxt, pyt, border=col, col=col)
        }
    }
