map.world <- function(add = FALSE, main = "", ...)
{
    #
    # read in the map data
    #
    lon <- gwm[, 1]
    lat <- gwm[, 2]

    if (!add){
        opar <- par(pty="m")
        on.exit(par(opar))
        plot.new()
        plot.window(xlim=c(-pi,pi), ylim=c(-pi/2,pi/2))
    }
    
    #
    # compute the startings and the endings of the polygons (lines)
    #
    rmp1 <- cumsum(rmp$inc)
    rmp1 <- append(rmp1, 0, 0)

    #
    # plot the polygons (lines)
    #
    for(j in 1:(length(rmp1) - 1)) {
        i1 <- rmp1[j] + 1
        i2 <- rmp1[j + 1]
        lines(lon[i1:i2]*pi/180, lat[i1:i2]*pi/180, ...)
    }

    if (!missing(main))
        title(main=main)
    
    if(!add){
        title(xlab="Longitude", ylab="Latitude")
        axis(1, at=seq(-pi, pi, by=pi/3), seq(-180, 180, by=60))
        axis(2, at=seq(-pi/2, pi/2, by=pi/6), seq(-90, 90, by=30))
        box()
    }    
}




