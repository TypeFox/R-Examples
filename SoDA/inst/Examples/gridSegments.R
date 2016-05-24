
require(grid)
require(SoDA)

load("Examples/with_gps1.rda")

pdf("Examples/gridSegments.pdf", width = 4, height=2)

with(gps1, {

    xy = geoXY(latitude, longitude)
    x = xy[,1]; y = xy[,2]
    vlim = range(elevation, na.rm=TRUE)
    v = 1.5 * (elevation - vlim[[1]])/diff(vlim) + 0.5
    xlim  = range(x, na.rm=TRUE); ylim  = range(y, na.rm=TRUE)
    ylim[1] <- ylim[1] + .5*  (diff(ylim))
    pushViewport(viewport(width = unit(4, "inches"), height = unit(2, "inches")))
    pushViewport(plotViewport(rep(2,4)))
    pushViewport(dataViewport(xscale = xlim, yscale = ylim))
    show <- y >= ylim[1]
    
    xn <- unit(x[show], "native"); yn <- unit(y[show], "native")
    segObject <- grid.segments(xn, yn, xn, yn +  unit(v[show], "lines"))
##   browser()

})

dev.off()
