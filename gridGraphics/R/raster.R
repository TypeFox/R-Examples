
# C_raster(image, xl, yb, xr, yt, angle, interpolate, ...)

C_raster <- function(x) {
    dev.set(recordDev())
    par <- currentPar(x[-(1:8)])
    dev.set(playDev())
    depth <- gotovp(par$xpd)
    image <- x[[2]]
    xl <- tx(x[[3]], par)
    yb <- ty(x[[4]], par)
    xr <- tx(x[[5]], par)
    yt <- ty(x[[6]], par)
    angle <- x[[7]]
    interpolate <- x[[8]]
    pushViewport(viewport(xl, yb, xr - xl, yt - yb, default.units="native",
                          just=c("left", "bottom"), angle=angle,
                          name=grobname("raster-vp")))
    grid.raster(image, width=1, height=1, interpolate=interpolate,
                name=grobname("raster"))
    upViewport(depth + 1)    
}

