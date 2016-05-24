"rgl.Map" <-
function(Map,which,...) {
  if (missing(which)) which <- TRUE
  if(!requireNamespace('rgl', quietly = TRUE)) stop("This function depends on the 'rgl' package which is not available")

  lapply(Map$Shapes[which], function(shape) {
    long <- shape$verts[,1] * pi/180
    lat  <- pi/2 - shape$verts[,2] * pi/180

    z <- cos(long)*sin(lat)
    y <- cos(lat)
    x <- sin(long)*sin(lat)


    mapply(function(pfrom, pto) {
        tmp.i <- rep( seq(along=x[pfrom:pto]), each=2 )
        tmp.i <- c(tmp.i[-1], 1)

        rgl::rgl.lines(x[tmp.i], y[tmp.i], z[tmp.i], ...)
    }, shape$Pstart + 1, c(shape$Pstart[-1], shape$nVerts + 1))

  })

  invisible()
}

