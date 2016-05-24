library(raadtools)
ice <- aggregate(readice("2014-10-30"), fact = 4, fun = mean)
library(graticule)
g <- graticule(c(xmin(ice) - res(ice)[1]/2, xFromCol(ice) + res(ice)[1]/2), c(ymin(ice) - res(ice)[2]/2, yFromRow(ice) + res(ice)[2]/2),  tiles = TRUE)
## fudge it
projection(g) <- projection(ice)
g$layer <- values(ice)
g <- g[!is.na(g$layer) & g$layer > 0, ]

## native
plot(g, col = sst.pal(g$layer /3))
p1 <- "+proj=laea +lon_0=147 +lat_0=-42  +ellps=WGS84"
plot(spTransform(g, CRS(p1)), col = sst.pal(g$layer /3), border = NA)
plot(graticule(proj = p1), add = TRUE)
p2 <- "+proj=laea +lon_0=0 +lat_0=42  +ellps=WGS84"
plot(spTransform(g, CRS(p2)), col = sst.pal(g$layer /3), border = NA)
plot(graticule(proj = p2), add = TRUE)
p3 <- "+proj=longlat +ellps=WGS84"
plot(spTransform(g, CRS(p3)), col = sst.pal(g$layer /3), border = NA)
plot(graticule(proj = p3), add = TRUE)

p4 <- "+proj=gnom +lon_0=147 +lat_0=-90 +ellps=WGS84"
plot(spTransform(g, CRS(p4)), col = sst.pal(g$layer /3), border = NA)
plot(graticule(lons = seq(-180, 180, by = 15), lats = seq(-85, -35, by = 5), proj = p4), add = TRUE)


p4 <- "+proj=lcc  +ellps=WGS84"
plot(spTransform(g, CRS(p4)), col = sst.pal(g$layer /3), border = NA)
plot(graticule(lons = seq(-180, 180, by = 15), lats = seq(-85, -35, by = 5), proj = p4), add = TRUE)


sst <- crop(aggregate(readsst(), fact = 8, fun = mean), extent(5, 175, -90, 90))

gs <- graticule(c(xmin(sst) - res(sst)[1]/2, xFromCol(sst) + res(sst)[1]/2),
                c(ymin(sst) - res(sst)[2]/2, yFromRow(sst) + res(sst)[2]/2),  tiles = TRUE)

gs$layer <- values(sst)
gs <- gs[!is.na(gs$layer),  ]



p5 <- projection(ice)
plot(spTransform(gs, CRS(p5)), col = sst.pal(gs$layer), border = NA)
plot(graticule(lons = seq(-180, 180, by = 15), lats = seq(-85, -35, by = 5), proj = p5), add = TRUE)



p6 <- "+proj=laea +lat_0=-90  +ellps=WGS84"
plot(spTransform(gs, CRS(p6)), col = sst.pal(gs$layer), border = NA)
plot(graticule(lons = seq(-180, 180, by = 15), lats = seq(-85, -35, by = 5), proj = p6), add = TRUE)


p7 <- "+proj=omerc +lat_0=-42 +lonc=147 +alpha=5  +ellps=WGS84"
x <- spTransform(gs, CRS(p7))
library(rgeos)
aa <- gArea(x, byid = TRUE)
plot(x[aa < quantile(aa, .975),  ], col = sst.pal(gs$layer), border = NA)
plot(graticule(lons = seq(-180, 180, by = 15), lats = seq(-85, -35, by = 5), proj = p7), add = TRUE)


p8 <- "+proj=ortho +lon_0=90  +ellps=WGS84"
plot(spTransform(gs, CRS(p8)), col = sst.pal(gs$layer), border = NA)
plot(graticule(lons = seq(0, 180, by = 15), lats = seq(-85, 85, by = 5), proj = p8), add = TRUE)

