
# make './inst/extdata/wrens/vector_combined/wrens'
    require(rgeos)
    require(rgdal)

    f = rangeMapper:::rangeFiles(new("rangeFiles", dir = './inst/extdata/wrens/vector/' ))

    o = list()
    for (i in 1:nrow(f)) {
        x  = readOGR(f[i,'dsn'], f[i,'layer'])
        dx = x@data[1, ]
        row.names(dx) = i
        x = gBuffer(x, width = 0)
        x  = gUnionCascaded(x, id = i)
        o[[i]] = SpatialPolygonsDataFrame(x, dx)
        }

    br = do.call(rbind, o)

    br$sciname = f$layer

    writeOGR(br, '/home/mihai/GIT/rangeMapper/inst/extdata/wrens/vector_combined/', 'wrens', driver="ESRI Shapefile", overwrite_layer = TRUE)


