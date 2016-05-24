"mcp.area" <- function(xy, id, percent = seq(20,100, by=5),
                       unin=c("m", "km"),
                       unout=c("ha", "km2", "m2"), plotit = TRUE)
{
    ## Verifications
    unin<-match.arg(unin)
    unout<-match.arg(unout)
    if (length(id) != nrow(xy))
        stop("xy and id should be of the same length")

    ## remove the missing values
    xy <- xy[!is.na(xy[, 1]), ]
    xy <- xy[!is.na(xy[, 2]), ]
    id <- id[!is.na(xy[, 1])]
    id <- id[!is.na(xy[, 2])]

    ## Bases
    lev<-percent
    res<-list()
    ar<-matrix(0,nrow=length(lev),
               ncol=nlevels(factor(id)))
    lixy <- split(xy, id)
    le <- names(lixy)

    ## For each home range level, computes the MCP, and its area
    for (i in 1:length(lev)) {

        ## Computes the MCP for each animal
        ar[i,] <- unlist(lapply(lixy, function(z) {
            res<-mcp(z, rep(1,nrow(z)), percent=lev[i])
            class(res)<-"data.frame"
            if (!all(unlist(res[1,])==unlist(res[nrow(res),])))
                res <- rbind(res,res[1,])
            pol <- Polygon(as.matrix(res[,2:3]))
            spdf <- SpatialPolygons(list(Polygons(list(pol), 1)))
            lar <- unlist(lapply(polygons(spdf)@polygons,
                                 function(x) unlist(lapply(x@Polygons, function(y)
                                                           .arcp(y@coords)))))
            lhol <- unlist(lapply(polygons(spdf)@polygons,
                                  function(x) unlist(lapply(x@Polygons, function(y)
                                                            y@hole))))
            sum(lar[!lhol])-sum(lar[lhol])
        }))

    }

    ar <- as.data.frame(ar)
    names(ar)<-le

    ## output units
    if (unin=="m") {
        if (unout=="ha")
            ar<-ar/10000
        if (unout=="km2")
            ar<-ar/1000000
    }
    if (unin=="km") {
        if (unout=="ha")
            ar<-ar*100
        if (unout=="m2")
            ar<-ar*1000000
    }

    ## output
    row.names(ar)<-lev
    class(ar)<-c("hrsize", "data.frame")
    attr(ar, "units")<-unout
    if (plotit)
        plot(ar)
    return(ar)
}

