"mcp" <- function(xy, percent=95, unin=c("m", "km"),
                  unout=c("ha", "km2", "m2"))
{

    ## Verifications
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should be of class SpatialPoints")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    pfs <- proj4string(xy)

    if (length(percent)>1)
        stop("only one value is required for percent")
    if (percent>100) {
        warning("The MCP is estimated using all relocations (percent>100)")
        percent<-100
    }
    unin <- match.arg(unin)
    unout <- match.arg(unout)



    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)!=1) {
            warning("xy should contain only one column (the id of the animals), id ignored")
            id <- factor(rep("a", nrow(as.data.frame(xy))))
        } else {
            id <- xy[[1]]
        }
    } else {
        id <- factor(rep("a", nrow(as.data.frame(xy))))
    }

    if (percent>100) {
	warning("The MCP is estimated using all relocations (percent>100)")
	percent<-100
    }



    if (min(table(id))<5)
        stop("At least 5 relocations are required to fit an home range")
    id<-factor(id)



    xy <- as.data.frame(coordinates(xy))


    ## Computes the centroid of the relocations for each animal
    r<-split(xy, id)
    est.cdg<-function(xy) apply(xy, 2, mean)
    cdg<-lapply(r,est.cdg)
    levid<-levels(id)

    res <- SpatialPolygons(lapply(1:length(r), function(i) {
        k<-levid[i]
	df.t<-r[[levid[i]]]
	cdg.t<-cdg[[levid[i]]]

        ## Distances from the relocations to the centroid: we keep
        ## the "percent" closest
	dist.cdg<-function(xyt) {
            d<-sqrt( ( (xyt[1]-cdg.t[1])^2 ) + ( (xyt[2]-cdg.t[2])^2 ) )
            return(d)
        }

	di<-apply(df.t, 1, dist.cdg)
	key<-c(1:length(di))

	acons<-key[di<=quantile(di,percent/100)]
	xy.t<-df.t[acons,]


        ## Coordinates of the MCP
	coords.t<-chull(xy.t[,1], xy.t[,2])
	xy.bord<-xy.t[coords.t,]
        xy.bord <- rbind(xy.bord[nrow(xy.bord),], xy.bord)
        so <- Polygons(list(Polygon(as.matrix(xy.bord))), k)
        return(so)
    }))


    are <- unlist(lapply(1:length(res), function(i) {
        .arcpspdf(res[i,])
    }))

    if (unin == "m") {
        if (unout == "ha")
            are <- are/10000
        if (unout == "km2")
            are <- are/1e+06
    }
    if (unin == "km") {
        if (unout == "ha")
            are <- are * 100
        if (unout == "m2")
            are <- are * 1e+06
    }

    df <- data.frame(id=unlist(lapply(1:nlevels(id), function(i) res[i]@polygons[[1]]@ID)), area=are)
    row.names(df) <- df[,1]
    res <- SpatialPolygonsDataFrame(res, df)
    if (!is.na(pfs))
        proj4string(res) <- CRS(pfs)
    return(res)
}

