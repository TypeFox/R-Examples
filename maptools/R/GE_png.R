# Copyright (c) 2007 by Duncan Golicher, David Forrest and Roger Bivand
#
# GE_SpatialGrid: function collecting and processing metadata for opening
# a PNG device to plot a Spatial* object for export to GE.
#
# arguments: obj: the Spatial* object; asp: if NA will be set to the latitude
# corrected value; maxPixels: the maximum dimension of the output PNG.
#
# values: returns a list containing:
# height and width (passed to png());
# SG (a SpatialGrid object with the grid topology of the output PNG);
# asp (the aspect value used);
# and xlim and ylim taken from SG. 
#
# These include the subcell inflation caused by rounding up the
# aspect-adjusted height or width, so that bbox(SG) is larger that
# bbox(obj) in one and only one value

GE_SpatialGrid <- function(obj, asp=NA, maxPixels=600) {
    if (!extends(class(obj), "Spatial")) 
        stop("GE_SpatialGrid only works for class(es extending) Spatial")
    p4s <- proj4string(obj)
    if (is.na(p4s) || is.projected(obj))
        stop("Spatial* object must be in geographical coordinates")
    xlim <- bbox(obj)[1,]
    ylim <- bbox(obj)[2,]
    s <- ifelse(is.na(asp), cos((mean(ylim) * pi)/180), asp)
    res <- Sobj_SpatialGrid(obj, asp=s, maxDim=maxPixels)
    class(res) <- "GE_SG"
    res
}

Sobj_SpatialGrid <- function(obj, asp=1, maxDim=100, n=NULL) {
    if (!extends(class(obj), "Spatial")) 
        stop("Sobj_SpatialGrid only works for class(es extending) Spatial")
    p4s <- proj4string(obj)
    xlim <- bbox(obj)[1,]
    ylim <- bbox(obj)[2,]
    m_asp <- (diff(ylim)/diff(xlim)) / asp
    names(m_asp) <- NULL
    if (!is.null(n)) {
        if (m_asp < 1) maxDim <- ceiling(sqrt(n/m_asp))
        else maxDim <- ceiling(sqrt(n*m_asp))
    }
    mywidth <- myheight <- maxDim
    if (m_asp < 1) {
	myheight1 <- mywidth * m_asp
        myheight <- ceiling(myheight1)
        cellsize <- c(diff(xlim)/mywidth, diff(ylim)/myheight1)

    } else {
        mywidth1 <- myheight / m_asp
        mywidth <- ceiling(mywidth1)
        cellsize <- c(diff(xlim)/mywidth1, diff(ylim)/myheight)
    }
    cells.dim <- c(mywidth, myheight)
    cellcentre.offset <- c(xlim[1]+(0.5*cellsize[1]), 
        ylim[1]+(0.5*cellsize[2]))
    names(cellcentre.offset) <- c("x", "y")
    grd <- GridTopology(cellcentre.offset, cellsize, cells.dim)
    mySG <- SpatialGrid(grd, proj4string=CRS(p4s))

    res <- list(height=as.integer(myheight), width=as.integer(mywidth),
        SG=mySG, asp=m_asp, xlim=bbox(mySG)[1,], ylim=bbox(mySG)[2,])
    res
}

# kmlOverlay: function to write image bounding box to GE GroundOverlay
# and link to image file to a kml file.
#
# arguments: obj: a GE_SG object from GE_SpatialGrid; kmlfile: If not NULL
# the name of the kml file to be written; imagefile: the name of the PNG
# file containing the image - this should be either relative (same
# directory as kml file) or abosolute (fully qualified); name: the name
# used to describe the image overlay in GE.
#
# values: x is a character vector containing the generated lines of the
# kml file

kmlOverlay <- function(obj, kmlfile=NULL, imagefile=NULL, name="R image") {
    if (class(obj) != "GE_SG") 
        stop("obj must be of class GE_SG from function GE_SpatialGrid")
    if (is.na(proj4string(obj$SG)) || is.projected(obj$SG))
        stop("Spatial* object must be in geographical coordinates")
    if (is.null(imagefile)) {
        imagefile <- "<fill_in_later>"
        warning("image file name missing, edit in manually")
    }
    bbox <- bbox(obj$SG)
    W <- bbox[1,1] ; E <- bbox[1,2]
    S <- bbox[2,1] ; N <- bbox[2,2]

    kmlheader <- c("<?xml version='1.0' encoding='UTF-8'?>",
        "<kml xmlns='http://earth.google.com/kml/2.0'>", "<GroundOverlay>")
    kmname <- paste("<name>", name, "</name>", sep="")
    icon <- paste("<Icon><href>", imagefile,
         "</href><viewBoundScale>0.75</viewBoundScale></Icon>", sep="")
    latlonbox <- paste("<LatLonBox><north>",
        N, "</north><south>",
        S, "</south><east>",
        E, "</east><west>",
        W, "</west></LatLonBox>", sep="")
    footer <- "</GroundOverlay></kml>"

    x <- (kmlheader)
    x <- append(x, kmname)
    x <- append(x, icon)
    x <- append(x, latlonbox)
    x <- append(x, footer)
    if (!is.null(kmlfile)) cat(paste(x, sep="", collapse="\n"), 
        file=kmlfile, sep="")
    x
}

