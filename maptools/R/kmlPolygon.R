
kmlPolygon <- function(obj = NULL, kmlfile = NULL, name = "R Polygon", 
    description = "", col = NULL, visibility = 1, lwd = 1, 
    border = 1, kmlname = "", kmldescription = "") {
    if (is.null(obj)) 
        return(list(header = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
            "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
            "<Document>", paste("<name>", kmlname, "</name>", 
              sep = ""), paste("<description><![CDATA[", 
              kmldescription, "]]></description>", sep = "")), 
            footer = c("</Document>", "</kml>")))
    if (class(obj) != "Polygons" && class(obj) != "SpatialPolygonsDataFrame") 
        stop("obj must be of class 'Polygons' or 'SpatialPolygonsDataFrame' [package 'sp']")
    if (class(obj) == "SpatialPolygonsDataFrame") {
        if (length(obj@polygons) > 1L) 
            warning(paste("Only the first Polygons object with the ID '", 
              obj@polygons[[1]]@ID, "' is taken from 'obj'", 
              sep = ""))
        obj <- obj@polygons[[1]]
    }

    col2kmlcolor <- function(col) paste(rev(sapply(col2rgb(col, 
        TRUE), function(x) sprintf("%02x", x))), collapse = "")

    kml <- kmlStyle <- ""

    kmlHeader <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
        "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
        "<Document>", paste("<name>", kmlname, "</name>", 
            sep = ""), paste("<description><![CDATA[", kmldescription, 
            "]]></description>", sep = ""))
    kmlFooter <- c("</Document>", "</kml>")

    kmlStyle <- append(kmlStyle, paste("<Style id=\"", obj@ID, 
        "\">", sep = ""))
    kmlStyle <- append(kmlStyle, "<LineStyle>")
    kmlStyle <- append(kmlStyle, paste("<width>", lwd, "</width>", 
        sep = ""))
    kmlStyle <- append(kmlStyle, paste("<color>", col2kmlcolor(border), 
        "</color>", sep = ""))
    kmlStyle <- append(kmlStyle, "</LineStyle>")
    kmlStyle <- append(kmlStyle, "<PolyStyle>")

    if (is.null(col)) {
        kmlStyle <- append(kmlStyle, "<fill>0</fill>")
    }
    else {
        kmlStyle <- append(kmlStyle, paste("<color>", col2kmlcolor(col), 
            "</color>", sep = ""))
        kmlStyle <- append(kmlStyle, "<fill>1</fill>")
    }

    kmlStyle <- append(kmlStyle, "</PolyStyle>")
    kmlStyle <- append(kmlStyle, "</Style>")

    kml <- append(kml, "<Placemark>")
    kml <- append(kml, paste("<name>", name, "</name>", sep = ""))
    kml <- append(kml, paste("<description><![CDATA[", description, 
        "]]></description>", sep = ""))
    kml <- append(kml, paste("<styleUrl>#", obj@ID, "</styleUrl>", 
        sep = ""))
    kml <- append(kml, paste("<visibility>", as.integer(visibility), 
        "</visibility>", sep = ""))
    kml <- append(kml, "<MultiGeometry>")
    holeFlag <- FALSE
    for (i in 1:length(obj@Polygons)) {
        if (!holeFlag) 
            kml <- append(kml, "<Polygon>")
        kml <- append(kml, ifelse(obj@Polygons[[i]]@hole, 
            "<innerBoundaryIs>", "<outerBoundaryIs>"))
        kml <- append(kml, c("<LinearRing>", "<coordinates>"))
        kml <- append(kml, paste(coordinates(obj@Polygons[[i]])[, 1], 
            coordinates(obj@Polygons[[i]])[, 2], sep = ","))
        kml <- append(kml, c("</coordinates>", "</LinearRing>"))
        kml <- append(kml, ifelse(obj@Polygons[[i]]@hole, 
            "</innerBoundaryIs>", "</outerBoundaryIs>"))
        # check if next polygon is an hole, if so export it as innerBoundaryIs;
        # this only works if the holes are following that polygon 
        # which contains these holes regardingless plotOrder;
        # TODO: rearrange holes according to 'their' polygons automatically via
        # hole.in.which.polygon?
        holeFlag <- ifelse((i + 1L) <= length(obj@Polygons), 
            obj@Polygons[[i + 1L]]@hole, FALSE)
        if (!holeFlag) 
            kml <- append(kml, "</Polygon>")
    }
    kml <- append(kml, "</MultiGeometry>")
    kml <- append(kml, "</Placemark>")
    
    if (!is.null(kmlfile)) 
        cat(paste(c(kmlHeader, kmlStyle, kml, kmlFooter), 
            sep = "", collapse = "\n"), "\n", file = kmlfile, sep = "")
    else list(style = kmlStyle, content = kml)
}
