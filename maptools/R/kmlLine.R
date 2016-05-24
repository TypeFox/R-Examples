
kmlLine <- function(obj = NULL, kmlfile = NULL, name = "R Line", 
    description = "", col = NULL, visibility = 1, lwd = 1, 
    kmlname = "", kmldescription = "") {
    if (is.null(obj)) 
        return(list(header = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
            "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
            "<Document>", paste("<name>", kmlname, "</name>", 
              sep = ""), paste("<description><![CDATA[", 
              kmldescription, "]]></description>", sep = "")), 
            footer = c("</Document>", "</kml>")))
    if (class(obj) != "Lines" && class(obj) != "SpatialLinesDataFrame") 
        stop("obj must be of class 'Lines' or 'SpatialLinesDataFrame' [package 'sp']")
    if (class(obj) == "SpatialLinesDataFrame") {
        if (length(obj@lines) > 1L) 
            warning(paste("Only the first Lines object with the ID '", 
              obj@lines[[1]]@ID, "' is taken from 'obj'", 
              sep = ""))
        obj <- obj@lines[[1]]
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
    kmlStyle <- append(kmlStyle, paste("<color>", col2kmlcolor(col), 
        "</color>", sep = ""))
    kmlStyle <- append(kmlStyle, "</LineStyle>")
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

    for (i in 1:length(obj@Lines)) {
        kml <- append(kml, "<LineString>")
        kml <- append(kml, "<tessellate>1</tessellate>")
        kml <- append(kml, "<coordinates>")
        kml <- append(kml, paste(coordinates(obj@Lines[[i]])[, 
            1], coordinates(obj@Lines[[i]])[, 2], sep = ","))
        kml <- append(kml, "</coordinates>")
        kml <- append(kml, "</LineString>")
    }
    kml <- append(kml, "</MultiGeometry>")
    kml <- append(kml, "</Placemark>")

    if (!is.null(kmlfile)) 
        cat(paste(c(kmlHeader, kmlStyle, kml, kmlFooter), 
            sep = "", collapse = "\n"), "\n", file = kmlfile, sep = "")
    else list(style = kmlStyle, content = kml)
}
