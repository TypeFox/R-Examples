
kmlPoints <- function (obj = NULL, kmlfile = NULL, kmlname = "", kmldescription = "",
                       name = NULL, description = "", 
                       icon = "http://www.gstatic.com/mapspro/images/stock/962-wht-diamond-blank.png") {
    # Handle null object
    if (is.null(obj)) 
        return(list(header = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
            "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
            "<Document>", paste("<name>", kmlname, "</name>", 
                sep = ""), paste("<description><![CDATA[", kmldescription, 
                "]]></description>", sep = "")), footer = c("</Document>", 
            "</kml>")))

    # Handle wrong data type
    if (class(obj) != "SpatialPointsDataFrame") 
        stop("obj must be of class 'SpatialPointsDataFrame' [package 'sp']")

    # Handle null name
    if (is.null(name)) {
        name = c()
        for (i in 1:nrow(obj))
            name <- append(name, paste("site", i))
    }

    # Handle single value name, description and icon
    if (length(name) < nrow(obj)) {
        if (length(name) > 1)
            warning("kmlPoints: length(name) does not match nrow(obj). The first name will be replicated.")
        name <- rep(name,nrow(obj))
    }
    if (length(description) < nrow(obj)) {
        if (length(description) > 1)
            warning("kmlPoints: length(description) does not match nrow(obj). The first description will be replicated.")
        description <- rep(description,nrow(obj))
    }
    if (length(icon) < nrow(obj)) {
        if (length(icon) > 1)
            warning("kmlPoints: length(icon) does not match nrow(obj). Only the first one will be used.")
        icon <- icon[1]
    }

    # Set up defaults for different sections
    kml <- kmlStyle <- ""
    kmlHeader <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
                   "<kml xmlns=\"http://earth.google.com/kml/2.2\">",
                   "<Document>", 
             paste("<name>", kmlname, "</name>", sep = ""),
             paste("<description><![CDATA[", kmldescription, "]]></description>", sep = ""))
    kmlFooter <- c("</Document>", "</kml>")

    # Create a set of styles
    # NOTE:  Available icons are at: http://sites.google.com/site/gmapsdevelopment/
    # NOTE:  No checking is done to make sure the icon href actually exists.
    for (i in 1:length(icon)) {
        pt_icon_href = icon[i]
        kmlStyle <- append(kmlStyle, paste("<Style id=\"style", i, "\">", sep = ""))
        kmlStyle <- append(kmlStyle, "  <IconStyle>")
        kmlStyle <- append(kmlStyle, "    <Icon>")
        kmlStyle <- append(kmlStyle, paste("      <href>", pt_icon_href, "</href>", sep = ""))
        kmlStyle <- append(kmlStyle, "    </Icon>")
        kmlStyle <- append(kmlStyle, "  </IconStyle>")
        kmlStyle <- append(kmlStyle, "</Style>")
    }

    # Create the sequential list of Placemarks
    for (i in 1:nrow(obj)) {
        point <- obj[i,]
        pt_name = name[i]
        pt_description = description[i]
        pt_style <- paste("#style",ifelse(length(icon) == 1, 1, i),sep="")
        kml <- append(kml, "<Placemark>")
        kml <- append(kml, paste("  <name>", pt_name, "</name>", sep = ""))
        kml <- append(kml, paste("  <description><![CDATA[", pt_description, "]]></description>", sep = ""))
        kml <- append(kml, paste("  <styleUrl>", pt_style, "</styleUrl>", sep = ""))
        kml <- append(kml, "  <Point>")
        kml <- append(kml, "    <coordinates>")
        kml <- append(kml, paste(point@coords[1], point@coords[2], sep = ","))
        kml <- append(kml, "    </coordinates>")
        kml <- append(kml, "  </Point>")
        kml <- append(kml, "</Placemark>")
    }
    
    # Write out the file or return the compontents
    if (!is.null(kmlfile)) 
        cat(paste(c(kmlHeader, kmlStyle, kml, kmlFooter), sep = "", 
            collapse = "\n"), "\n", file = kmlfile, sep = "")
    else list(style = kmlStyle, content = kml)
}

