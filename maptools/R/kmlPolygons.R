# Generates kml file of polygon(s) for SpatialPolygonDataFrame object

kmlPolygons <- function(obj = NULL, kmlfile = NULL, name = "KML Polygons", 
                        description = "", col = NULL, visibility = 1, lwd = 1, 
                        border = "white", kmlname = "", kmldescription = "") {
  
  # Handle NULL object
  if (is.null(obj)) 
    return(list(header = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
                           "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
                           "<Document>", paste("<name>", kmlname, "</name>", 
                           sep = ""), paste("<description><![CDATA[", 
                           kmldescription, "]]></description>", sep = "")), 
                footer = c("</Document>", "</kml>")))

  # Handle wrong data type
  if (class(obj) == "list") {
    if (!is.null(obj[[1]])) {
      if ( class(obj[[1]]) != "Polygons" && class(obj[[1]]) != "SpatialPolygonsDataFrame" ) {
        stop("obj must be of class 'Polygons' or 'SpatialPolygonsDataFrame' or a 'list' of either [package 'maptools']")      
      }
    }
  } else {
    if (class(obj) != "Polygons" && class(obj) != "SpatialPolygonsDataFrame") {
      stop("obj must be of class 'Polygons' or 'SpatialPolygonsDataFrame' or a 'list' of either [package 'maptools']")
    } else { 
      # Put the "Polygons" or "SpatialPolygonsDataFrame" inside a list for consistent code further down
      obj <- list(obj)
    }
  }
  
  # Color conversion
  col2kmlcolor <- function(col) paste(rev(sapply(col2rgb(col,TRUE), function(x) sprintf("%02x", x))), collapse = "")

  # Set up defaults for different sections
  kml <- kmlStyle <- ""

  kmlHeader <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
                 "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
                 "<Document>", paste("<name>", kmlname, "</name>", 
                 sep = ""), paste("<description><![CDATA[", kmldescription, 
                 "]]></description>", sep = ""))
  kmlFooter <- c("</Document>", "</kml>")

  # Replicate items as neded so as to avoid "NA" if the user provided a list of SpatialPolygons
  # but only a single name or description
  name <- rep(name,length.out=length(obj))
  description <- rep(description,length.out=length(obj))
  visibility <- rep(visibility,length.out=length(obj))
  lwd <- rep(lwd,length.out=length(obj))
  col <- rep(col,length.out=length(obj))
  border <- rep(border,length.out=length(obj))
  
  # Create a <PlaceMark> for each SpatialPolygons object
  # NOTE:  Use character vectors to store text and then paste it all together.
  # NOTE:  Do not use 'kml <- append(kml, next_section)' as pass-by-value greatly degrades performance
  
  placemarks <- vector("character", length(obj))

  for (list_index in seq(length(obj))) {
    
    spatialPolygons <- obj[[list_index]]
  
    if (!is.null(spatialPolygons)) {
      
      # Create a Placemark with a MultiGeometry
      placemarkHeader <- paste("<Placemark>\n",
                               "<name>",name[list_index],"</name>\n",
                               "<description><![CDATA[",description[list_index],"]]></description>\n",
                               "<visibility>",as.integer(visibility[list_index]),"</visibility>",sep="")

      # Create a style
      if (is.null(col)) {
        PolyStyle <- "<PolyStyle><fill>0</fill></PolyStyle>"
      } else {
        PolyStyle <- paste("<PolyStyle><color>",col2kmlcolor(col[list_index]),
                           "</color><fill>1</fill></PolyStyle>",sep="")
      }
      polygonStyle <- paste("<Style><LineStyle><width>",lwd[list_index],"</width>",
                            "<color>",col2kmlcolor(border[list_index]),"</color></LineStyle>",
                            PolyStyle,"</Style>",sep="")
                           
      # NOTE:  Insert the kmlStyle in the <Placemark> instead of having a separate <Style> section 
      # NOTE:  at the <Document> level with <styleUrl> at the the <Placemark> level.
      # TODO:  May want to use <styleUrl> with a separate <Style> section.

      # NOTE:  The <Placemark> we are currently working on can have multiple <Polygon> elements

      # Create a polygonString for each element of obj@polygons
      
      polygonStrings <- vector("character", length(spatialPolygons@polygons))
      for (i in 1:length(spatialPolygons@polygons)) {
        polygon <- spatialPolygons@polygons[[i]]@Polygons[[1]]
        coordinates <- paste(coordinates(polygon)[,1], coordinates(polygon)[,2], sep=",", collapse="\n")
        polygonStrings[i] <- paste("<Polygon>",
                                   ifelse(polygon@hole,"<innerBoundaryIs>","<outerBoundaryIs>"),
                                   "<LinearRing><coordinates>",coordinates,"</coordinates></LinearRing>",
                                   ifelse(polygon@hole,"</innerBoundaryIs>","</outerBoundaryIs>"),
                                   "</Polygon>",sep="\n")
      }
      
      multiGeometry <- paste("<MultiGeometry>",paste(polygonStrings,collapse="\n"),"</MultiGeometry>",sep="\n")
      placemarkFooter <- "</Placemark>\n"

      placemarks[list_index] <- paste(placemarkHeader,polygonStyle,multiGeometry,placemarkFooter,sep="\n")
      
    }

  }

  kml <- paste(placemarks,sep="",collapse="\n")

  # Write out the file or return the components
  if (!is.null(kmlfile)) 
    cat(paste(c(kmlHeader, "", kml, kmlFooter),  # placeholder for separate kmlStyle section
              sep = "", collapse = "\n"), "\n", file = kmlfile, sep = "")
  else list(style = kmlStyle, content = kml)

}
