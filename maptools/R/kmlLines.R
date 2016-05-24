# Generates kml file of lines for SpatialLinesDataFrame object

kmlLines <- function(obj = NULL, kmlfile = NULL, name = "R Lines", 
                     description = "", col = NULL, visibility = 1, lwd = 1, 
                     kmlname = "", kmldescription = "") {
  
  # Handle NULL object
  if (is.null(obj)) 
    return(list(header = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
                           "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
                           "<Document>", paste("<name>", kmlname, "</name>", sep = ""),
                           paste("<description><![CDATA[", kmldescription, "]]></description>", sep = "")), 
                footer = c("</Document>", "</kml>")))
  
  # Handle wrong data type
  if (class(obj) == "list") {
    if (!is.null(obj[[1]])) {
      if ( class(obj[[1]]) != "Lines" && class(obj[[1]]) != "SpatialLinesDataFrame" ) {
        stop("obj must be of class 'Lines' or 'SpatialLinesDataFrame' or a 'list' of either [package 'maptools']")      
      }
    }
  } else {
    if (class(obj) != "Lines" && class(obj) != "SpatialLinesDataFrame") {
      stop("obj must be of class 'Lines' or 'SpatialLinesDataFrame' or a 'list' of either [package 'maptools']")
    } else {
      # Put the "Lines" or "SpatialLinesDataFrame" inside a list for consistent code further down
      obj <- list(obj)
    }
  }
  
  # Color conversion
  col2kmlcolor <- function(col) paste(rev(sapply(col2rgb(col, TRUE), function(x) sprintf("%02x", x))), collapse = "")
  
  # Set up defaults for different sections
  kml <- kmlStyle <- ""  
  kmlHeader <- c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
                 "<kml xmlns=\"http://earth.google.com/kml/2.2\">", 
                 "<Document>",
                 paste("<name>", kmlname, "</name>", sep = ""),
                 paste("<description><![CDATA[", kmldescription, "]]></description>", sep = ""))
  kmlFooter <- c("</Document>", "</kml>")
  
  # Replicate items as neded so as to avoid "NA" if the user provided a list of SpatialLines
  # but only a single name or description
  name <- rep(name,length.out=length(obj))
  description <- rep(description,length.out=length(obj))
  visibility <- rep(visibility,length.out=length(obj))
  lwd <- rep(lwd,length.out=length(obj))
  col <- rep(col,length.out=length(obj))
  
  # Create a <PlaceMark> for each SpatialLines object
  # NOTE:  Use character vectors to store text and then paste it all together.
  # NOTE:  Do not use 'kml <- append(kml, next_section)' as pass-by-value greatly degrades performance
  
  placemarks <- vector("character", length(obj))

  for (list_index in seq(length(obj))) {
    
    spatialLines <- obj[[list_index]]
  
    if (!is.null(spatialLines)) {
      
      # Create a Placemark with a MultiGeometry
      placemarkHeader <- paste("<Placemark>\n",
                               "<name>",name[list_index],"</name>\n",
                               "<description><![CDATA[",description[list_index],"]]></description>\n",
                               "<visibility>",as.integer(visibility[list_index]),"</visibility>",sep="")
      
      # Create a style
      lineStyle <- paste("<Style><LineStyle><width>",lwd[list_index],"</width>",
                         "<color>",col2kmlcolor(col[list_index]),"</color></LineStyle></Style>",sep="")
                           
      # NOTE:  Insert the kmlStyle in the <Placemark> instead of having a separate <Style> section 
      # NOTE:  at the <Document> level with <styleUrl> at the the <Placemark> level.
      # TODO:  May want to use <styleUrl> with a separate <Style> section.
      
      # Create a LineString for each element of obj@lines
      lineStrings <- vector("character", length(spatialLines@lines))
      for (i in 1:length(spatialLines@lines)) {
        line <- spatialLines@lines[[i]]
        coordinates <- paste(coordinates(line@Lines[[1]])[, 1], coordinates(line@Lines[[1]])[, 2], sep = "," ,collapse="\n")
        lineStrings[i] <- paste("<LineString><coordinates>",coordinates,"</coordinates></LineString>",sep="\n")
      }
      
      multiGeometry <- paste("<MultiGeometry>",paste(lineStrings,collapse="\n"),"</MultiGeometry>",sep="\n")
      placemarkFooter <- "</Placemark>\n"

      placemarks[list_index] <- paste(placemarkHeader,lineStyle,multiGeometry,placemarkFooter,sep="\n")
      
    }
  }

  kml <- paste(placemarks,sep="",collapse="\n")

  # Write out the file or return the components
  if (!is.null(kmlfile)) 
    cat(paste(c(kmlHeader, "", kml, kmlFooter),  # placeholder for separate kmlStyle section
              sep = "", collapse = "\n"), "\n", file = kmlfile, sep = "")
  else list(style = kmlStyle, content = kml)
}
