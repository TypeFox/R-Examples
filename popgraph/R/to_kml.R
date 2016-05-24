#' Converts population graph to KML file
#' 
#' This function takes a population graph that has been 'decorated' with
#'  sufficient spatial data to make a KML file from it for viewing in 
#'  GoogleEarth.
#' @param graph A \code{popgraph} object.
#' @param file The location to save the kml to (if not passed, it returns it)
#' @return The text of the KML file to be saved or viewed in the appropriate editor.
#' @author Rodney J. Dyer <rjdyer@@vcu.edu>
#' @export
to_kml <- function( graph, file ) {
  if( !inherits( graph, "popgraph") )
    stop("Cannot save a kml file from a popgraph that is not made from a popgraph...")

  if( !("Latitude" %in% list.vertex.attributes(graph) ))
    stop("Cannot save to a kml file without a 'Latitude' property for each node.")
  if( !("Longitude" %in% list.vertex.attributes(graph)))
    stop("Cannot save to a kml file without a 'Longitude' property for each node.")
  
  kml.placemark <- function(name,lat,lon){
    r <- "\t<Placemark>"
    r <- c(r,paste("\t\t<name>",name,"</name>",sep=""))
    r <- c(r,"\t\t<description>Graph Node</description>")
    r <- c(r,"\t\t<styleUrl>#nodeIcon</styleUrl>")
    r <- c(r,"\t\t<altitudeMode>relativeToGround</altitudeMode>")
    r <- c(r,"\t\t<Point>")
    r <- c(r, paste("\t\t\t<coordinates>",lat,",",lon,",50</coordinates>",sep=""))
    r <- c(r,"\t\t</Point>")
    r <- c(r,"\t</Placemark>")
    r <- paste(r,collapse="\n")
    return(r)
  }
  
  kml.line <- function(name,lat1,lat2,lon1,lon2){
    r <- "\t<Placemark>"
    r <- c(r,paste("\t\t<name>",name,"</name>",sep=""))
    r <- c(r,"\t\t<description>Graph Edge</description>")
    r <- c(r,"\t\t<styleUrl>#edgeStyle</styleUrl>")
    r <- c(r,"\t\t<altitudeMode>clampToGround</altitudeMode>")
    r <- c(r,"\t\t<LineString>")
    r <- c(r,"\t\t\t<tessellate>1</tessellate>")
    r <- c(r, paste("\t\t\t<coordinates>\n\t\t\t\t",lon1,",",lat1,",0\n\t\t\t\t",lon2,",",lat2,",0\n\t\t\t</coordinates>",sep=""))
    r <- c(r,"\t\t</LineString>")
    r <- c(r,"\t</Placemark>")
    r <- paste(r,collapse="\n")
    return(r)  
  }
  
  ret <- "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
  ret <- c(ret, "<kml xmlns=\"http://www.opengis.net/kml/2.2\">")
  ret <- c(ret, "\t<Document>")
  ret <- c(ret, paste("\t<name>PopGraph: ",deparse(substitute(graph)),"</name>",sep=""))
  ret <- c(ret, "\t<description>A spacial representation of a Population Graph</description>")
  
  # put in the node style
  ret <- c(ret, "\t<Style id=\"nodeIcon\">" )
  ret <- c(ret, "\t\t<IconStyle>")
  ret <- c(ret, "\t\t\t<Icon>")
  ret <- c(ret, "\t\t\t\t<href>http://maps.google.com/mapfiles/kml/pal2/icon18.png</href>")
  ret <- c(ret, "\t\t\t</Icon>")
  ret <- c(ret, "\t\t<LineStyle>")
  ret <- c(ret, "\t\t\t<width>3</width>")
  ret <- c(ret, "\t\t</LineStyle>")
  ret <- c(ret, "\t\t</IconStyle>")
  ret <- c(ret, "\t</Style>")
  
  # put in the edge style
  ret <- c(ret, "\t<Style id=\"edgeStyle\">")
  ret <- c(ret, "\t\t<LineStyle>\n\t\t\t<color>ffdd7f4c</color>\n\t\t\t<width>4</width>\n\t\t</LineStyle>")
  ret <- c(ret, "\t</Style>")
  
  # make the nodes
  ret <- c(ret, "\t<Folder>")
  ret <- c(ret, "\t<name>Nodes</name>")
  ret <- c(ret, "\t<description>All the nodes from the graph file.</description>")
  
  K <- length(V(graph))
  for( i in 1:K) {
    place <- kml.placemark( V(graph)$name[i], V(graph)$Longitude[i],V(graph)$Latitude[i] )
    ret <- c(ret, place)
    
  }
  ret <- c(ret,"\t</Folder>")
  
  
  # make the edges
  ret <- c(ret, "\t<Folder>")
  ret <- c(ret, "\t<name>Edges</name>")
  ret <- c(ret, "\t<description>All the edges from the graph file.</description>")
  
  A <- get.adjacency(graph)
  for(i in 1:K){
    for(j in i:K){
      if( A[i,j] > 0 ){
        name <- paste(V(graph)$name[i],V(graph)$name[j],sep="-")
        lat1 <- V(graph)$Latitude[i]
        lat2 <- V(graph)$Latitude[j]
        lon1 <- V(graph)$Longitude[i]
        lon2 <- V(graph)$Longitude[j]
        line <- kml.line(name,lat1,lat2,lon1,lon2)
        ret <- c(ret, line)
      }
    }
  }
  
  ret <- c(ret,"\t</Folder>")
  
  ret <- c(ret, "\t</Document>\n</kml>")
  ret <- paste( ret, collapse="\n")
  
  if( !missing(file))
    write(ret,file)
  else
    return("the graph contents")
}