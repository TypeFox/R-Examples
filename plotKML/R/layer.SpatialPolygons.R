# Purpose        : Write polygon-objects (SpatialPolygons) to KML;
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : this operation can be time consuming for large grids!;

kml_layer.SpatialPolygons <- function(
  obj,
  subfolder.name = paste(class(obj)),
  extrude = TRUE,
  tessellate = FALSE,
  outline = TRUE,
  plot.labpt = FALSE,
  z.scale = 1,
  LabelScale = get("LabelScale", envir = plotKML.opts),
  metadata = NULL,
  html.table = NULL,
  TimeSpan.begin = "",
  TimeSpan.end = "",
  colorMode = "normal",
  ...
  ){
  
  # invisible file connection
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  # Checking the projection is geo
  prj.check <- check_projection(obj, control = TRUE)

  # Trying to reproject data if the check was not successful
  if (!prj.check) {  obj <- reproject(obj)  }

  # Parsing the call for aesthetics
  aes <- kml_aes(obj, ...)

  # Read the relevant aesthetics
  poly_names <- aes[["labels"]]
  colours <- aes[["colour"]]
  sizes <- aes[["size"]]
  shapes <- aes[["shape"]]
  altitude <- aes[["altitude"]]  # this only works if the altitudes have not been defined in the original sp class
  altitudeMode <- aes[["altitudeMode"]]
  balloon <- aes[["balloon"]]

  # Parse ATTRIBUTE TABLE (for each placemark):
  if((balloon==TRUE | class(balloon) %in% c('character','numeric')) & ("data" %in% slotNames(obj))){
      html.table <- .df2htmltable(obj@data)
  }
  
  # Folder and name of the points folder
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", subfolder.name, parent = pl1)

  if(plot.labpt==TRUE){
    pl1b = newXMLNode("Folder", parent=kml.out[["Document"]])
    pl2b <- newXMLNode("name", "labpt", parent = pl1b)
  }

  # Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txt <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txt, parent=pl1)
  }
  message("Writing to KML...")  

  # Prepare data for writing
  # ==============

  # number of polygons:
  pv <- length(obj@polygons)
  # number of Polygons:
  pvn <- lapply(lapply(obj@polygons, slot, "Polygons"), length)
  # parse coordinates:
  coords <- rep(list(NULL), pv)
  hole <- rep(list(NULL), pv)
  labpt <- rep(list(NULL), pv)
  for(i.poly in 1:pv) { 
    for(i.Poly in 1:pvn[[i.poly]]){
    # get coordinates / hole definition:
    xyz <- slot(slot(obj@polygons[[i.poly]], "Polygons")[[i.Poly]], "coords")
    cxyz <- slot(slot(obj@polygons[[i.poly]], "Polygons")[[i.Poly]], "labpt")
    # if altitude is missing, add the default altitudes:
    if(ncol(xyz)==2){  xyz <- cbind(xyz, rep(altitude[i.poly], nrow(xyz)))  }
    # format coords for writing to KML [https://developers.google.com/kml/documentation/kmlreference#polygon]:
    hole[[i.poly]][[i.Poly]] <- slot(slot(obj@polygons[[i.poly]], "Polygons")[[i.Poly]], "hole")
    coords[[i.poly]][[i.Poly]] <- paste(xyz[,1], ',', xyz[,2], ',', xyz[,3], collapse='\n ', sep = "")
    labpt[[i.poly]][[i.Poly]] <- paste(cxyz[1], ',', cxyz[2], ',', altitude[i.poly], collapse='\n ', sep = "")
    }
  }

  # reformatted aesthetics (one "polygons" can have multiple "Polygons"):
  poly_names.l <- list(NULL)
  for(i.poly in 1:pv){ poly_names.l[[i.poly]] <- as.vector(rep(poly_names[i.poly], pvn[[i.poly]])) }
  # polygon times (if applicable)
  TimeSpan.begin.l <- list(NULL)
  TimeSpan.end.l <- list(NULL)
  when.l <- list(NULL)
  # check if time span has been defined:
  if(all(nzchar(TimeSpan.begin))&all(nzchar(TimeSpan.end))){
      if(identical(TimeSpan.begin, TimeSpan.end)){
        if(length(TimeSpan.begin)==1){ 
          when.l = rep(TimeSpan.begin, sum(unlist(pvn))) 
        } else {
          for(i.poly in 1:pv){ when.l[[i.poly]] <- as.vector(rep(TimeSpan.begin[i.poly], pvn[[i.poly]])) }
          }} else {
            for(i.poly in 1:pv){ TimeSpan.begin.l[[i.poly]] <- as.vector(rep(TimeSpan.begin[i.poly], pvn[[i.poly]])) }
            for(i.poly in 1:pv){ TimeSpan.end.l[[i.poly]] <- as.vector(rep(TimeSpan.end[i.poly], pvn[[i.poly]])) }
    }
  }          

  # Polygon styles
  # ==============
  if(!length(unique(colours))==1|colorMode=="normal"){
    colours.l <- list(NULL)
    for(i.poly in 1:pv){ colours.l[[i.poly]] <- as.vector(rep(colours[i.poly], pvn[[i.poly]])) }    
    txts <- sprintf('<Style id="poly%s"><PolyStyle><color>%s</color><outline>%s</outline><fill>%s</fill></PolyStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', 1:sum(unlist(pvn)), unlist(colours.l), rep(as.numeric(outline), sum(unlist(pvn))), as.numeric(!(unlist(hole))))
    parseXMLAndAdd(txts, parent=pl1)
  } else {
    # random colours:
    txts <- sprintf('<Style id="poly%s"><PolyStyle><colorMode>random</colorMode><outline>%s</outline><fill>%s</fill></PolyStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', 1:sum(unlist(pvn)), rep(as.numeric(outline), sum(unlist(pvn))), as.numeric(!(unlist(hole))))
    parseXMLAndAdd(txts, parent=pl1)
  }

  # Point styles
  # ==============
  if(plot.labpt == TRUE){
    sizes.l <- list(NULL)
    shapes.l <- list(NULL)
    # reformat size / shapes:
    for(i.poly in 1:pv){sizes.l[[i.poly]] <- as.vector(rep(sizes[i.poly], pvn[[i.poly]])) }
    for(i.poly in 1:pv){shapes.l[[i.poly]] <- as.vector(rep(shapes[i.poly], pvn[[i.poly]])) }    
    txtsp <- sprintf('<Style id="pnt%s"><LabelStyle><scale>%.1f</scale></LabelStyle><IconStyle><color>ffffffff</color><scale>%s</scale><Icon><href>%s</href></Icon></IconStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', 1:sum(unlist(pvn)), rep(LabelScale, sum(unlist(pvn))), unlist(sizes.l), unlist(shapes.l))
    parseXMLAndAdd(txtsp, parent=pl1b)

  # Writing labpt
  # ================  
  if(all(is.null(unlist(TimeSpan.begin.l))) & all(is.null(unlist(TimeSpan.end.l)))){
    if(all(is.null(unlist(when.l)))){
    # time span undefined:
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></Point></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), rep(as.numeric(extrude), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(labpt)))
    } else {
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></Point></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), unlist(when.l), rep(as.numeric(extrude), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(labpt)))  
    } } else{
      # fixed begin/end times:
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></Point></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), unlist(TimeSpan.begin.l), unlist(TimeSpan.end.l), rep(as.numeric(extrude), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(labpt)))
  }

  parseXMLAndAdd(txtc, parent=pl1b)
  }
  # finished writing the labels

  # Writing polygons
  # ================
  
  if(length(html.table)>0){   
    html.table.l <- list(NULL)
    for(i.poly in 1:pv){ html.table.l[[i.poly]] <- as.vector(rep(html.table[i.poly], pvn[[i.poly]])) }    
  
  # with attributes:
  if(all(is.null(unlist(TimeSpan.begin.l))) & all(is.null(unlist(TimeSpan.end.l)))){
    if(all(is.null(unlist(when.l)))){
    # time span undefined:
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><description><![CDATA[%s]]></description><Polygon><extrude>%.0f</extrude><tessellate>%.0f</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), unlist(html.table.l), rep(as.numeric(extrude), sum(unlist(pvn))), rep(as.numeric(tessellate), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(coords)))
     } else { 
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><description><![CDATA[%s]]></description><Polygon><extrude>%.0f</extrude><tessellate>%.0f</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), unlist(when.l), unlist(html.table.l), rep(as.numeric(extrude), sum(unlist(pvn))), rep(as.numeric(tessellate), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(coords)))
     }} else {
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><description><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><![CDATA[%s]]></description><Polygon><extrude>%.0f</extrude><tessellate>%.0f</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), unlist(TimeSpan.begin.l), unlist(TimeSpan.end.l), unlist(html.table.l), rep(as.numeric(extrude), sum(unlist(pvn))), rep(as.numeric(tessellate), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(coords)))
    }
  }

  # without attributes:
  else{
  if(all(is.null(unlist(TimeSpan.begin.l))) & all(is.null(unlist(TimeSpan.end.l)))){
    if(all(is.null(unlist(when.l)))){
    # time span undefined:
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><Polygon><extrude>%.0f</extrude><tessellate>%.0f</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), rep(as.numeric(extrude), sum(unlist(pvn))), rep(as.numeric(tessellate), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(coords)))
      } else {
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><Polygon><extrude>%.0f</extrude><tessellate>%.0f</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), unlist(when.l), rep(as.numeric(extrude), sum(unlist(pvn))), rep(as.numeric(tessellate), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(coords)))  
     }} else {   
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#poly%s</styleUrl><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><Polygon><extrude>%.0f</extrude><tessellate>%.0f</tessellate><altitudeMode>%s</altitudeMode><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>', unlist(poly_names.l), 1:sum(unlist(pvn)), TimeSpan.begin, TimeSpan.end, rep(as.numeric(extrude), sum(unlist(pvn))), rep(as.numeric(tessellate), sum(unlist(pvn))), rep(altitudeMode, sum(unlist(pvn))), paste(unlist(coords)))     
  }
  }

  parseXMLAndAdd(txt, parent=pl1)

  # save results: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)

}

setMethod("kml_layer", "SpatialPolygons", kml_layer.SpatialPolygons)

# end of script;