# Purpose        : Parsing SpatialPoints layer to KML
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : Pre-alpha
# Note           : This file gathers the layer() methods. kml.compress(), kml.open() and kml.close();


kml_layer.SpatialPoints <- function(
  obj,
  subfolder.name = paste(class(obj)),
  extrude = TRUE,
  z.scale = 1,
  LabelScale = get("LabelScale", envir = plotKML.opts),
  metadata = NULL,
  html.table = NULL,
  TimeSpan.begin = "",
  TimeSpan.end = "",
  points_names,
  ...
  ){
  
  # invisible file connection
  kml.out <- get('kml.out', envir=plotKML.fileIO)
  
  # Checking the projection
  prj.check <- check_projection(obj, control = TRUE)

  # Trying to reproject data if the check was not successful
  if(prj.check==FALSE) {   obj <- reproject(obj)  }

  # Parsing the call for aesthetics
  aes <- kml_aes(obj, ...)

  # Read the relevant aesthetics
  if(missing(points_names)){ points_names <- aes[["labels"]] }
  colours <- aes[["colour"]]
  shapes <- aes[["shape"]]
  sizes <- aes[["size"]]
  altitude <- aes[["altitude"]]
  altitudeMode <- aes[["altitudeMode"]]
  balloon <- aes[["balloon"]]

  # Parse ATTRIBUTE TABLE (for each placemark):
  if(is.null(html.table)){
    if((balloon==TRUE | class(balloon) %in% c('character','numeric')) & ("data" %in% slotNames(obj))){
      html.table <- .df2htmltable(obj@data) 
    }}

  # Folder and name of the points folder
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", subfolder.name, parent = pl1)

  # Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txt <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txt, parent=pl1)
  }
  message("Writing to KML...")
  
  # Writing points styles
  # =====================
  txts <- sprintf('<Style id="pnt%s"><LabelStyle><scale>%.1f</scale></LabelStyle><IconStyle><color>%s</color><scale>%s</scale><Icon><href>%s</href></Icon></IconStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', 1:length(obj), rep(LabelScale, length(obj)), colours, sizes, shapes)
  parseXMLAndAdd(txts, parent=pl1)
  
  # Writing points coordinates
  # ==========================
  
  # with attributes:
  if(length(html.table)>0){
  if(nzchar(TimeSpan.begin[1])&nzchar(TimeSpan.end[1])){
      if(identical(TimeSpan.begin, TimeSpan.end)){
      when = TimeSpan.begin
      if(length(when)==1){ when = rep(when, length(obj)) }
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><description><![CDATA[%s]]></description><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', points_names, 1:length(obj), when, html.table, rep(as.numeric(extrude), length(obj)), rep(altitudeMode, length(obj)), coordinates(obj)[, 1], coordinates(obj)[, 2], altitude)
      } 
      else{
      if(length(TimeSpan.begin)==1){ TimeSpan.begin = rep(TimeSpan.begin, length(obj)) }
      if(length(TimeSpan.end)==1){ TimeSpan.end = rep(TimeSpan.end, length(obj)) }      
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><description><![CDATA[%s]]></description><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', points_names, 1:length(obj), TimeSpan.begin, TimeSpan.end, html.table, rep(as.numeric(extrude), length(obj)), rep(altitudeMode, length(obj)), coordinates(obj)[, 1], coordinates(obj)[, 2], altitude)      
      }
  }
  else{
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><description><![CDATA[%s]]></description><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', points_names, 1:length(obj), html.table, rep(as.numeric(extrude), length(obj)), rep(altitudeMode, length(obj)), coordinates(obj)[, 1], coordinates(obj)[, 2], altitude)  
  }
  }
  
  # without attributes:
  else{
      if(nzchar(TimeSpan.begin[1])&nzchar(TimeSpan.end[1])){
      if(identical(TimeSpan.begin, TimeSpan.end)){
      when = TimeSpan.begin
      if(length(when)==1){ when = rep(when, length(obj)) }
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', points_names, 1:length(obj),  when, rep(as.numeric(extrude), length(obj)), rep(altitudeMode, length(obj)), coordinates(obj)[, 1], coordinates(obj)[, 2], altitude)
      }
      else {
      if(length(TimeSpan.begin)==1){ TimeSpan.begin = rep(TimeSpan.begin, length(obj)) }
      if(length(TimeSpan.end)==1){ TimeSpan.end = rep(TimeSpan.end, length(obj)) }      
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', points_names, 1:length(obj), TimeSpan.begin, TimeSpan.end, rep(as.numeric(extrude), length(obj)), rep(altitudeMode, length(obj)), coordinates(obj)[, 1], coordinates(obj)[, 2], altitude)    
      }     
  }
      else{
      txtc <- sprintf('<Placemark><name>%s</name><styleUrl>#pnt%s</styleUrl><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', points_names, 1:length(obj), rep(as.numeric(extrude), length(obj)), rep(altitudeMode, length(obj)), coordinates(obj)[, 1], coordinates(obj)[, 2], altitude)      
      }
  }

  parseXMLAndAdd(txtc, parent=pl1)

  # save results: 
  assign('kml.out', kml.out, envir=plotKML.fileIO)

}

setMethod("kml_layer", "SpatialPoints", kml_layer.SpatialPoints)

# end of script;
