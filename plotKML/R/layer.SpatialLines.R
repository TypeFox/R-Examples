# Purpose        : Write line-objects (SpatialLines) to KML;
# Maintainer     : Pierre Roudier (pierre.roudier@landcare.nz);
# Contributions  : Tomislav Hengl (tom.hengl@wur.nl); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : only 2-3 aesthetics can be used - color, width and names;

kml_layer.SpatialLines <- function(
  obj,
  subfolder.name = paste(class(obj)),
  extrude = FALSE,
  z.scale = 1,
  metadata = NULL,
  html.table = NULL,
  TimeSpan.begin = "",
  TimeSpan.end = "",
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
  lines_names <- aes[["labels"]]
  colours <- aes[["colour"]]
  widths <- aes[["width"]]
  altitude <- aes[["altitude"]]
  altitudeMode <- aes[["altitudeMode"]]
  balloon <- aes[["balloon"]]

  # Parse ATTRIBUTE TABLE (for each placemark):
  if(balloon & ("data" %in% slotNames(obj))){
      html.table <- .df2htmltable(obj@data)
  }
  
  message("Writing to KML...")
  # Folder name / name of the points folder
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", subfolder.name, parent = pl1)

  # Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txtm <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txtm, parent=pl1)
  }  

  # process lines:
  lv <- length(obj@lines) 
  coords <- NULL
  for (i.line in 1:lv) { 
    xyz <- slot(slot(obj@lines[[i.line]], "Lines")[[1]], "coords")
      if(ncol(xyz)==2){ xyz <- cbind(xyz, rep(altitude[i.line], nrow(xyz))) }
    coords[[i.line]] <- paste(xyz[,1], ',', xyz[,2], ',', xyz[,3], collapse='\n ', sep = "")
  }
  
  # Line styles
  # ============
  txts <- sprintf('<Style id="line%s"><LineStyle><color>%s</color><width>%.0f</width></LineStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', 1:lv, colours, widths)
  parseXMLAndAdd(txts, parent=pl1)  
  
  # Writing lines
  # =============

  # with attributes:
  if(length(html.table)>0){
  if(nzchar(TimeSpan.begin[1])&nzchar(TimeSpan.end[1])){
      if(identical(TimeSpan.begin, TimeSpan.end)){
      when = TimeSpan.begin
      if(length(when)==1){ when = rep(when, lv) }
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#line%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><description><![CDATA[%s]]></description><LineString><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', lines_names, 1:lv, when, html.table, rep(as.numeric(extrude), lv), rep(altitudeMode, lv), paste(unlist(coords)))
      } 
      else{
      if(length(TimeSpan.begin)==1){ TimeSpan.begin = rep(TimeSpan.begin, lv) }
      if(length(TimeSpan.end)==1){ TimeSpan.end = rep(TimeSpan.end, lv) }
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#line%s</styleUrl><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><description><![CDATA[%s]]></description><LineString><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', lines_names, 1:lv, TimeSpan.begin, TimeSpan.end, html.table, rep(as.numeric(extrude), lv), rep(altitudeMode, lv), paste(unlist(coords)))    
      }
  }
  else{      
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#line%s</styleUrl><description><![CDATA[%s]]></description><LineString><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', lines_names, 1:lv, html.table, rep(as.numeric(extrude), lv), rep(altitudeMode, lv), paste(unlist(coords)))
  }
  }
  
  # without attributes:
  else{
      if(nzchar(TimeSpan.begin[1])&nzchar(TimeSpan.end[1])){
      if(identical(TimeSpan.begin, TimeSpan.end)){
      when = TimeSpan.begin
      if(length(when)==1){ when = rep(when, lv) }
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#line%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><LineString><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', lines_names, 1:lv, when, rep(as.numeric(extrude), lv), rep(altitudeMode, lv), paste(unlist(coords)))
     }
      else {
      if(length(TimeSpan.begin)==1){ TimeSpan.begin = rep(TimeSpan.begin, lv) }
      if(length(TimeSpan.end)==1){ TimeSpan.end = rep(TimeSpan.end, lv) }   
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#line%s</styleUrl><TimeSpan><begin>%s</begin><end>%s</end></TimeSpan><LineString><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', lines_names, 1:lv, TimeSpan.begin, TimeSpan.end, rep(as.numeric(extrude), lv), rep(altitudeMode, lv), paste(unlist(coords)))
      }     
  }
      else{
      txt <- sprintf('<Placemark><name>%s</name><styleUrl>#line%s</styleUrl><LineString><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', lines_names, 1:lv, rep(as.numeric(extrude), lv), rep(altitudeMode, lv), paste(unlist(coords))) 
      }
  }

  parseXMLAndAdd(txt, parent=pl1)

  # save results: 
  assign('kml.out', kml.out, envir=plotKML.fileIO)

}

setMethod("kml_layer", "SpatialLines", kml_layer.SpatialLines)

# end of script;