# Purpose        : Writing of Trajectory-type objects to KML
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Pierre Roudier (pierre.roudier@landcare.nz); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : Alpha
# Note           : This method works only with the Space time irregular data frame class objects from the spacetime package; see also how Time Stamps work at [http://kml-samples.googlecode.com]

kml_layer.STTDF <- function(
  obj,
  id.name = names(obj@data)[which(names(obj@data)=="burst")],   # trajectory ID column 
  ## TH: Normally we should be able to pass the ID column via "labels"
  dtime, # time support
  extrude = FALSE,
  # tessellate = FALSE,
  start.icon = paste(get("home_url", envir = plotKML.opts), "3Dballyellow.png", sep=""),
  end.icon = paste(get("home_url", envir = plotKML.opts), "golfhole.png", sep=""),
  LabelScale = .8*get("LabelScale", envir = plotKML.opts),
  z.scale = 1,
  metadata = NULL,
  html.table = NULL,
  ...
  ){
  
  # Get our invisible file connection from custom environment
  kml.out <- get('kml.out', envir=plotKML.fileIO)
  
  # Checking the projection is geo
  prj.check <- check_projection(obj@sp, control = TRUE)

  # Trying to reproject data if the check was not successful
  if (!prj.check) { obj@sp <- reproject(obj@sp) }

  # Parsing the call for aesthetics
  aes <- kml_aes(obj, ...)   

  # Read the relevant aesthetics:
  lines_names <- aes[["labels"]]  
  colours <- aes[["colour"]]   
  width <- aes[["width"]]
  altitudeMode <- aes[["altitudeMode"]]
  balloon <- aes[["balloon"]]

  # object ID names:
  lv <- levels(as.factor(obj@data[,id.name]))
  line.colours <- hex2kml(brewer.pal(n=2+length(lv), name = "Set1"))
  ## names of the coordinate columns:
  nc <- lapply(obj@traj, FUN=function(x){attr(x@sp@coords, "dimname")[[2]]})
  ## strip times:
  xt <- as.POSIXct(unlist(lapply(lapply(obj@traj, slot, "time"), time)), origin="1970-01-01")
  
  ## Format the time slot for writing to KML:
  if(missing(dtime)) {
    when <- format(xt, "%Y-%m-%dT%H:%M:%SZ")
    dtime = 0
  } else {
    ## Begin end times:
    TimeSpan.begin <- format(xt, "%Y-%m-%dT%H:%M:%SZ")
    TimeSpan.end <- format(as.POSIXct(unclass(xt) + dtime, origin="1970-01-01"), "%Y-%m-%dT%H:%M:%SZ")
  }

  # Parse ATTRIBUTE TABLE (for each placemark):
  if (balloon & ("data" %in% slotNames(obj))){
      html.table <- .df2htmltable(obj@data)
  }

  # Name of the object
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  
  # Insert metadata:
  if(!is.null(metadata)){
    md.txt <- kml_metadata(metadata, asText = TRUE)
    txt <- sprintf('<description><![CDATA[%s]]></description>', md.txt)
    parseXMLAndAdd(txt, parent=pl1)
  }
    
  # Sorting lines
  # =============
 
  current.line.coords <- NULL
  ldist <- NULL
  coords <- NULL
  
  for (i.line in 1:length(lv)) {  # for each line

    cfd <- data.frame(coordinates(obj@traj[[i.line]]))
    # convert to line objects (this assumes that the points are sorted chronologically!)
    cl <- Line(cfd)
    # line length:
    ldist[[i.line]] <- LineLength(cl, longlat=TRUE, sum=TRUE) 
    current.line.coords[[i.line]] <- cfd[,nc[[i.line]]]
    if(length(nc[[i.line]])<3){
      current.line.coords[[i.line]][,3] <- rep(0, length(cfd[,1])) 
    } 
    current.line.coords[[i.line]][,3] <- current.line.coords[[i.line]][,3] * z.scale
    # parse coordinates:
    coords[[i.line]] <- paste(current.line.coords[[i.line]][, 1], ',', current.line.coords[[i.line]][, 2], ',', current.line.coords[[i.line]][,3], collapse='\n ', sep = "")
  }

  # Styles - lines:
  # ======
  message("Writing to KML...")
  txts <- sprintf('<Style id="line_%s"><LineStyle><color>%s</color><width>%.1f</width></LineStyle><BalloonStyle><text>$[description]</text></BalloonStyle></Style>', 1:length(lv), line.colours[1:length(lv)], width[1:length(lv)])
  parseXMLAndAdd(txts, parent=pl1)
  
  # Styles - points:
  # ======
  nt <- sapply(current.line.coords, nrow)
  
  nx <- 0
  for (i.line in 1:length(lv)) { 
    # for each line:
    nx <- nx + nt[i.line]
    n1 <- nx - nt[i.line]
    txtsp <- sprintf('<Style id="pnt_%s"><IconStyle><color>%s</color><scale>%.1f</scale><Icon><href>%s</href></Icon></IconStyle></Style>', (n1+1):(nx-1), colours[(n1+1):(nx-1)], rep(LabelScale, nt[i.line]-1), rep(start.icon, nt[i.line]-1))
    parseXMLAndAdd(txtsp, parent=pl1)
    # the last point:
    txtspl <- sprintf('<Style id="pnt_%s"><IconStyle><color>%s</color><scale>%.1f</scale><Icon><href>%s</href></Icon></IconStyle></Style>', nx, colours[nx], LabelScale*2.5, end.icon)
    parseXMLAndAdd(txtspl, parent=pl1)
  }   

  # Writing observed vertices
  # =============

  # for each line:  
  nx <- 0
  
  for (i.line in 1:length(lv)) {
  pl2 = newXMLNode("Folder", parent=pl1)
  pl3 <- newXMLNode("name", lv[i.line], parent = pl2)

  nx <- nx + nt[i.line]
  n1 <- nx - nt[i.line]
  # Parse point coordinates:
  ## TH: I do not know how to make the following code more compact.
  if(length(html.table)>0 & all(dtime==0)){  # with attributes / point temporal support  
    txtl <- sprintf('<Placemark><styleUrl>#pnt_%s</styleUrl><description><![CDATA[%s]]></description><TimeStamp><when>%s</when></TimeStamp><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', (n1+1):nx, html.table[(n1+1):nx], when[(n1+1):nx], rep(as.integer(extrude), nt[i.line]), rep(altitudeMode, nt[i.line]), current.line.coords[[i.line]][,1], current.line.coords[[i.line]][,2], current.line.coords[[i.line]][,3])
  }
  else {
  if(length(html.table)>0 & any(!dtime==0)){  # with attributes / block temporal support 
    txtl <- sprintf('<Placemark><styleUrl>#pnt_%s</styleUrl><description><![CDATA[%s]]></description><TimeStamp><begin>%s</begin><end>%s</end></TimeStamp><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', (n1+1):nx, html.table[(n1+1):nx], TimeSpan.begin[(n1+1):nx], TimeSpan.end[(n1+1):nx], rep(as.integer(extrude), nt[i.line]), rep(altitudeMode, nt[i.line]), current.line.coords[[i.line]][,1], current.line.coords[[i.line]][,2], current.line.coords[[i.line]][,3])
  }
  else {
  if(is.null(html.table) & any(!dtime==0)){   # no attributes / block temporal support 
    txtl <- sprintf('<Placemark><styleUrl>#pnt_%s</styleUrl><TimeStamp><begin>%s</begin><end>%s</end></TimeStamp><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', (n1+1):nx, TimeSpan.begin[(n1+1):nx], TimeSpan.end[(n1+1):nx], rep(as.integer(extrude), nt[i.line]), rep(altitudeMode, nt[i.line]), current.line.coords[[i.line]][,1], current.line.coords[[i.line]][,2], current.line.coords[[i.line]][,3])
  }
  else {  # no attributes / point temporal support 
    txtl <- sprintf('<Placemark><styleUrl>#pnt_%s</styleUrl><TimeStamp><when>%s</when></TimeStamp><Point><extrude>%.0f</extrude><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point></Placemark>', (n1+1):nx, when[(n1+1):nx], rep(as.integer(extrude), nt[i.line]), rep(altitudeMode, nt[i.line]), current.line.coords[[i.line]][,1], current.line.coords[[i.line]][,2], current.line.coords[[i.line]][,3])
  }}}
  
  parseXMLAndAdd(txtl, parent=pl2)
  }

  # Writing Lines
  # =============
  txtl <- sprintf('<Placemark><name>length: %.2f</name><styleUrl>#line_%s</styleUrl><LineString><altitudeMode>%s</altitudeMode><coordinates>%s</coordinates></LineString></Placemark>', unlist(ldist), 1:length(lv), rep(altitudeMode, length(lv)), paste(unlist(coords)))
  
  parseXMLAndAdd(txtl, parent=pl1)
   
  # save results: 
  assign('kml.out', kml.out, envir=plotKML.fileIO)
}

setMethod("kml_layer", "STTDF", kml_layer.STTDF)

# end of script;
