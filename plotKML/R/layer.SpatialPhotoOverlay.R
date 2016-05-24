# Purpose        : Parsing "SoilPhotoOverlay" objects to KML
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ;
# Status         : alpha version
# Note           : plots either a PhotoOverlay or a monolith (COLLADA model);


kml_layer.SpatialPhotoOverlay <- function(
  obj,
  method = c("PhotoOverlay", "monolith")[1],
  PhotoOverlay.shape = obj@PhotoOverlay$shape,
  href = obj@filename,
  coords,
  dae.name = "",
  heading = obj@PhotoOverlay$heading,
  tilt = obj@PhotoOverlay$tilt,
  roll = obj@PhotoOverlay$roll,    
  near = obj@PhotoOverlay$near,
  range = obj@PhotoOverlay$range,
  leftFov = obj@PhotoOverlay$leftFov, 
  rightFov = obj@PhotoOverlay$rightFov, 
  bottomFov = obj@PhotoOverlay$bottomFov, 
  topFov = obj@PhotoOverlay$topFov, 
  altitudeMode = "clampToGround",
  block.size = 100, # about 3-arcsecs
  max.depth = 300,
  scale.x = 1,
  scale.y = 1,
  scale.z = 1,
  refreshMode = "once",
  html.table = NULL, 
  ...
  ){

  # get our invisible file connection from custom evnrionment
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  # check the projection:
  prj.check <- check_projection(obj@sp, control = TRUE) 

  # Trying to reproject data if the check was not successful
  if (!prj.check) { obj@sp <- reproject(obj@sp) }  
  
  when <- obj@exif.info[["DateTime"]]
  LON <- as.vector(coordinates(obj@sp)[,1])
  LAT <- as.vector(coordinates(obj@sp)[,2])
  ALT <- as.vector(coordinates(obj@sp)[,3])
  if(ALT == 0|is.na(ALT)) { ALT = 10 }
    else {
    altitudeMode = "absolute"
    } 

  # unique image ID:
  x <- strsplit(obj@filename, "/")[[1]]
  image.id <- x[length(x)]

  # Parsing the call for aesthetics
  aes <- kml_aes(obj, ...)

  # Read the relevant aesthetics
  # altitudeMode <- aes[["altitudeMode"]]
  balloon <- aes[["balloon"]]

  # Parse ATTRIBUTE TABLE (EXIF data):
  if ((is.logical(balloon) | class(balloon) %in% c('character','numeric')) & ("exif.info" %in% slotNames(obj))){
     html.table <- .df2htmltable(data.frame(obj@exif.info)) 
  }

  # Simple PhotoOverlay
  if(method == "PhotoOverlay"){
  
  # Object name:
  pl1 = newXMLNode("PhotoOverlay", attrs = c("id" = image.id), parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", image.id, parent = pl1)
  pl2b <- newXMLNode("TimeStamp", parent = pl1)
  pl3b <- newXMLNode("when", when, parent = pl2b)
  # description:
  txtd <- sprintf('<description><![CDATA[%s]]></description>', html.table) 
  parseXMLAndAdd(txtd, parent=pl1) 

  # Add a Camera view location:
  # ==========================    
  txtv <- sprintf('<Camera><longitude>%.6f</longitude><latitude>%.6f</latitude><altitude>%f</altitude><heading>%.1f</heading><tilt>%.1f</tilt><roll>%.1f</roll></Camera>', LON, LAT, ALT, heading, tilt, roll) 
  parseXMLAndAdd(txtv, parent=pl1)
  
  pl3 <- newXMLNode("Icon", parent = pl1)
  pl4 <- newXMLNode("href", href, parent = pl3)  

  # Add ViewVolume info:
  # ==========================    
  txtw <- sprintf('<ViewVolume><leftFov>%s</leftFov><rightFov>%s</rightFov><bottomFov>%s</bottomFov><topFov>%s</topFov><near>%s</near></ViewVolume>', leftFov, rightFov, bottomFov, topFov, near) 
  parseXMLAndAdd(txtw, parent=pl1)

  # Location of the camera
  # ========================== 
  txtc <- sprintf('<Point><altitudeMode>%s</altitudeMode><coordinates>%.5f,%.5f,%.0f</coordinates></Point>', altitudeMode, LAT, LON, ALT)   
  parseXMLAndAdd(txtc, parent=pl1)
  
  # Display type:
  pl5 <- newXMLNode("shape", PhotoOverlay.shape, parent = pl1)  
  }

  # Block model
  if(method == "monolith"){

  # derive required input pars: 
  asp = obj@exif.info$ImageWidth / obj@exif.info$ImageHeight
  block.height = block.size / asp
  if(dae.name==""){
    dae.name = gsub(x=image.id, "\\.", "_")
  }

  # make a COLLADA file:
  if(missing(coords)){
  X1 = block.size/2; X2 = -block.size/2; X3 = block.size/2; X4 = -block.size/2
  Y1 = 0; Y2 = 0; Y3 = 0; Y4 = 0
  ALT1 = max.depth; ALT2 = max.depth; ALT3 = max.depth-block.height; ALT4 = max.depth-block.height 
  coords <- matrix(c(X=c(X1, X2, X3, X4), Z=c(ALT1, ALT2, ALT3, ALT4), Y=c(Y1, Y2, Y3, Y4)), ncol=3)
  }
  
  makeCOLLADA.rectangle(filename = set.file.extension(dae.name, ".dae"), coords = coords, href = href)
  
  # Object name:
  pl1 = newXMLNode("Folder", parent=kml.out[["Document"]])
  pl2 <- newXMLNode("name", image.id, parent = pl1)
  # description:
  txtd <- sprintf('<description><![CDATA[%s]]></description>', html.table) 
  parseXMLAndAdd(txtd, parent=pl1)  
  txtl <- sprintf('<LookAt><longitude>%.5f</longitude><latitude>%.5f</latitude><altitude>%.1f</altitude><heading>%.1f</heading><tilt>%.1f</tilt><range>%.5f</range><altitudeMode>relativeToGround</altitudeMode></LookAt>', LON, LAT, max.depth/2, heading, tilt, range) 
  parseXMLAndAdd(txtl, parent=pl1)

  # Parse the placemark:
  txtm <- sprintf('<Placemark><name>monolith</name><TimeStamp><when>%s</when></TimeStamp><Model id="%s"><altitudeMode>relativeToGround</altitudeMode><Location><longitude>%.5f</longitude><latitude>%.5f</latitude><altitude>%.0f</altitude></Location><Orientation><heading>%.1f</heading><tilt>%.1f</tilt><roll>%.1f</roll></Orientation><Scale><x>%.1f</x><y>%.1f</y><z>%.1f</z></Scale><Link><href>%s</href><refreshMode>%s</refreshMode></Link><ResourceMap><Alias><targetHref>%s</targetHref><sourceHref>%s</sourceHref></Alias></ResourceMap></Model></Placemark>', when, dae.name, LON, LAT, max.depth, heading, tilt, roll, scale.x, scale.y, scale.z, set.file.extension(dae.name, ".dae"), refreshMode, href, href)
  parseXMLAndAdd(txtm, parent=pl1)
  }

  # save results: 
  assign("kml.out", kml.out, envir=plotKML.fileIO)
}
    
setMethod("kml_layer", "SpatialPhotoOverlay", kml_layer.SpatialPhotoOverlay)


## Create a block model:
makeCOLLADA.rectangle <- function(coords, filename, href, DateTime, up_axis = "Z_UP", authoring_tool = "plotKML", technique_profile = "GOOGLEEARTH", double_sided = TRUE){

  # coordinates of the bbox:
  pnts = paste(signif(as.vector(t(coords)), 5), collapse=" ")
  if(missing(DateTime)) { DateTime <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ") }
  
  # file header:
  doc = newXMLDoc()
  dae.out <- newXMLNode("COLLADA", attrs=c(version="1.5.0"), namespaceDefinitions = c("xmlns"="http://www.collada.org/2008/03/COLLADASchema"), parent=doc)
  txta <- sprintf('<asset><contributor><authoring_tool>%s</authoring_tool></contributor><created>%s</created><modified>%s</modified><unit meter="1" name="meter" /><up_axis>%s</up_axis></asset>', authoring_tool, DateTime, DateTime, up_axis)
  parseXMLAndAdd(txta, parent=dae.out)

  # library visual scenes:
  txts <- sprintf('<library_visual_scenes><visual_scene id="ID1"><node name="rectangle">         <instance_geometry url="#ID2"><bind_material><technique_common><instance_material symbol="Material2" target="#ID3"><bind_vertex_input semantic="UVSET0" input_semantic="TEXCOORD" input_set="0" /></instance_material></technique_common></bind_material></instance_geometry></node></visual_scene></library_visual_scenes>', authoring_tool, DateTime, DateTime, up_axis)    
  parseXMLAndAdd(txts, parent=dae.out)
     
  # library geometries:
  txtg <- sprintf('<library_geometries><geometry id="ID2"><mesh><source id="ID8"><float_array id="ID12" count="12">%s</float_array><technique_common><accessor count="4" source="#ID12" stride="3"><param name="X" type="float" /><param name="Y" type="float" /><param name="Z" type="float" /></accessor></technique_common></source><source id="ID11"><float_array id="ID14" count="8">0 0 -1 0 0 1 -1 1</float_array><technique_common><accessor count="4" source="#ID14" stride="2"><param name="S" type="float" /><param name="T" type="float" /></accessor></technique_common></source><vertices id="ID10"><input semantic="POSITION" source="#ID8" /><input semantic="NORMAL" source="#ID9" /></vertices><triangles count="2" material="Material2"><input offset="0" semantic="VERTEX" source="#ID10" /><input offset="1" semantic="TEXCOORD" source="#ID11" /><p>0 0 1 1 2 2 3 3 2 2 1 1</p></triangles></mesh></geometry></library_geometries>', pnts)    
  parseXMLAndAdd(txtg, parent=dae.out)  

  # library materials:
  txtm <- sprintf('<library_materials><material id="ID3" name="material_0"><instance_effect url="#ID4" /></material></library_materials>')    
  parseXMLAndAdd(txtm, parent=dae.out)
    
  # library effects:
  txte <- sprintf('<library_effects><effect id="ID4"><profile_COMMON><newparam sid="ID6"><surface type="2D"><init_from>ID5</init_from></surface></newparam><newparam sid="ID7"><sampler2D><source>ID6</source></sampler2D></newparam><technique sid="COMMON"><lambert><diffuse><texture texture="ID7" texcoord="UVSET0" /></diffuse></lambert></technique><extra><technique profile="%s"><double_sided>%s</double_sided></technique></extra></profile_COMMON></effect></library_effects>', technique_profile, as.numeric(double_sided))
  parseXMLAndAdd(txte, parent=dae.out)  
  
  # library images:
  txti <- sprintf('<library_images><image id="ID5"><init_from>%s</init_from></image></library_images><scene><instance_visual_scene url="#ID1" /></scene> ', href)
  parseXMLAndAdd(txti, parent=dae.out)

  saveXML(doc, filename)
  message(paste("Creating", filename))  

}

# end of script;