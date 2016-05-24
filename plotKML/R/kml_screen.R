# Purpose        : Adds different types of screen overlays to KML;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Pierre Roudier (pierre.roudier@landcare.nz); Dylan Beaudette (debeaudette@ucdavis.edu); 
# Status         : pre-alpha
# Note           : for more info see [http://code.google.com/apis/kml/documentation/kmlreference.html#screenoverlay];
 
 
kml_screen <- function(
  image.file,
  sname = "", 
  position = c("UL","ML","LL","BC","LR","MR","UR","TC")[1],
  overlayXY,
  screenXY,
  xyunits = c("fraction", "pixels", "insetPixels")[1],
  rotation = 0,
  size = c(0,0)
  ){

  # get our invisible file connection from custom evnrionment
  kml.out <- get("kml.out", envir=plotKML.fileIO)
  
  # if nothing is specified estimate position of the overlay based on "position":
  if(position=="UL"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="0" y="1"'; screenXY = 'x="0" y="1"' }
  if(position=="ML"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="0" y="0.5"'; screenXY = 'x="0" y="0.5"' }
  if(position=="LL"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="0" y="0"'; screenXY = 'x="0" y="0"' }
  if(position=="BC"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="0.5" y="0"'; screenXY = 'x="0.5" y="0"' }
  if(position=="LR"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="1" y="0"'; screenXY = 'x="1" y="0"' }
  if(position=="MR"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="1" y="0.5"'; screenXY = 'x="1" y="0.5"' }
  if(position=="UR"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="1" y="1"'; screenXY = 'x="1" y="1"' }
  if(position=="TC"&missing(screenXY)&missing(overlayXY)){ overlayXY = 'x="0.5" y="1"'; screenXY = 'x="0.5" y="1"' }
  
  # Parse KML code:
  if(size[1] == 0&size[2] == 0&rotation == 0){
  screen_txt <- sprintf('<ScreenOverlay><name>%s</name><Icon><href>%s</href></Icon><overlayXY %s xunits="%s" yunits="%s"/><screenXY %s xunits="%s" yunits="%s"/></ScreenOverlay>', sname, image.file, overlayXY, xyunits, xyunits, screenXY, xyunits, xyunits)
  }
  else { if(size[1] == 0&size[2] == 0&rotation != 0){
  screen_txt <- sprintf('<ScreenOverlay><name>%s</name><Icon><href>%s</href></Icon><overlayXY %s xunits="%s" yunits="%s"/><screenXY %s xunits="%s" yunits="%s"/><rotation>%.1f</rotation></ScreenOverlay>', sname, image.file, overlayXY, xyunits, xyunits, screenXY, xyunits, xyunits, rotation)
  }
  else { if(size[1] != 0&size[2] != 0&rotation == 0){
  screen_txt <- sprintf('<ScreenOverlay><name>%s</name><Icon><href>%s</href></Icon><overlayXY %s xunits="%s" yunits="%s"/><screenXY %s xunits="%s" yunits="%s"/><size x="%s" y="%s" xunits="%s" yunits="%s"/></ScreenOverlay>', sname, image.file, overlayXY, xyunits, xyunits, screenXY, xyunits, xyunits, size[1], size[2], xyunits, xyunits)
  }
  else {
  screen_txt <- sprintf('<ScreenOverlay><name>%s</name><Icon><href>%s</href></Icon><overlayXY %s xunits="%s" yunits="%s"/><screenXY %s xunits="%s" yunits="%s"/><size x="%s" y="%s" xunits="%s" yunits="%s"/></ScreenOverlay><rotation>%.1f</rotation>', sname, image.file, overlayXY, xyunits, xyunits, screenXY, xyunits, xyunits, size[1], size[2], xyunits, xyunits, rotation)
  }}}
  
  parseXMLAndAdd(screen_txt, parent=kml.out[["Document"]])
  
  # save to a connection: 
  assign('kml.out', kml.out, envir=plotKML.fileIO)
  
}