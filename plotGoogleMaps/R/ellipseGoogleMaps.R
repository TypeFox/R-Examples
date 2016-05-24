ellipseGoogleMaps <-
function(SP,
                            filename="",
                            zcol=1:3,
                            add=F,
                            previousMap=NULL,
                            scale_e=10,
                            colPalette=NULL,
                            strokeColor="#FFAA00",
                            strokeOpacity=1,
                            fillOpacity=0.7,
                            strokeWeight=1,
                            geodesic=TRUE,
                            clickable=TRUE,
                            zIndex="null",
                            map.width="100%",
                            map.height="100%",
                            layerName="",
                            control.width="100%",
                            control.height="100%",
                            zoom=15,
                            fitBounds=TRUE,
                            mapTypeId = "HYBRID",
                            disableDoubleClickZoom =FALSE,
                            draggable= TRUE ,
                            keyboardShortcuts=TRUE,
                            mapTypeControlOptions='DEFAULT',
                            navigationControl=TRUE,
                            navigationControlOptions='DEFAULT',
                            scaleControlOptions= 'STANDARD',
                            noClear=FALSE,
                            scrollwheel =TRUE     ,
                            streetViewControl= FALSE,
                            legend=TRUE,
                            control=TRUE,
                            InfoWindowControl=list(map="map", event="click",position="event.latLng",
                                disableAutoPan=FALSE, maxWidth=330,pixelOffset="null",
                                zIndex="null") ,
                            api="https://maps.google.com/maps/api/js?sensor=false&v=3.18",
                            openMap=TRUE){
 

    
 ################################ ellipseGoogleMaps ###############################
 ###############################################################################
 ###############################################################################
nameOfSP<-sapply(as.list(substitute({SP})[-1]), deparse)
nameOfSP<-gsub("\\s","", nameOfSP)
nameOfSP<-gsub('[!,",#,$,%,&,(,),*,+,-,_,.,/,:,;,<,=,>,?,@,^,`,|,~]', "x", nameOfSP)
nameOfSP<-gsub('[[]', "X", nameOfSP)
nameOfSP<-gsub('[]]', "X", nameOfSP)


a<- elSPDF(SP,zcol,scale_e)
SP<-a
 
 

   SP.ll <- spTransform(SP, CRS("+proj=longlat +datum=WGS84"))

 
 
disableDefaultUI=FALSE
Centar=c(mean(SP.ll@bbox[1,]),mean(SP.ll@bbox[2,]))
sw<-c(SP.ll@bbox[2,1],SP.ll@bbox[1,1])
ne<-c(SP.ll@bbox[2,2],SP.ll@bbox[1,2])
attribute=SP@data[,zcol[1]]
temporary = FALSE
if(filename==""){
  filename <- tempfile("map", fileext = c(".html"))
  temporary = TRUE
}
randNum = sample(1:10000, 1)

if(layerName==""){
layerName=nameOfSP}

if(strokeColor!=""){
rgbc<-col2rgb(strokeColor)
strokeColor<-rgb(rgbc[1],rgbc[2],rgbc[3],maxColorValue=255) }

if(!is.null(colPalette)){
rgbc<-col2rgb(colPalette)
colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}

if (!is.list(previousMap)) {
functions<-""
# Creating functions for checkbox comtrol, Show , Hide and Toggle control
# Set of JavaScript functionalities
functions<-paste(functions,' function setOpacL(MLPArray,textname) {
opacity=0.01*parseInt(document.getElementById(textname).value) \n
for (var i = 0; i < MLPArray.length; i++) { MLPArray[i].setOptions
({strokeOpacity: opacity}); } } \n',sep="")
functions<-paste(functions,'function showO(MLPArray,boxname) { \n
for (var i = 0; i < MLPArray.length; i++) { \n MLPArray[i].setMap(map); } \n
 document.getElementById(boxname).checked = true; } \n ',sep="")
functions<-paste(functions,'function hideO(MLPArray,boxname) { \n
for (var i = 0; i < MLPArray.length; i++) { \n MLPArray[i].setMap(null);} \n
 document.getElementById(boxname).checked = false; } \n ',sep="")
functions<-paste(functions,'function boxclick(box,MLPArray,boxname)
{ \n if (box.checked) { \n showO(MLPArray,boxname); \n }
 else { \n hideO(MLPArray,boxname);} } \n',sep="")
functions<-paste(functions,' function setOpac(MLPArray,textname)
 {opacity=0.01*parseInt(document.getElementById(textname).value) \n for
  (var i = 0; i < MLPArray.length; i++) { MLPArray[i].setOptions({strokeOpacity:
  opacity, fillOpacity: opacity}); } } \n',sep="")
functions<-paste(functions,' function setLineWeight(MLPArray,textnameW)
 {weight=parseInt(document.getElementById(textnameW).value) \n
 for (var i = 0; i < MLPArray.length; i++)
  { MLPArray[i].setOptions({strokeWeight: weight}); } } \n',sep="")
functions<-paste(functions,'function legendDisplay(box,divLegendImage){ \n
element = document.getElementById(divLegendImage).style; \n if (box.checked)
 { element.display="block";} else {  element.display="none";}} \n',sep="")
functions<-paste(functions,'function showR(R,boxname) { R.setMap(map);
 \n document.getElementById(boxname).checked = true; } \n ',sep="")
functions<-paste(functions,' function hideR(R,boxname) { R.setMap(null);
 \n document.getElementById(boxname).checked = false; } \n ',sep="")
functions<-paste(functions,'function boxclickR(box,R,boxname) { \n if (box.checked)
{ \n showR(R,boxname);\n } else { \n hideR(R,boxname);} }  \n',sep="")
functions<-paste(functions,'function legendDisplay(box,divLegendImage){
\n element = document.getElementById(divLegendImage).style; \n if (box.checked)
 { element.display="block";} else {  element.display="none";}} \n',sep="")
functions<-paste(functions,"google.maps.event.addListener(map,'mousemove',function(event) {
  document.getElementById('latlong').innerHTML = event.latLng.lng() + ', ' + event.latLng.lat()}); \n", sep="")
 
init<-createInitialization(SP.ll,
                               add=T,
                               zoom=zoom,
                               fitBounds=fitBounds,
                               mapTypeId = mapTypeId,
                               disableDefaultUI=disableDefaultUI,
                               disableDoubleClickZoom =disableDoubleClickZoom,
                               draggable= draggable ,
                               keyboardShortcuts=keyboardShortcuts,
                               mapTypeControlOptions=mapTypeControlOptions,
                               scaleControlOptions=scaleControlOptions,
                               navigationControl=navigationControl,
                               navigationControlOptions=navigationControlOptions,
                               noClear=noClear,
                               scrollwheel =scrollwheel    ,
                               streetViewControl= streetViewControl)
# Put all functions together
functions<-paste( functions,init, sep="")  }else{ functions<- previousMap$functions}

fjs=""

fjs<-paste(fjs,'\n USGSOverlay.prototype = new google.maps.OverlayView(); \n',sep="")
fjs<-paste(fjs,'function USGSOverlay(bounds, image, map) {\n      this.bounds_ = bounds;\n      this.image_ = image;\n      this.map_ = map;\n      this.div_ = null;\n      this.setMap(map); }\n',sep="")
fjs<-paste(fjs, 'USGSOverlay.prototype.onAdd = function() {\n      var div = document.createElement("DIV");\n      div.style.border = "none";\n      div.style.borderWidth = "0px";\n      div.style.position = "absolute";\n      var img = document.createElement("img");\n      img.src = this.image_;\n      img.style.width = "100%";\n      img.style.height = "100%";\n      div.appendChild(img);\n      this.div_ = div;\n      this.div_.style.opacity = 0.5;\n      var panes = this.getPanes();\n      panes.overlayImage.appendChild(this.div_);}\n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.draw = function() {\n        var overlayProjection = this.getProjection();\n        var sw = overlayProjection.fromLatLngToDivPixel(this.bounds_.getSouthWest());\n        var ne = overlayProjection.fromLatLngToDivPixel(this.bounds_.getNorthEast());\n        var div = this.div_;\n        div.style.left = sw.x + "px";\n        div.style.top = ne.y + "px";\n        div.style.width = (ne.x - sw.x) + "px";\n        div.style.height = (sw.y - ne.y) + "px";} \n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.onRemove = function() { \n this.div_.parentNode.removeChild(this.div_);} \n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.hide = function() { if (this.div_) { this.div_.style.visibility = "hidden";} } \n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.show = function() {if (this.div_) {  this.div_.style.visibility = "visible";}} \n' ,sep="")
fjs<-paste(fjs,'       USGSOverlay.prototype.toggle = function() { \n if (this.div_) { \n  if (this.div_.style.visibility == "hidden") {  \n   this.show(); \n  } else { \n  this.hide(); } } } \n' ,sep="")
fjs<-paste(fjs,' USGSOverlay.prototype.toggleDOM = function() {\n          if (this.getMap()) {\n            this.setMap(null);\n          } else {\n            this.setMap(this.map_);}}\n' ,sep="")
fjs<-paste(fjs,' function setOpacR(Raster,textname) { \n  opac=0.01*parseInt(document.getElementById(textname).value) \n    Raster.div_.style.opacity= opac } \n' ,sep="")

 if(map.width!=control.width){
   css= paste('\n #map_canvas { float: left;
 width:', map.width,';
              height:' , map.height,'; }
              \n #cBoxes {float: left;
              width:', control.width,';
              height: ', control.height,';
              overflow:auto} \n') 
}else{
  css=' #map_canvas {min-height: 100%;
height:auto; }
  
  #cBoxes {position:absolute;
  right:5px;
  top:50px;
  background:white}'
}
 
 
 starthtm=paste('<html> \n <head> \n <meta name="viewport" content="initial-scale=1.0, user-scalable=no" />
\n <style type="text/css">  \n html { height: 100% } \n body { height: 100%; margin: 0px; padding: 0px }
                ',css,'
                </style> \n
                <script type="text/javascript" src="http://maps.google.com/maps/api/js?sensor=false"> </script>  \n
                <script language="javascript"> \n ')
 
starthtm<-paste(starthtm, fjs)

################################################################################
################################################################################



                  polyName<-paste('poly',nameOfSP,randNum, sep="")
                  boxname<-paste(nameOfSP,'box',sep="")
                  textname<- paste(nameOfSP,'text',sep="")
                  divLegendImage<-tempfile("Legend")  
                  divLegendImage<-substr(divLegendImage, start=regexpr("Legend",divLegendImage),stop=nchar(divLegendImage))
                  legendboxname<-paste('box',divLegendImage,sep="")
                  textnameW<-paste(textname,'W',sep="")
                  
                  
                  
                  for(i in 1:length(SP.ll@data)) {
                  if( identical(attribute,SP.ll@data[,i])){
                   attributeName<-names(SP.ll@data)[i]  }
                  }
                  
                  att<-rep(NA,length(slot(SP.ll,"polygons")))
                  att1=""
                  
                  if (!is.list(previousMap)) {
                  var<-""
                  # Declare variables in JavaScript marker and map
                  var<-c(' \n var map \n')
                  # Create all markers and store them in markersArray - PointsName
                  }else{ var<-previousMap$var}
                  
                  var<-paste(var,'var ',polyName,'=[] ; \n')
                  var1=""
                  cxx<-PolyCol(attribute,colPalette)
                   xx<-cxx$cols
                  
#  pp<-legendbar(cxx$brks,colPalette=cxx$col.uniq,legendName=divLegendImage)

 var1<- paste( lapply(as.list(1:length(SP.ll@polygons)), function(i) paste(var1,createPolygon(SP.ll@polygons[[i]],
                                                                                              fillColor=xx[i],
                                                                                              strokeColor=strokeColor,
                                                                                              strokeOpacity=strokeOpacity,
                                                                                              strokeWeight=strokeWeight,
                                                                                              geodesic=geodesic,
                                                                                              clickable=clickable,
                                                                                              fillOpacity=fillOpacity,
                                                                                              zIndex=zIndex),'\n',sep="") 
 )
               ,polyName,'.push(polygon); \n',sep="",collapse='\n')

 k = 1:length(names(SP.ll@data))
 att<- paste ( lapply(as.list(1:length(SP.ll@polygons)), function(i) paste(names(SP.ll@data),':',sapply(k ,function(k) as.character(SP.ll@data[i,k]))
                                                                           ,'<br>', collapse="")  )   )
 
 
                  
                  # Put all variables together
                  var<-paste(var,var1)
                  
                  infW<-""
 
infW<- paste ( lapply(as.list(1:length(SP.ll@polygons)), function(i) 
  paste(infW,createInfoWindowEvent(Line_or_Polygon=
                                     paste(polyName,'[',i-1,'] ',sep=""),
                                   content=att[i],
                                   map=InfoWindowControl$map,
                                   event=InfoWindowControl$event,
                                   position= InfoWindowControl$position,
                                   disableAutoPan = InfoWindowControl$disableAutoPan,
                                   maxWidth=InfoWindowControl$maxWidth,
                                   pixelOffset=InfoWindowControl$pixelOffset,
                                   zIndex=InfoWindowControl$zIndex)
        ,' \n')  )  ,collapse='\n' )                
 
                  
#                   for(i in 1:length(SP.ll@polygons)){
#                   infW<-paste(infW,createInfoWindowEvent(Line_or_Polygon=
#                       paste(polyName,'[',i-1,'] ',sep=""),content=att[i]),' \n')}
                  
                  functions<-paste(functions,infW,'showO(',polyName,',"',boxname,'");',sep="")
                  
                  
                  
                  if (!is.list(previousMap)) {
                    endhtm<-c('</script> \n </head> \n <body onload="initialize()"> \n 
                            <div id="map_canvas"></div>  \n
                           \n <div id="cBoxes"> \n')
                  } else { endhtm<- previousMap$endhtm }
                  
if(control){
  endhtm<- paste(endhtm,'<table border="0"> \n <tr> \n  <td> 
                                 <input type="checkbox" id="',boxname,'" 
                 onClick=\'boxclick(this,',polyName,',"',boxname,
                 '");\' /> <b> ', layerName,'<b> </td> </tr> \n',sep="")
  endhtm<- paste(endhtm,' \n <tr> <td> \n <input type="text" id="',
                 textname,'" value="50" onChange=\'setOpac(',
                 polyName,',"',textname,'")\' size=3 /> 
                 Opacity (0-100 %) </td> </tr> \n',sep="")
  endhtm<- paste(endhtm,' \n <tr>  <td> \n <input type="text" 
                 id="',textnameW,'" value="1" onChange=\'
                 setLineWeight(',polyName,',"',textnameW,'")\' 
                 size=3 /> Line weight (pixels) </td> </tr> \n ',sep="")
}
if(legend){
  
  pp <-legendbar(cxx$brks,colPalette=cxx$col.uniq,legendName=divLegendImage, temp=temporary)  
  
  endhtm<- paste(endhtm,' \n <tr>  <td> <input type="checkbox"  checked="checked" id="'
                 ,legendboxname,'" onClick=\'legendDisplay(this,"',
                 divLegendImage,'");\' /> LEGEND </td> </tr>  <tr> <td>',
                 attributeName,'</td> </tr>
                 <tr> <td> <div style="display:block;" id="',
                 divLegendImage,'"> <img src="',divLegendImage,
                 '.png" alt="Legend" height="70%"> </div>
                 </td> </tr> \n </table> \n  <hr> \n',sep="") 
}else{ endhtm<- paste(endhtm, '</tr> \n </table> \n <hr>  \n') }

 if (add==F){functions<- paste(functions," google.maps.event.addListener(map, 'rightclick', function(event) {
    var lat = event.latLng.lat();
    var lng = event.latLng.lng();
    alert('Lat=' + lat + '; Lng=' + lng);}); " , " \n }" )
endhtm<-paste(endhtm,'</div> \n </body>  \n  </html>')
write(starthtm, filename,append=F)
write(var, filename,append=TRUE)
write(functions, filename,append=TRUE)
write(endhtm, filename,append=TRUE)
if(openMap){browseURL(filename)}
             }


x <- list(starthtm=starthtm,var=var, functions=functions,endhtm=endhtm)
return(x)

}
