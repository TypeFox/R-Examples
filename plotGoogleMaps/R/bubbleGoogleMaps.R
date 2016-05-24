bubbleGoogleMaps <-
function(SP,
                     filename="",
                     zcol=1,
                     max.radius=100,
                     key.entries = quantile(SP@data[,zcol],(1:5)/5),
                     do.sqrt = TRUE,
                     add=FALSE,
                     previousMap=NULL,
                     colPalette=NULL,
                     strokeColor="",
                     strokeOpacity=1,
                     fillOpacity=0.7,
                     strokeWeight=1,
                     geodesic=TRUE,
                     clickable=TRUE,
                     zIndex="null",
                     shape="c",
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
                     InfoWindowControl=list(map=map, event="click",position="event.latLng",
                                disableAutoPan=FALSE, maxWidth=330,pixelOffset="null",
                                zIndex="null") ,
                     map="map",
                     mapCanvas="map_canvas",
                     css = "",
                     api="https://maps.google.com/maps/api/js?sensor=false&v=3.18",
                     openMap=TRUE) {


################################################################################
################################################################################
################################################################################

disableDefaultUI=FALSE

if(min(key.entries) <0 ) { key.entries <- sort(c(key.entries, 0))}

	obj = as(SP, "SpatialPointsDataFrame")
	data = obj@data
	if (NCOL(data) == 1){
		z = as.numeric(data[,1])
	         }else {
    	z = as.numeric( data[, zcol] )  }
    	# avoid negative values

kkk=c(min(z),key.entries)

kkk=sapply(2:length(kkk), function(i) mean( c(kkk[i],kkk[i-1]) ) )
    
ke<-abs(kkk) + mean(abs(kkk))    # no zeros as max for radius vecor
# creating a vector for subgroups
if(do.sqrt){scale.level<- sqrt(ke/(max(ke)) ) }else{scale.level<-ke/(max(ke))}
radius.level<-max.radius*scale.level

# new
if (key.entries[1] > min(z) ) {
  breakss <- factor(  c(min(z), key.entries) )
} else {
  breakss <- factor(  c(key.entries) )
}

break_unique <- as.numeric(levels(breakss))

# new
if (break_unique[length(break_unique)] < max(z) ) {
  break_unique[length(break_unique)] <- max(z)
}

#setting of auxiliary variable				
breaks_used=1:(length(break_unique))

if (length(unique(z)) == length(key.entries)) {
  zz = factor(z, labels = radius.level)
  radius.vector <- floor(as.numeric(as.vector(zz)))
}  else {
  break_unique2=break_unique
  break_unique2[length(break_unique2)] <- max(z)+1
  #now find out which categories are used ...
  #findInterval(z, break_unique2)
  #data.frame(z, int=findInterval(z, break_unique2), classes=t(break_unique2))
  breaks_used=sort(unique(findInterval(z, break_unique2)))
  #NEW CODE
  zz = factor(cut(z, break_unique, include.lowest = TRUE), labels = radius.level[breaks_used])        
  radius.vector <- floor(as.numeric(as.vector((zz))))
  
}

               
               
###############


SP.ll <- spTransform(SP, CRS("+proj=longlat +datum=WGS84"))


Centar=c(mean(SP.ll@bbox[1,]),mean(SP.ll@bbox[2,]))
sw<-c(SP.ll@bbox[2,1],SP.ll@bbox[1,1])
ne<-c(SP.ll@bbox[2,2],SP.ll@bbox[1,2])
###################################################

nameOfSP<-sapply(as.list(substitute({SP})[-1]), deparse)
nameOfSP<-gsub('[!,",#,$,%,&,(,),*,+,-,_,.,/,:,;,<,=,>,?,@,^,`,|,~]', "x", nameOfSP)
nameOfSP<-gsub('[[]', "x", nameOfSP)
nameOfSP<-gsub('[]]', "x", nameOfSP)
temporary = FALSE
if(filename==""){
  filename <- tempfile("map", fileext = c(".html"))
  temporary = TRUE
}
randNum = sample(1:10000, 1)
attribute=SP@data[,zcol]
polyName<-paste('poly',nameOfSP,randNum, sep="")
boxname<-paste(polyName,'box',sep="")
textname<- paste(polyName,'text',sep="")
divLegendImage<-tempfile("Legend")  
divLegendImage<-substr(divLegendImage, start=regexpr("Legend",divLegendImage),stop=nchar(divLegendImage))
legendboxname<-paste('box',divLegendImage,sep="")
textnameW<-paste(textname,'W',sep="")
if(layerName==""){
layerName=nameOfSP}



if(strokeColor!=""){
rgbc<-col2rgb(strokeColor)
strokeColor<-rgb(rgbc[1],rgbc[2],rgbc[3],maxColorValue=255) }

if(is.null(colPalette) & min(key.entries)<0){
  colPalette=rep("#99000D",length(key.entries))
  colPalette[which(key.entries<=0)]="#084594"
                                       }

if(!is.null(colPalette)){
rgbc<-col2rgb(colPalette)
colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}

if(length(colPalette)==1){
  colPalette=rep(colPalette,length(key.entries))
  rgbc<-col2rgb(colPalette)
  colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}

if(length(colPalette)==2 & min(key.entries)<0){
  cop=rep(colPalette[2],length(key.entries))
  cop[which(key.entries<0)]=colPalette[1]
  rgbc<-col2rgb(cop)
  colPalette<-apply(rgbc,2,function(x) rgb(x[1],x[2],x[3],maxColorValue=255))}



for(i in 1:length(SP.ll@data)) {
if( identical(attribute,SP.ll@data[,i])){
 attributeName<-names(SP.ll@data)[i]  }
}

att<-rep(NA,length(SP.ll@coords[,1]))
att1=""

if (!is.list(previousMap)) {
  var<-""
  # Declare variables in JavaScript marker and map
  var<- paste(' \n var', map, ' \n')
  # Create all markers and store them in markersArray - PointsName
  }else{ 
    var<-previousMap$var}

var<-paste(var,'var ',polyName,'=[] ; \n')
var1=""


# 
# use only labels for existing intervals
if (length(breaks_used)==length(key.entries)) {
  print("using original PolyCol")
  cxx <- PolyCol(factor(zz, labels = key.entries), colPalette)
} else {
  cxx <- PolyCol(factor(zz, labels = key.entries[breaks_used]), colPalette[breaks_used])
}										   
xx <- cxx$cols

if (length(key.entries) == 1) {
  att_L <- key.entries
} else {
  att_L <- factor(cut(z, break_unique, include.lowest = TRUE))
}



var1<- paste( lapply(as.list(1:length(SP.ll@coords[,1])), function(i) paste(var1, createSphereShape(shape=shape,center=c(SP.ll@coords[i,1],
                                                                                                                        SP.ll@coords[i,2]),
                                                                                                   radius=radius.vector[i],
                                                                                                   fillColor=xx[i],
                                                                                                   map=map,
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
att<- paste ( lapply(as.list(1:length(SP.ll@coords[,1])), function(i) paste(names(SP.ll@data),':',sapply(k ,function(k) as.character(SP.ll@data[i,k]))
                                                                            ,'<br>', collapse="")  )   )



# Put all variables together
var<-paste(var,var1)
if (!is.list(previousMap)) {
  functions<-""
  # Creating functions for checkbox comtrol, Show , Hide and Toggle control
  # Set of JavaScript functionalities
  funs ='function showR(R,boxname, map) {
  R.setMap(map);
  document.getElementById(boxname).checked = true; }
  
  function hideR(R,boxname) {
  R.setMap(null);
  document.getElementById(boxname).checked = false; }
  
  function showO(MLPArray,boxname, map ) { 
  for (var i = 0; i < MLPArray.length; i++) { 
  MLPArray[i].setMap(map); } 
  document.getElementById(boxname).checked = true; }
  
  function hideO(MLPArray,boxname) { 
  for (var i = 0; i < MLPArray.length; i++) { 
  MLPArray[i].setMap(null);} 
  document.getElementById(boxname).checked = false; } 
  
  function boxclick(box,MLPArray,boxname, map) { 
  if (box.checked) { showO(MLPArray,boxname, map); 
  }else {  hideO(MLPArray,boxname);} }
  
  function setOpac(MLPArray,textname){
  opacity=0.01*parseInt(document.getElementById(textname).value) 
  for(var i = 0; i < MLPArray.length; i++) {
  MLPArray[i].setOptions({strokeOpacity: opacity, fillOpacity: opacity}); }}
  
  function setOpacL(MLPArray,textname) {
  opacity=0.01*parseInt(document.getElementById(textname).value) 
  for (var i = 0; i < MLPArray.length; i++) {
  MLPArray[i].setOptions({strokeOpacity: opacity});}}
  
  function setLineWeight(MLPArray,textnameW){
  weight=parseInt(document.getElementById(textnameW).value)
  for (var i = 0; i < MLPArray.length; i++){
  MLPArray[i].setOptions({strokeWeight: weight}); } }
  
  function legendDisplay(box,divLegendImage){
  element = document.getElementById(divLegendImage).style;
  if (box.checked){ element.display="block";} else {  element.display="none";}}
  
  function boxclickR(box,R,boxname, map) {
  if (box.checked){
  showR(R,boxname,map); } else { hideR(R,boxname);} }
  
  function legendDisplay(box,divLegendImage){
  element = document.getElementById(divLegendImage).style; 
  if (box.checked){ element.display="block";} else {  element.display="none";}}  \n '

functions<-paste(functions,funs,sep="")

init<- createInitialization(SP.ll,
                           add=T,
                           name=map,
                           divname=mapCanvas,
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
infW<-""

infW<- paste ( lapply(as.list(1:length(SP.ll@coords[,1])), function(i) 
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


functions<-paste(functions,infW,'showO(',polyName,',"',boxname,'",',map,');',sep="")
fjs=""

fjs<-paste(fjs,'\n USGSOverlay.prototype = new google.maps.OverlayView(); \n',sep="")
fjs<-paste(fjs,'function USGSOverlay(bounds, image, map) {\n      this.bounds_ = bounds;\n      this.image_ = image;\n      this.map_ = map;\n      this.div_ = null;\n      this.setMap(map); }\n',sep="")
fjs<-paste(fjs, 'USGSOverlay.prototype.onAdd = function() {\n      var div = document.createElement("DIV");\n      div.style.border = "none";\n      div.style.borderWidth = "0px";\n      div.style.position = "absolute";\n      var img = document.createElement("img");\n      img.src = this.image_;\n      img.style.width = "100%";\n      img.style.height = "100%";\n      div.appendChild(img);\n      this.div_ = div;\n      this.div_.style.opacity = ',fillOpacity,';\n      var panes = this.getPanes();\n      panes.overlayImage.appendChild(this.div_);}\n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.draw = function() {\n        var overlayProjection = this.getProjection();\n        var sw = overlayProjection.fromLatLngToDivPixel(this.bounds_.getSouthWest());\n        var ne = overlayProjection.fromLatLngToDivPixel(this.bounds_.getNorthEast());\n        var div = this.div_;\n        div.style.left = sw.x + "px";\n        div.style.top = ne.y + "px";\n        div.style.width = (ne.x - sw.x) + "px";\n        div.style.height = (sw.y - ne.y) + "px";} \n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.onRemove = function() { \n this.div_.parentNode.removeChild(this.div_);} \n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.hide = function() { if (this.div_) { this.div_.style.visibility = "hidden";} } \n' ,sep="")
fjs<-paste(fjs,'USGSOverlay.prototype.show = function() {if (this.div_) {  this.div_.style.visibility = "visible";}} \n' ,sep="")
fjs<-paste(fjs,'       USGSOverlay.prototype.toggle = function() { \n if (this.div_) { \n  if (this.div_.style.visibility == "hidden") {  \n   this.show(); \n  } else { \n  this.hide(); } } } \n' ,sep="")
fjs<-paste(fjs,' USGSOverlay.prototype.toggleDOM = function() {\n          if (this.getMap()) {\n            this.setMap(null);\n          } else {\n            this.setMap(this.map_);}}\n' ,sep="")
fjs<-paste(fjs,' function setOpacR(Raster,textname) { \n  opac=0.01*parseInt(document.getElementById(textname).value) \n    Raster.div_.style.opacity= opac } \n' ,sep="")



if(map.width!=control.width & css==""){
  css= paste('\n #',mapCanvas,' { float: left;
 width:', map.width,';
 height:' , map.height,'; }
\n #cBoxes {float: left;
width:', control.width,';
height: ', control.height,';
overflow:auto} \n', sep='') 
}else if (css==""){
  css=paste(' #',mapCanvas,' {min-height: 100%;height:auto; } \n #cBoxes {position:absolute;right:5px; top:50px; background:white}',sep='')
}



starthtm=paste('<!DOCTYPE html> \n <html> \n <head> \n <meta name="viewport" content="initial-scale=1.0, user-scalable=no" />
 <meta charset="utf-8"> \n <style type="text/css">  \n html { height: 100% ; font-size: small} \n body { height: 100%; margin: 0px; padding: 0px }
',css,'
</style> \n
 <script type="text/javascript" src="',api,'"> </script>  \n
 <script language="javascript"> \n ',sep='')

starthtm<-paste(starthtm, fjs)

if (!is.list(previousMap)) {
  endhtm<-paste('</script> \n </head> \n <body onload="initialize()"> \n  <div id="',mapCanvas,'"></div>  \n
                           \n <div id="cBoxes"> \n', sep='') 
} else { endhtm<- previousMap$endhtm }

if(control){
  endhtm<- paste(endhtm,'<table border="0"> \n <tr> \n  <td> 
                                 <input type="checkbox" id="',boxname,'" 
                                 onClick=\'boxclick(this,',polyName,',"',boxname,
                 '",',map,');\' /> <b> ', layerName,'<b> </td> </tr> \n',sep="")
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
  
#   pp<-bubbleLegend(shape=shape,attribute=att_L,colPalette=colPalette
#                    ,legendName=divLegendImage,scale.level=scale.level,strokeColor=strokeColor, temp=temporary)  
pp <- bubbleLegend(shape = shape, attribute = att_L, 
                        colPalette = colPalette, legendName = divLegendImage, 
                        scale.level = scale.level, strokeColor = strokeColor, 
                        temp = temporary , dirname=dirname(filename) ) 

    endhtm<- paste(endhtm,' \n <tr>  <td> <input type="checkbox"  checked="checked" id="'
                   ,legendboxname,'" onClick=\'legendDisplay(this,"',
                   divLegendImage,'");\' /> LEGEND </td> </tr>  <tr> <td>',
                   attributeName,'</td> </tr>
                                    <tr> <td> <div style="display:block;" id="',
                   divLegendImage,'"> <img src="',divLegendImage,
                   '.png" alt="Legend" height="70%"> </div>
                           </td> </tr> \n </table> \n  <hr> \n',sep="") 
}else{ endhtm<- paste(endhtm, '</tr> \n </table> \n <hr>  \n') }




if (add==F){functions<- paste(functions," google.maps.event.addListener( " ,map,", 'rightclick', function(event) {
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
