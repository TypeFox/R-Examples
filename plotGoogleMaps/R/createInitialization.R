createInitialization <-
function(SP,
                               name="map",
                               zoom=15,
                               mapTypeId = "HYBRID",
                               fitBounds=TRUE,
                               add=FALSE,
                               divname="map_canvas",
                               disableDefaultUI=FALSE,
                               disableDoubleClickZoom =FALSE,
                               draggable= TRUE ,
                               keyboardShortcuts=TRUE,
                               mapTypeControlOptions='DEFAULT',
                               navigationControl=TRUE,
                               navigationControlOptions='DEFAULT',
                               noClear=FALSE,
                               scaleControl=TRUE,
                               scaleControlOptions= 'STANDARD',
                               scrollwheel =TRUE     ,
                               streetViewControl= FALSE) {


                                if (scaleControl!=FALSE)
                                   {scaleControl='true'}else{ scaleControl='false' }

                               if (disableDefaultUI!=FALSE)
                                   {disableDefaultUI='true'}else{ disableDefaultUI='false' }

                               if (disableDoubleClickZoom!=FALSE)
                                   {disableDoubleClickZoom='true'}else{ disableDoubleClickZoom='false' }

                               if (draggable!=FALSE)
                                   {draggable='true'}else{ draggable='false' }

                               if (keyboardShortcuts!=FALSE)
                                   {keyboardShortcuts='true'}else{ keyboardShortcuts='false' }

                               if (navigationControl!=FALSE)
                                   {navigationControl='true'}else{ navigationControl='false' }

                               if (noClear!=FALSE)
                                   {noClear='true'}else{ noClear='false' }

                               if (scaleControl!=FALSE)
                                   {scaleControl='true'}else{ scaleControl='false'
                                    scaleControlOptions='null'}

                               if (scrollwheel!=FALSE)
                                   {scrollwheel='true'}else{ scrollwheel='false' }

                               if (streetViewControl!=FALSE)
                                   {streetViewControl='true'}else{ streetViewControl='false' }


mapTypeId<-paste('google.maps.MapTypeId.',mapTypeId,sep="")
mapTypeControlOptions<-paste('{style: google.maps.MapTypeControlStyle.',
                                      mapTypeControlOptions,'}',sep="")
navigationControlOptions<-paste('{style: google.maps.NavigationControlStyle.',
                                navigationControlOptions ,'}',sep="")
if(!is.null(scaleControlOptions)){
scaleControlOptions<-paste('{style: google.maps.ScaleControlStyle.',
scaleControlOptions ,'}',sep="")  }

SP.ll <- spTransform(SP, CRS("+proj=longlat +datum=WGS84"))
Centar=c(mean(SP.ll@bbox[1,]),mean(SP.ll@bbox[2,]))
sw<-c(SP.ll@bbox[2,1],SP.ll@bbox[1,1])
ne<-c(SP.ll@bbox[2,2],SP.ll@bbox[1,2])



x<-'function initialize() { \n'

x<-paste(x,'var latlng = new google.maps.LatLng(',Centar[2],',',Centar[1],') ; \n')
x<-paste(x,'\n var myOptions = { zoom:',zoom,', \n center: latlng',', \n mapTypeId:',mapTypeId,
' ,\n disableDefaultUI:',disableDefaultUI,' ,\n disableDoubleClickZoom:',disableDoubleClickZoom,
' ,\n  draggable:',draggable,' ,\n  keyboardShortcuts: ', keyboardShortcuts,
' ,\n mapTypeControlOptions:', mapTypeControlOptions,  ' ,\n  navigationControl:',navigationControl,
' ,\n navigationControlOptions:',navigationControlOptions,' ,\n noClear:',noClear,
 ' ,\n scaleControl:', scaleControl ,' ,\n scaleControlOptions:',scaleControlOptions,
' ,\n  scrollwheel:', scrollwheel, ' ,\n streetViewControl:',streetViewControl,'} ; \n')



x<-paste(x,'\n',name)
x<-paste(x,'= new google.maps.Map(document.getElementById("',divname,'"),myOptions); \n',sep="")

if( fitBounds==T){
x<-paste(x,' ',name,'.fitBounds(new google.maps.LatLngBounds(
 new google.maps.LatLng(',sw[1],',',sw[2],'),
 new google.maps.LatLng( ',ne[1],',',ne[2],'))); ',sep="")  }

if (add==F){x<- paste(x,'}')}

return(x)


}
