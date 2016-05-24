createInfoWindowEventM <-
function(Marker,
                                map="map",
                                event="click",
                                content="The content",
                                position="event.latLng",
                                disableAutoPan=FALSE,
                                maxWidth="null",
                                pixelOffset="null",
                                zIndex="null") {

             if (disableAutoPan!=FALSE)
                  {disableAutoPan='true'}else{ disableAutoPan='false' }


          if (is.numeric(position)) {
          position=paste('position=new google.maps.LatLng(',position[2],',',position[1],')')
          }
           if (is.numeric(position)) {
          pixelOffset=paste('new google.maps.Size',pixelOffset[1],',',pixelOffset[1],')')
          }
x=paste( ' var infowindow = new google.maps.InfoWindow({ content: "", disableAutoPan:', disableAutoPan,
                                                        ', maxWidth :',maxWidth,
                                                        ', pixelOffset :',pixelOffset,
                                                        '  }); ',
         '\n google.maps.event.addListener(',Marker,
',"',event,'",function(event){ 
\n infowindow.setContent("',content,'"); 
\n infowindow.setPosition(event.latLng);
infowindow.open(',map,',' ,Marker,')}); ',sep="")

return(x)

}
