createLine <-
function(oneSPlinesSlot,
                     name="line",
                     map="map",
                     strokeColor="#FFAA00",
                     strokeOpacity=1,
                     strokeWeight=1,
                     geodesic=TRUE,
                     clickable=TRUE,
                     zIndex="null") {

              if (clickable!=FALSE)
                  {clickable='true'}
              else {
                  clickable='false'
             }

              if (geodesic!=FALSE)
                { geodesic='true'}
              else{
                 geodesic='false'
              }

# lonlatmatrix<-slot(slot(oneSPlinesSlot,"Lines")[[1]],"coords")              
              
lonlatmatrix<-  oneSPlinesSlot@Lines[[1]]@coords
              
              
pts=paste('new google.maps.LatLng(',lonlatmatrix[1:((length(lonlatmatrix)/2-1)),2],
          ',',lonlatmatrix[1:((length(lonlatmatrix)/2-1)),1],'),\n',collapse="")
pts=paste('[',pts,'new google.maps.LatLng(',lonlatmatrix[length(lonlatmatrix)/2,2],
          ',',lonlatmatrix[length(lonlatmatrix)/2,1],')]')

x<-paste('var ',name,'= new google.maps.Polyline({ \n path:',pts,', \n','map:',
          map,', \n clickable:',clickable,',\n strokeColor: "',strokeColor,
          '", \n strokeOpacity:',strokeOpacity,',\n strokeWeight:',strokeWeight,
          ',\n geodesic:',geodesic,',\n zIndex:',zIndex,'});')
return(x)

}
