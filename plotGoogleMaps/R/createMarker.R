createMarker <-
function(lonlat,
                       name="marker",
                       title="Point" ,
                       map="map",
                       clickable=TRUE,
                       draggable=FALSE,
                       flat=TRUE,
                       visible=TRUE,
                       zIndex="null",
                       icon="", ...) {
# ...  shape="" , icon="",shadow="",cursor=""

             p1<-paste('position: new google.maps.LatLng(',lonlat[2],',',lonlat[1],'),\n',sep="")
             m1<-paste('map:',map,',\n',sep="")
             t1<-paste('title:"',title,'",\n',sep="")
             c1<-ifelse(clickable!=FALSE,'clickable: true,\n', 'clickable: false,\n')
             d1<-ifelse(draggable!=FALSE,'draggable: true,\n', 'draggable: false,\n')
             f1<- ifelse(flat!=FALSE,'flat: true,\n', 'flat: false,\n')
             v1<- ifelse(visible!=FALSE,'visible: true,\n', 'visible: false,\n')
             i1<-ifelse (icon==""," " ,paste('icon: ','new google.maps.MarkerImage("',icon,'"), \n',sep=""))


x<-paste('var ',name,'= new google.maps.Marker({ \n',p1,m1,t1,c1,d1,f1,v1,i1,
         ' zIndex:',zIndex)

argsList <- list(...)
vectorNames<-names(argsList)

 if (length(argsList)>1){

     for(i in 1:(length(argsList)-1)) {
     x<-paste(x,', \n',vectorNames[i],':"',argsList[[i]],'"',sep="")
     }
     x<-paste(x,', \n',vectorNames[i+1],':"',argsList[[i+1]],'"}); ',sep="")
  } else { if( length(argsList)==1 ){
          x<-paste(x,', \n',vectorNames[1],':"',argsList[[1]],'"',sep="")}

   x<-paste(x,'}); ') }
   return(x)

}
