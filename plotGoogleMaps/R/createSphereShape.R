createSphereShape <-
function( shape="t",
                             center,
                             radius=20,  #m
                             name="polygon",
                             fillColor="#00AAFF",
                             fillOpacity=0.7,
                             map="map",
                             strokeColor="#FDA0C0C",
                             strokeOpacity=1,
                             strokeWeight=1,
                             geodesic="null",
                             clickable=TRUE,
                             zIndex=1) {

if (shape=="t"){
x<-createSphereTriangle(center=center,
                             radius=radius,  #m
                             name=name,
                             fillColor=fillColor,
                             fillOpacity=fillOpacity,
                             map=map,
                             strokeColor=strokeColor,
                             strokeOpacity=strokeOpacity,
                             strokeWeight=strokeWeight,
                             geodesic=geodesic,
                             clickable=clickable,
                             zIndex=zIndex)
        }else if (shape=="q"){
          x<-createSphereQuadrangle(center=center,
                             radius=radius,  #m
                             name=name,
                             fillColor=fillColor,
                             fillOpacity=fillOpacity,
                             map=map,
                             strokeColor=strokeColor,
                             strokeOpacity=strokeOpacity,
                             strokeWeight=strokeWeight,
                             geodesic=geodesic,
                             clickable=clickable,
                             zIndex=zIndex)

        }else{
                  x<-createSphereCircle(center=center,
                             radius=radius,  #m
                             name=name,
                             fillColor=fillColor,
                             fillOpacity=fillOpacity,
                             map=map,
                             strokeColor=strokeColor,
                             strokeOpacity=strokeOpacity,
                             strokeWeight=strokeWeight,
                             geodesic=geodesic,
                             clickable=clickable,
                             zIndex=zIndex)


      }

      return(x)

}
