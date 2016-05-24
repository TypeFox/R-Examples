PolarPlot <-
function(theta,r,degrees=FALSE,method=plot,geographical=FALSE,...)
  {
    if(degrees){theta=theta*pi/180}
    if(geographical){
      x = r * sin(theta)
      y = r * cos(theta)
    }else{
      x = r * cos(theta)
      y = r * sin(theta)
    }
    method(x,y,asp=1,...) 
  }

