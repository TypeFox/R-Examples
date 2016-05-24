pitch2 <-
function(x,y,z,roll){
  pitch = atan((-x)/((y*sin(roll))+(z*cos(roll)))  )
  return(pitch)  
}

