pitch <-
function(x, y, z){  
  pitch90 = -atan(x/sqrt(y^2 + z^2))  
  return(pitch90)  
}

