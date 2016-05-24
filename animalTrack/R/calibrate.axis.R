calibrate.axis <-
function(x,max,min,scale=1){
  total <- abs(min- max)
  dist <- total/2
  offset <- (max+min)/2
  scaled <- ((x - offset)/dist)*scale  
  return(scaled)
}

