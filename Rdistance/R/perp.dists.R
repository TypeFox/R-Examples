perp.dists <- function(obs.dist, obs.angle, digits=1){
  
  # Stop and print error if dist or angle is NA
  if(any(is.na(obs.dist)))
    stop("Please remove detections for which obs.dist and/or obs.angle is NA.")
  if(any(is.na(obs.angle)))
    stop("Please remove detections for which obs.dist and/or obs.angle is NA.")
  
  # Convert degrees to radians, then apply sine function and round off decimal places
  round(obs.dist * sin(obs.angle*(pi/180)), digits)
    
}