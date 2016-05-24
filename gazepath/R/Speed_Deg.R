Speed_Deg <-
function(X, Y, distance, height_mm, width_mm, height_px, width_px, Hz){
  hor <- atan((width_mm / 2) / distance) * (180 / pi) * 2 / width_px * X
  ver <- atan((height_mm / 2) / distance) * (180 / pi) * 2 / height_px * Y
  speed <- sqrt( diff(hor,2) ^ 2 + diff(ver,2) ^ 2) * (Hz/2)
  return(c(.001, speed, 0.001))
}
