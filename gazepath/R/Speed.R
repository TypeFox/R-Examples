Speed <-
function(X, Y, distance, height_mm, width_mm, height_px, width_px, res_x = 1280, res_y = 1024, Hz){
  d2p_px <- sqrt(abs(X - (res_x / 2))**2 + abs(Y - (res_y / 2))**2 + (distance * width_px / width_mm)**2)
  dbp_px <- c(.001, sqrt(diff(X, 2)**2 + diff(Y, 2)**2), .001)
  speed <- atan((dbp_px / 2) / d2p_px) * (180 / pi) * 2 * (Hz / 2)
}
