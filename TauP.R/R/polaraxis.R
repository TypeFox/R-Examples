polaraxis <-
function(rp = 6371, at = 0:17 * 20){
  text( sin(at * pi/180) * rp * 1.1, cos(at * pi/180) * rp * 1.1, at)
}

