angle <-
function(cut.slope, ref.slope, tick.dist){
  cut.slope <- cut.slope/tick.dist
  ref.slope <- ref.slope/tick.dist
  acos((1 + cut.slope * ref.slope)/(sqrt((1 + cut.slope^2)*(1 + ref.slope^2)))) #/
  #tick.dist #/ 4 * pi
}
