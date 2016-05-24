CalcCorners = function(PZ, f = 1:1000000/1000, PLOT = FALSE){
  R = abs(PZ2Resp(PZ, f, PLOT))
  bounds = which(diff(sign(R - PZ$Sense/sqrt(2))) != 0)
  if(length(bounds) == 0){
    return(NULL)
  }
  if(bounds[1] == 1){
    stop('Low corner not detected: min(f) is too high')
  }
  fc = (f[bounds] + f[bounds+1])/2
  if(PLOT){
    abline(v = fc)
  }
  return(fc)
}
