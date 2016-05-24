MortSmooth_tpower <-
function(x, t, p){
  ## Input:
  ## x = abcissae of data
  (x - t) ^ p * (x > t)
  ## (x-t)^p gives the curve
  ## (x>t) is an indicator function; it is 1 when x>t
  ## and 0 when x<=t, i.e. before each knot
}
