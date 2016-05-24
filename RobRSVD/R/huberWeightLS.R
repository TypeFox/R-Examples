huberWeightLS <-
function(data, k){
  w = k / abs(data)
  w[w>1] = 1
  return(w)
}
