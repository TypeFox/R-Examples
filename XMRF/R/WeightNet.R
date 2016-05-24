WeightNet <-
function(net, weight, thw=0){
  wtemp = weight > thw
  nnet = net * wtemp
  return(nnet)
}
