ANDNet <-
function(net, th=0){
  net = abs(net)
  tmp = net > th
  tmp = tmp * t(tmp)
  return(tmp)
}
