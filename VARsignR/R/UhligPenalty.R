UhligPenalty <-
function(g, first, last, constrained, impulses, scales, pen){
  func <- 0.0
  q <- matrix(stereo(v=g))
  for(k in first:last){
    ik <- (impulses[k, , ]%*%q) / scales
    for(i in 1:length(constrained)){
      if(constrained[i]<0){
        value <- ik[-1.0*constrained[i]]
      }else{
        value <- -1.0 * ik[constrained[i]]
      }
      if(value<0){
        func <- func +  value
      }else{
        func <- func + pen * value
      }
    }
  }
  acc <- func
  return(acc)
}
