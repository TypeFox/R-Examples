delta <-
function(SpA, SpB){
  sum(apply(cbind(SpA, SpB), 1, function(x){abs(diff(x))/2}))
}

