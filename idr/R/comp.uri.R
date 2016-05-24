comp.uri <-
function(tv, x){
  n <- length(x)
  qt <- quantile(x, prob=1-tv[1]) # tv[1] is t
  sum(x[1:ceiling(n*tv[2])] >= qt)/n

}

