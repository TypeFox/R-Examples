
.make_lookup = function(cx) {
  
  f = function(j, x, cx) {
    if(x==1) {
      -1
    }
    else if(j==1) {
      x
    } else {
      sum(cx[1:(j-1)]-1)+x
    }
  }
  
  ta = do.call(rbind, lapply(1:length(cx), function(x) data.frame(j=x, c=1:cx[x])))
  cbind(ta, ix = apply(ta, 1, function(x) f(x[1], x[2], cx)))
  
}


.make_lookup_centered = function(cx) {
  
  f = function(j, x, cx) {
    if(j==1) {
      x
    } else {
      sum(cx[1:(j-1)])+x
    }
  }
  
  ta = do.call(rbind, lapply(1:length(cx), function(x) data.frame(j=x, c=1:cx[x])))
  cbind(ta, ix = apply(ta, 1, function(x) f(x[1], x[2], cx)))
  
}
