finner.prod.iter <-
function(x) {
  f=.C("fast_inner_prod",as.double(x), as.integer(length(x)), as.double(rep(0,(length(x)-1))), as.double(rep(0,(length(x)-1))), as.double(rep(0,(length(x)-1))))[[3]]
  return(f)
}
