computePSS <-
function(combs, ra){
  sum(apply(ra[combs,],2,sum)^2)
}
