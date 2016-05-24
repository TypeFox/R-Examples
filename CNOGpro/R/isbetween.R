isbetween <-
function(x, lower, upper){
  return(x >= as.integer(lower) & x <= as.integer(upper))
}
