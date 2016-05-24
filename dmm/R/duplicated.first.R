duplicated.first <-
function(x){
# duplicated.first() - sets T for the first of each set of duplicate entries
  y <- rep(F,length(x))
  y[match(x[duplicated(x)],x)] <- T
  return(y)
}
