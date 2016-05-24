insert <-
function(v,e,pos){
  return(c(v[1:(pos-1)],e,v[(pos):length(v)]))
}
