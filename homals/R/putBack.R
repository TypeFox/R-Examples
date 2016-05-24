`putBack` <-
function(n,blocklist,blockvalues){
  x<-rep(0,n);nb<-length(blockvalues)
  for (i in 1:nb) {
  		x[blocklist[i,1]:blocklist[i,2]]<-blockvalues[i]}
  return(x)
}

