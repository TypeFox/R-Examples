nxt <-
function(g,max.indeg) {
  j<-1
  s<-sum(g)
  cont<-TRUE
  while (cont) {
  	if (g[j]==0 && s<max.indeg) {
  		g[j]<-1
  		cont<-FALSE
  	} else {
  		if (g[j]==1) {s<-s-1}
  		g[j]<-0
  		j<-j+1
  	}
  	if (j>length(g)) {g<-rep(0,length(g));cont<-FALSE}
  }
  return(g)
}
