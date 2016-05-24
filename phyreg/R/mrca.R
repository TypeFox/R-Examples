mrca <-
function(x,y,phy) {while(x!=y){if (x<y) {x<-phy[x]} else {y<-phy[y]}}; return(x)}
