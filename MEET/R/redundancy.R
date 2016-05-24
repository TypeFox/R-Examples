redundancy <-function(HX,HXmax,Herror,finite,correction=TRUE){
	
	R<-1-(HX/HXmax)
 	if (correction) { R<-correction.redundancy(r=R,HXmax,Herror,finite)}
R
}

