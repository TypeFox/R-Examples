getAIC <-
function(L,np,n,AICc=TRUE){
	result<-(-2)*(L)+2*np
	if(AICc)result<-result+2*np*((np+1)/(n-np-1))
	result
	}
