mdist <-
function(x,y){  
	m=length(x)
	n=length(y)
	if(m!=n){
		cat("Error: dimensions do not match","\n")
		return (-1)
	}
	out=0
	for (i in 1:n){
	deltax=as.numeric(I(x[i]!=0))
	deltay=as.numeric(I(y[i]!=0))
	out=out+as.numeric(I(deltax!=deltay))+(x[i]-y[i])^2
	}
	out=sqrt(out)
	return(out)
}
