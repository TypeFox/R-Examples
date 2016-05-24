detail=function(x){
res <- list(x=x,
	    isS4=isS4(x),
	    is.object=is.object(x),
	    class=class(x),
	    attributes=attributes(x))
res
}