crossmemb<-function(x,y,relativize=TRUE) {
	if(class(x)=="vegclust" || class(x)=="vegclass") x = x$memb
	if(class(y)=="vegclust" || class(y)=="vegclass") y = y$memb
	c=t(x)%*%as.matrix(y)
	if(relativize) c = sweep(c,1,colSums(x),"/")
    return(c)
}