
test.clusters<-function(x, ids.cluster1, ids.cluster2)    
{
	#class(x) <- "hclust"
	
	d.sub1 = rep(0, factorial(x$m))
	d.sub2 = rep(0, factorial(x$m))
	
	for (i in ids.cluster1)
	{
		d.sub1 <- d.sub1 + codebook(x$data[,i], m=x$m, normalized=F)
	}
	for (i in ids.cluster2)
	{
		d.sub2 <- d.sub2 + codebook(x$data[,i], m=x$m, normalized=F)
	}
	
	d.parent = d.sub1+d.sub2
	
	N.sub1 <- sum(d.sub1)
	N.sub2 <- sum(d.sub2)
	N.parent <- sum(d.parent)
	
	llr <- 2*( log(d.sub1/d.parent*N.parent/N.sub1)*d.sub1+log(d.sub2/d.parent*N.parent/N.sub2)*d.sub2 )
	
	lambda <- sum(llr, na.rm=T)
	
	p.value <-  1- pchisq(lambda, df=factorial(x$m)-1)
	#
	r<- c()
	r$id1 <- ids.cluster1
	r$id2 <- ids.cluster2
	r$lambda <- lambda
	r$p.value <- p.value
	r$sub1 <- d.sub1
	r$sub2 <- d.sub2
	
	return(r)
	
}


plot.add.pvalues <- function(x) {
clusts <- c()
ctrs <- c()
for (i in 1:dim(x$merge)[1]) {
	cur.merge<-x$merge[i,]
	
	if (cur.merge[1] < 0) { id1 <- -cur.merge[1]; ctr1<-which(x$order==id1)}	else {id1 <- clusts[[cur.merge[1]]]; ctr1<-ctrs[[cur.merge[1]]]}
	if (cur.merge[2] < 0) { id2 <- -cur.merge[2]; ctr2<-which(x$order==id2) } else {id2 <- clusts[[cur.merge[2]]]	
	ctr2 <- ctrs[[cur.merge[2]]]}
	

	result <- test.clusters(x, id1, id2)

	cat(result$p.value,"\n")

	clusts[[i]] <- c(id1,id2)
	ctrs[[i]] <- (ctr1+ctr2)/2
	#if (result$p.value <= 0.991)
	{ ptext( ctrs[[i]], x$height[i], paste("p=",signif( result$p.value,3),sep=""))}
}

}

ptext <- function(x,y,labels,...) 
{
	xoffset <- strwidth('p=XX')
	yoffset <- strheight('X')
	
	text(x+xoffset,y+yoffset,labels)
	
}