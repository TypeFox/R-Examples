stability <-
function(tree,B=100,graph=TRUE)
{
	if (!inherits(tree, "hclustvar")) 
        	stop("use only with \"hclustvar\" objects")
	cl <- match.call()
	X.quali <- tree$X.quali
	X.quanti <- tree$X.quanti
	clusmat <- tree$clusmat
	nmax <- ncol(clusmat)-1
	matRandCorrige <- matrix(,nrow=B,ncol=nmax-1)
	matRandNonCorrige <- matrix(,nrow=B,ncol=nmax-1)
	for (b in 1:B)
	{
		#print(b)
		res1<-bootvar(X.quanti,X.quali)
		Xboot.quanti<-res1$Xboot.quanti
		Xboot.quali<-res1$Xboot.quali
      	# a ameliorer peut-etre
		if (!is.null(Xboot.quanti)&& !is.null(Xboot.quali)) 
			{clusmatboot <- hclustvar(Xboot.quanti,Xboot.quali)$clusmat}
		if (is.null(Xboot.quanti)&& !is.null(Xboot.quali)) 
			{clusmatboot <- hclustvar(X.quali=Xboot.quali)$clusmat}
		if (!is.null(Xboot.quanti)&& is.null(Xboot.quali)) 
			{clusmatboot <- hclustvar(X.quanti=Xboot.quanti)$clusmat}	
		for (i in 2:nmax) 
			{matRandCorrige[b,i-1] <- rand(clusmat[,i],clusmatboot[,i],adj=TRUE)}
		for (i in 2:nmax) 
			{matRandNonCorrige[b,i-1] <- rand(clusmat[,i],clusmatboot[,i],adj=FALSE)}
	}
	colnames <- paste("P", 2:nmax,sep = "")
	rownames <- paste("B", 1:B, sep = "")
	colnames(matRandCorrige) <- colnames
	rownames(matRandCorrige) <- rownames
	critRandmoyenCorrige<-apply(matRandCorrige,2,mean)
	if (graph!=FALSE) {
	plot(x=seq(1,nmax-1),critRandmoyenCorrige,xaxt="n",ylim=c(0,1),xlab="number of clusters",ylab="mean adjusted Rand criterion",main="Stability of the partitions",type="b")
	axis(side=1,at=seq(1,nmax-1),labels=paste(2:nmax))
	}
	#dev.new()
	#colnames(matRandNonCorrige) <- colnames
	#rownames(matRandNonCorrige) <- rownames
	#critRandmoyenNonCorrige<-apply(matRandNonCorrige,2,mean)
	#plot(x=seq(1,nmax-1),critRandmoyenNonCorrige, xaxt="n",ylim=c(0,1),xlab="number of clusters",ylab="mean Rand criterion",main="Stability of the partitions",type="b")
	#axis(side=1,at=seq(1,nmax-1),labels=paste(2:nmax))

	retlist <- list(call=cl,matCR=matRandCorrige,meanCR=critRandmoyenCorrige)
	class(retlist) <- "clustab"
   	retlist

}

