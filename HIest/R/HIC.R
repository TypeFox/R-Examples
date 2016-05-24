HIC <-
function(G){
	n <- nrow(G)
	k <- ncol(G)
	if(is.null(names(G))){Locus <- 1:k}else{Locus=names(G)}
	P <- data.frame(Locus=Locus,Allele="1",P1=1,P2=0)
	G <- replace(G,G<0,NA)
	out <- matrix(nrow=n,ncol=3)
	out[,1] <- rowMeans(G,na.rm=TRUE)/2
	out[,2] <- rowMeans(G==1,na.rm=TRUE)
	G <- replace(G,is.na(G),-9)
	for(i in 1:n){
		out[i,3] <- HILL(out[i,1:2],G[i,],P,type="allele.count")
	}
	colnames(out) <- c("S","H","logLik")
	out
}
