wupdate <-
function(x,n,G,mug,sigmainv,dfnewg,p,delta,w){
#	delta <- matrix(0,n,G)
	for(g in 1:G){
#		delta[,g] <- mahalanobis(x, mug[g,], sigmainv[,,g], inverted=TRUE)
		w[,g] <- (dfnewg[g]+p)/(dfnewg[g]+delta[,g])
	}
	w
}
