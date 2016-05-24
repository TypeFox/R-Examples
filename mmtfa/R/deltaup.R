deltaup <-
function(x,mug,sigma,sigmainv,G,n,delta){
	for(g in 1:G){
#		delta[,g] <- mahalanobis(x, mug[g,],sigmainv[,,g], inverted=TRUE)
	  delta[,g] <- maha(x, mug[g,], sigma[,,g])
	}
	delta
}

