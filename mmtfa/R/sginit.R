sginit <-
function(p,G,x,mug,zmat,n,ng){
	sg <- array(0, dim=c(p, p, G))
#	for(g in 1:G){
#		for(i in 1:n){
#			xminus <- x[i,]-mug[g,]
#			sg[,,g] <- sg[,,g] + (zmat[i,g]/ng[g])*(xminus%*%t(xminus))
#		}
#	}
	for(g in 1:G){
		sg[,,g] <- cov.wt(x,wt=zmat[,g],method="ML")$cov
	}
	sg 
}
