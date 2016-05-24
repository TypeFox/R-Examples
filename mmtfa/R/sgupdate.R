sgupdate <-
function(p,G,n,x,mug,zmat,w,ng,mod,pig,sg){
	#for(g in 1:G){
#			for(i in 1:n){
#				xminus <- x[i,]-mug[g,]
#				sg[,,g] <- sg[,,g] + (zmat[i,g]*w[i,g]/ng[g])*
#						(xminus%*%t(xminus))
#			} 
#	}
	for(g in 1:G){
		#sg[,,g] <- cov.wt(x,wt=zmat[,g]*w[,g],method="ML",center=mug[g,])$cov
	  sg[,,g] <- cov.wt(x,wt=zmat[,g]*w[,g],method="ML")$cov
	  #sg[,,g] <- cov.wt(x,wt=c(zmat[,g])*c(w[,g]),method="ML")$cov
	}
	if(substring(mod,1,2)=="CC"){
	  sgc <- matrix(0,p,p)
		#for(g in 1:G){
#			sgc <- sgc + pig[g]*sg[,,g]
#		}
		for(g in 1:G){
			#sgc <- sgc + pig[g]*cov.wt(x,wt=zmat[,g]*w[,g],method="ML")$cov
			sgc <- sgc + pig[g]*sg[,,g]
		}
		for(g in 1:G){
			sg[,,g] <- sgc
		}
	}
	sg	
}
