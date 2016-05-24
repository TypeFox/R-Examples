zupdate <-
function(x,G,pig,dfnewg,p,yg,q,betag,lg,mug,sigmainv,n,clas,kno,known,unkno,delta){
	#num <- matrix(0,n,G)
	log_num <- matrix(0,n,G)
	for(g in 1:G){
		log_num[,g]<-log(pig[g])+lgamma((dfnewg[g]+p)/2)-(1/2)*
			#log((prod(diag(yg[,,g]))/det(diag(q)-betag[,,g]%*%lg[,,g])))
			(sum(log(diag(yg[,,g])))-log(det(diag(q)-betag[,,g]%*%lg[,,g])))-
			((p/2)*(log(pi)+log(dfnewg[g]))+
			lgamma(dfnewg[g]/2)+((dfnewg[g]+p)/2)*(log(1+ (delta[,g]/dfnewg[g]))))
#		log_num[,g]<-(log(pig[g])-(1/2)*	(sum(log(diag(yg[,,g])))-log(det(diag(q)-betag[,,g]%*%lg[,,g])))-
#			(((dfnewg[g]+p)/2)*(log(1+ (mahalanobis(x, mug[g,], sigmainv[,,g], inverted=TRUE)/dfnewg[g])))))
	}

	# CATCH NEGATIVE DETERMINANTS
# 			if(any(is.nan(log_num))){
# 				message("NEGATIVE DETERMINANT (probably)");
# 				#break
# 			}
	#print(num)
	#log_num <- log_num - apply(log_num,1,max)
	kcon <- -rowMaxs(log_num)
  log_num <- log_num + kcon
	num <- exp(log_num)
	zmat <- num/rowSums(num)
	
	if(clas>0){
	### REPLACE KNOWNS
		zmat <- unkno*zmat
		for(i in 1:n){
			if(kno[i]==1){
				zmat[i, known[i]] <- 1
			}
		}
	}
	
	#message("HEYYYYY")
	store <- list()
	store[["zmat"]] <- zmat
	store[["num"]] <- num
	store[["kcon"]] <- kcon
	store
}
