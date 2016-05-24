ICLcalc <-
function(conv,n,zmat,bic,modnum,q,G){
	if(conv==0){
		ICL <- -Inf
	}
	if(conv==1){
		ent <- matrix(0, n, 1)
		zmap <- apply(zmat,1,which.max)
		for(i in 1:n){
			ent[i] <- zmat[i, zmap[i]]
		}
		ICL <- bic[modnum,q,G] + 2*sum(log(ent))
	}
	ICL
}
