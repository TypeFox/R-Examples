trans <-
function(s,w){
		#	added following line to avoid package notes
		m=get("m");
  		sddev=apply(s,2,sd)
		for(k in 1:m){s[,k]=(s[,k]-mean(s[,k]))/(sddev[k])}
		sthird=apply(s^3,2,mean)
		ssign=sign(sthird)
		for(k in 1:m){s[,k]=s[,k]*ssign[k]}
		sthird=sthird*ssign
		or=order(sthird)
		s=t(t(s)[or,])
		for(k in 1:m){w[,k]=w[,k]/(sddev[k])*ssign[k]}
		w=w[,or]
		return(c(list(s),list(w)))	
  }
