dfupdatefun <- function(mod,dfnewg,ng,zmat,w,G,p,n,x,mug,sigmainv){
	if(substring(mod,4,4)=="U"){
		dfoldg <- dfnewg
		for(g in 1:G){
      constn <-  1+ (1/ng[g]) * sum(zmat[,g]*(log(w[,g])-w[,g])) + digamma((dfoldg[g]+p)/2) - log((dfoldg[g]+p)/2)
			dfnewg[g] <- uniroot( function(v) log(v/2) - digamma(v/2) + constn, 
					lower=.0001, upper=1000 )$root
			if(dfnewg[g]>200){
				dfnewg[g]<-200
			}
			if(dfnewg[g]<2){
				dfnewg[g]<-2
			}
		}
	}
	if(substring(mod,4,4)=="C"){
		dfoldg <- dfnewg[1]
		constn <- 1+ (1/n) * sum(zmat * (log(w)-w)) + digamma((dfoldg+p)/2) - log((dfoldg+p)/2)
		dfsamenewg <- uniroot( function(v) log(v/2) - digamma(v/2) + constn , lower=.0001, upper=1000 )$root
		if(dfsamenewg>200){
			dfsamenewg<-200
		}
		if(dfsamenewg<2){
			dfsamenewg<-2
		}
		dfnewg <- c(rep(dfsamenewg,G))
	}
	dfnewg
}
