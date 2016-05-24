

cindexes  <- function(lp, stime, status, ties, tau, Wi, cindex) { 
	n <- length(lp);
	sumCPE <- 0;
	comparable  <- unusable  <- concordante  <- discordante  <- tiedcomp <- tiedtot  <- tiedtime  <- sumCPEj <- uno.comparable <- uno.concordante <- uno.tiedcomp <- rep(NA,(n-1)) 
	Npairs <- n*(n-1)/2
  
	for (i in 1:(n-1)) {
		# C Harrell
		if(cindex==1){
			tiedtimej     <- (stime[i]==stime & status[i]==1 & status==1)
			tiedtime[i]   <- sum(tiedtimej[(i+1):n])
			unusablej     <- tiedtimej | ((stime[i]<stime & status[i]==0) | (stime<stime[i] & status==0)
						| (stime[i]==stime & status[i]==0 & status==0))
			unusable[i]   <- sum(unusablej[(i+1):n])
			comparablej   <- ! unusablej
			comparable[i] <- sum(comparablej[(i+1):n])
			concordantej  <- (comparablej & ((stime[i]<stime & lp[i]>lp) | (stime<stime[i] & lp>lp[i])
					| (stime[i]==stime & status[i]==0 & status==1 & lp>lp[i])
					| (stime[i]==stime & status[i]==1 & status==0 & lp<lp[i])))
			concordante[i]<- sum(concordantej[(i+1):n])
			discordantej  <- (comparablej & ((stime[i]>stime & lp[i]>lp) | (stime>stime[i] & lp>lp[i]) 
					| (stime[i]==stime & status[i]==0 & status==1 & lp<lp[i])
					| (stime[i]==stime & status[i]==1 & status==0 & lp>lp[i])))
			discordante[i]<- sum(discordantej[(i+1):n])
			tiedcompj         <- (comparablej & lp[i]==lp)
			tiedcomp[i]       <- sum(tiedcompj[(i+1):n])  
		}
		tiedtotj         <- (lp[i]==lp)
		tiedtot[i]       <- sum(tiedtotj[(i+1):n])
	
		#CPE Gonen & Heller
		bxjxi <- lp - lp[i]
		bxixj <- - bxjxi
		if (ties==1) sumCPEj[i] <- sum(((bxjxi<=0)/(1+exp(bxjxi)) + (bxixj<0)/(1+exp(bxixj)))[(i+1):n])
		if (ties==0) sumCPEj[i] <- sum(((bxjxi<0)/(1+exp(bxjxi)) + (bxixj<0)/(1+exp(bxixj)))[(i+1):n])
		
		#C-Uno
		uno.comparablej   <- status[i] * Wi[i] * (stime[i]<stime & stime[i] < tau)
		uno.comparable[i] <- sum(as.numeric(uno.comparablej))
		uno.concordantej  <- (status[i] * Wi[i] * (stime[i]<stime & stime[i] < tau) * (lp[i]>lp))
		uno.concordante[i]<- sum(as.numeric(uno.concordantej))
		uno.tiedcompj     <- (status[i] * Wi[i] * (stime[i]<stime & stime[i] < tau) * (lp[i]==lp))
		uno.tiedcomp[i]   <- sum(as.numeric(uno.tiedcompj))
	}
	# patient n pour Uno
	uno.comparablej   <- status[n] * Wi[n] * (stime[n]<stime & stime[n] < tau)
	uno.comparable[n] <- sum(as.numeric(uno.comparablej))
	uno.concordantej  <- (status[n] * Wi[n] * (stime[n]<stime & stime[n] < tau) * (lp[n]>lp))
	uno.concordante[n]<- sum(as.numeric(uno.concordantej))
	uno.tiedcompj     <- (status[n] * Wi[n] * (stime[n]<stime & stime[n] < tau) * (lp[n]==lp))
	uno.tiedcomp[n]   <- sum(as.numeric(uno.tiedcompj))
	
	if (ties==1) {
		if(cindex==1) cindex_global <- (sum(concordante + tiedcomp/2))/sum(comparable)
		res.cpe <- (2/(n*(n-1))) * sum(sumCPEj)
		uno.cindex <- (sum(as.numeric(uno.concordante) + as.numeric(uno.tiedcomp)/2)/sum(as.numeric(uno.comparable)))
	}
	if (ties==0) {
		if(cindex==1) cindex_global <- sum(concordante)/(sum(comparable)-sum(tiedcomp))
		res.cpe <- (1/((n*(n-1)/2)-sum(tiedtot))) * sum(sumCPEj)
		uno.cindex <- sum(as.numeric(uno.concordante))/sum(as.numeric(uno.comparable))
	}

	out <- list(CPE=res.cpe,Npairs=Npairs,comparable=sum(comparable),tiedtot=sum(tiedtot),c.uno=uno.cindex)
	if(cindex==1){
		out <- c(out,concordante=sum(concordante),discordante=sum(discordante),
		tiedcomp=sum(tiedcomp),tiedtime=sum(tiedtime),unusable=sum(unusable),cindex=cindex_global)
	
	}
	return(out)
	
}
