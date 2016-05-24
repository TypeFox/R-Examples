
cindexes.W  <- function(lp, stime, status, groupe, ties, cindex, tau) {

	cindexg <- cpeg <- Npairsg <- comparableg <- unusableg <- concordanteg  <- discordanteg  <- tiedcompg <- tiedtotg  <- tiedtimeg <- uno.cindexg <- rep(0,length(unique(groupe)))   

	##### Calcul des poids selon Uno 
	if(tau==0) tau=max(stime[status==1])
	cens <- kmcens(stime, status, tau)
	GXi <- cens$surv[match(stime, cens$distinct, nomatch = 1)]
	Wipop <- 1/GXi/GXi * status * as.numeric(stime < tau)
	
	for(g in sort(unique(groupe)) ) {
		indiceg                 <- which(sort(unique(groupe))==g)
		cindexesg               <- cindexes(lp[groupe==g], stime[groupe==g], status[groupe==g], ties, tau, Wipop[groupe==g], cindex)
		Npairsg[indiceg]        <- cindexesg$Npairs
		cpeg[indiceg]           <- cindexesg$CPE
		tiedtotg[indiceg]       <- cindexesg$tiedtot
		uno.cindexg[indiceg]    <- cindexesg$c.uno
		
		if(cindex==1){	
			cindexg[indiceg]        <- cindexesg$cindex
			comparableg[indiceg]    <- cindexesg$comparable
			concordanteg[indiceg]   <- cindexesg$concordante
			discordanteg[indiceg]   <- cindexesg$discordante
			tiedcompg[indiceg]      <- cindexesg$tiedcomp
			tiedtimeg[indiceg]      <- cindexesg$tiedtime
			unusableg[indiceg]      <- cindexesg$unusable
		}

	}
	res.cpe <- mean(cpeg, na.rm=TRUE)
	uno.cindex <- mean(uno.cindexg, na.rm=TRUE)
	
	out <- list(Npairs=sum(Npairsg, na.rm=T),comparable=sum(comparableg, na.rm=T),tiedtot=sum(tiedtotg, na.rm=T),CPE=res.cpe, CPE.by.group=cpeg,c.uno=uno.cindex,c.uno.by.group=uno.cindexg)
	if(cindex==1){
		cindex_global <- mean(cindexg, na.rm=TRUE)
		out <- c(out,concordante=sum(concordanteg, na.rm=T),discordante=sum(discordanteg, na.rm=T),tiedcomp=sum(tiedcompg, na.rm=T),tiedtime=sum(tiedtimeg, na.rm=T),unusable=sum(unusableg, na.rm=T),cindex=cindex_global,cindex.by.group=cindexg) 
	}
	return(out)
}
              
             
