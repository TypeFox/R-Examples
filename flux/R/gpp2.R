gpp2 <- function(NEE, PAR, ts.NEE, oot, oot.id = c("D", "T"), method = "Michaelis-Menten", allow.offset = FALSE, virtual = FALSE, start.par = max(PAR), ...){
	
	# check method
	METHODS <- c("Michaelis-Menten", "Falge", "Smith", "Misterlich")
	method <- pmatch(method, METHODS)
	
	# make artificial NEE and PAR data for falling back to a "mean" model
	if(virtual){
		NEE <- c(seq(max(NEE, na.rm=TRUE), mean(NEE, na.rm=TRUE), length.out=5), rep(mean(NEE, na.rm=TRUE), 25))
		PAR <- c(seq(0, 200, length.out=10), seq(250, 1550, length.out=20))
		ts.NEE <- seq(min(ts.NEE, na.rm=TRUE), max(ts.NEE, na.rm=TRUE), length.out=30)
	}

	# compile data
	sel <- apply(as.matrix(dist(as.numeric(julian(ts.NEE)))*24*60)[oot==oot.id[1],oot==oot.id[2]], 2, which.min)
	mins <- apply(as.matrix(dist(as.numeric(julian(ts.NEE)))*24*60)[oot==oot.id[1],oot==oot.id[2]], 2, min)
	GPP <- NEE[as.numeric(names(sel))] - NEE[sel]
	PAR <- PAR[as.numeric(names(sel))]
	dat <- data.frame(NEE = NEE[as.numeric(names(sel))], GPP = GPP, Reco = NEE[sel], PAR = PAR, timestamp = ts.NEE[as.numeric(names(sel))], mins = mins)
	# correct offset when wanted (default)
	offset <- 0
	if(!allow.offset){
		#offset <- as.numeric(coefs[1])
		mGPP <- max(GPP)
		if(mGPP > 0){offset <- mGPP}
		GPP <- GPP - offset
	}
	# derive start values for the non-linear fitting
	# the start for alpha is taken from the linear regression
	# of GPP against PAR for PAR <= start.par
	sel <- PAR <= start.par
	coefs <- coef(lm(GPP ~ PAR, subset=sel))
	s.alpha <- coefs[2]
	if(is.na(s.alpha) | (s.alpha >= 0)){
		s.alpha <- -0.005
	}
	# the start value for GPmax is obtained by averaging the fluxes
	# at the 5 highest PAR values (if there are less than 5 PAR values take all)
	nopv <- 5
	if(length(PAR)<5){nopv <- length(PAR)}
	PAR.ord <- order(PAR)
	PAR.sel <- PAR.ord[(length(PAR.ord)-nopv) : length(PAR.ord)]
	s.GPmax <- mean(GPP[PAR.sel])
	# s.GPmax <- mean(GPP)
	# compile the start value list
	s.list <- list(GPmax = s.GPmax, alpha = s.alpha)
	# do the modeling
	tryGPP <- function(s.alpha, ...){
		if(method==1){
			try(nls(GPP ~ (GPmax * alpha * PAR)/(alpha * PAR + GPmax), start=list(GPmax = s.GPmax, alpha = s.alpha), model=TRUE, ...), silent=TRUE)
			}
		else if(method==2){
			try(nls(GPP ~ (alpha * PAR)/(1 - (PAR/2000) + (alpha*PAR/GPmax)), start=list(GPmax = s.GPmax, alpha = s.alpha), model=TRUE, ...), silent=TRUE)
			}		
		else if(method==3){
			try(nls(GPP ~ (alpha * PAR * GPmax)/sqrt(GPmax^2 + (alpha*PAR)^2), start=list(GPmax = s.GPmax, alpha = s.alpha), model=TRUE, ...), silent=TRUE)
			}		
		else if(method==4){
			try(nls(GPP ~ GPmax * (1 - exp((alpha*PAR) / GPmax)), start=list(GPmax = s.GPmax, alpha = s.alpha), model=TRUE, ...), silent=TRUE)
			}				
	}
	mgs <- lapply(seq(-0.001,-1,length.out=20), function(x) tryGPP(x))
	if(sum(sapply(mgs, class) != "try-error")==0){
		mg <- NA
		class(mg) <- "try-error"
	}
	else{
		mgs <- mgs[sapply(mgs, class) != "try-error"]
		mg <- mgs[which.min(sapply(mgs, function(x) x[]$convInfo$finIter))][[1]]
	}
	res <- list(mg = mg, data = list(dat = dat, offset = offset, start=s.list))
	class(res) <- "gpp2"
	return(res)
}