gpp <- function(NEE, PAR, ts.NEE, PAR.Temp, Reco.m, ts.Reco = NULL, method = "Michaelis-Menten", units = "30mins", allow.offset = FALSE, virtual = FALSE, start.par = max(PAR), ...){	
	# check method
	METHODS <- c("Michaelis-Menten", "Falge", "Smith", "Misterlich")
	method <- pmatch(method, METHODS)
	
	# make artificial NEE and PAR data for falling back to a "mean" model
	if(virtual){
		NEE <- c(seq(max(NEE, na.rm=TRUE), mean(NEE, na.rm=TRUE), length.out=5), rep(mean(NEE, na.rm=TRUE), 25))
		PAR <- c(seq(0, 200, length.out=10), seq(250, 1550, length.out=20))
		ts.NEE <- seq(min(ts.NEE, na.rm=TRUE), max(ts.NEE, na.rm=TRUE), length.out=30)
	}
	# check whether ts.Reco is provided. If so, it is assumed that Reco.m contains  
	# modelled reco values and PAR.Temp contains the relevant temperatures and at least contains the 
	# data for dates/times where the values in NEE and PAR were obtained
	if(!is.null(ts.Reco)){
		ts.NEEr <- round.POSIXlt(ts.NEE, digits=units)
		dat.NEE <- data.frame(NEE = NEE, PAR = PAR, dtr = ts.NEEr, dtrchar = as.character(ts.NEEr))
		dat.Reco <- data.frame(Reco = Reco.m, dt = ts.Reco, dtchar = as.character(ts.Reco))
		dat.all <- merge(dat.Reco, dat.NEE, by.x="dtchar", by.y="dtrchar")
		dat.all <- dat.all[order(dat.all$PAR),]
		NEE <- dat.all$NEE
		Reco <- dat.all$Reco
		PAR <- dat.all$PAR
		PAR.Temp <- NULL
		Reco.m <- NULL
	}
	else{
		# predict Reco values using the Reco model(s) provided
		# first extract and assign when there is a budget.reco results structure (class = "breco") provided
		if(class(Reco.m)=="breco"){
			t.diff <- abs(sapply(Reco.m, function(x) julian(x$ts)) - mean(julian(ts.NEE)))
			Reco.m <- Reco.m[which.min(t.diff)]
			which.Temp <- Reco.m[[1]]$which.Temp
			Reco.m <- Reco.m[[1]]$mod[[1]]
			PAR.Temp <- PAR.Temp[,which.Temp]
		}
		# then predict
		Reco <- predict(Reco.m, newdata = data.frame(Temp = PAR.Temp))
	}
	# calculate GPP
	GPP <- NEE - Reco
	# correct offset when wanted (default)
	offset <- 0
	if(!allow.offset){
		#offset <- as.numeric(coefs[1])
		offset <- max(GPP)
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
	# prepare output
	dat <- data.frame(NEE = NEE, GPP = GPP, Reco = Reco, PAR = PAR, timestamp = ts.NEE, mins = 0)
	if(!is.null(ts.Reco)){
			res <- list(mg = mg, data = list(dat = dat, offset = offset, start=s.list))
			class(res) <- "gpp2"
	}
	else{
		res <- list(mg = mg, mr = Reco.m, data = list(dat = dat, offset = offset, start=s.list, PAR.Temp = PAR.Temp))
		class(res) <- "gpp"
	}
	return(res)
}