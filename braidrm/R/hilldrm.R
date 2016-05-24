findBestHill <- function(model,data,defaults,startparv=NULL,llims=NULL,ulims=NULL,...) UseMethod("findBestHill")

findBestHill.default <- function(model,data,defaults,startparv=NULL,llims=NULL,ulims=NULL,...) {
	conc <- model
	act <- data
	if (is.null(startparv)) { startparv <- c(NA,NA,NA,NA) }
	drdef <- dr_start(conc,act)
	def <- startparv
	def[which(is.na(def))] <- drdef[which(is.na(def))]
	# medAct <- min(max(median(act),min(def[c(2,3)])+0.001),max(def[c(2,3)])-0.001)
	medAct <- (def[2]+def[3])/2
	if (defaults[1]<defaults[2]) {
		if (drdef[2]>drdef[3]) { drdef[c(2,3)] <- drdef[c(3,2)] }
		if (def[2]>def[3]) { def[c(2,3)] <- def[c(3,2)] }
		mlims <- array(c(0,2*min(drdef[2],def[2])-max(drdef[3],def[3]),medAct+0.001,1.5*min(log(conc))-0.5*max(log(conc)),
						0,medAct-0.001,2*max(drdef[3],def[3])-min(drdef[2],def[2]),1.5*max(log(conc))-0.5*min(log(conc))),dim=c(4,2))
	} else {
		if (drdef[2]<drdef[3]) { drdef[c(2,3)] <- drdef[c(3,2)] }
		if (def[2]<def[3]) { def[c(2,3)] <- def[c(3,2)] }
		mlims <- array(c(0,medAct+0.001,2*max(drdef[3],def[3])-min(drdef[2],def[2]),1.5*min(log(conc))-0.5*max(log(conc)),
						0,2*min(drdef[2],def[2])-max(drdef[3],def[3]),medAct-0.001,1.5*max(log(conc))-0.5*min(log(conc))),dim=c(4,2))
	}
	def[4] <- max(min(def[4],mlims[4,2]),mlims[4,1])
	if (def[1]<0) { mlims[1,] <- c(-10,-0.01) }
	else { mlims[1,] <- c(0.01,10) }
	if (!is.null(llims)) { mlims[!is.na(llims),1] <- pmax(mlims[,1],llims)[!is.na(llims)] }
	if (!is.null(ulims)) { mlims[!is.na(ulims),2] <- pmax(mlims[,2],ulims)[!is.na(ulims)] }
	
	pindl <- list(c(1,4),c(1,2,4),c(1,3,4),c(1,2,3,4),c(3))
	allfits <- list()
	for (i in 1:4) {
		fixed <- c(0,defaults[1],defaults[2],0)
		fixed[pindl[[i]]] <- NA
		opdef <- drdef
		opdef[which(!is.na(fixed))] <- def[which(!is.na(fixed))]
		mod <- fitHillDrm(conc,act,fixed=fixed,llims=mlims[,1],ulims=mlims[,2],def=opdef[which(is.na(fixed))])
		coefs <- fixed
		coefs[pindl[[i]]] <- mod$par
		allfits[[i]] <- list(coefficients=coefs,parv=mod$par,pinds=pindl[[i]])
	}
	allfits[[5]] <- list(coefficients=c(mean(act)),parv=c(mean(act)),pinds=c(3))
	names(allfits) <- c("m2P","m3Puc","m3Plc","m4P","Lin")
	
	modAIC <- rep(0,times=5)
	for (i in 1:5) {
		modAIC[i] <- getHillAIC(allfits[[i]],conc,act)
		allfits[[i]]$AIC <- modAIC[i]
	}
	bmi <- modelSelect(modAIC, c(2,3,3,4,1))
	fit_l <- list(conc=conc,act=act,bestModIdx=bmi,bestModName=names(allfits)[bmi],allfits=allfits,mlims=mlims)
	return(fit_l)
}
findBestHill.formula <- function(model,data,...) {
	mf <- model.frame(formula=model, data=data)
	conc <- model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(conc,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { conc <- conc[,-i] }
	}
	conc <- as.vector(conc)
	act <- model.response(mf)
	return(findBestHill.default(conc,act,...))
}
getHillBootstrap <- function(hfitl,ciLevs=c(0.025,0.975),mi=NULL) {
	if (is.null(mi)) { mi <- hfitl$bestModIdx }
	nhfitl <- hfitl
	if (!is.null(nhfitl$ciPass)) {
		nhfitl$ciPass <- NULL
		nhfitl$ciLevs <- NULL
		nhfitl$ciMInd <- NULL
		nhfitl$bCoefs <- NULL
		nhfitl$ciVec <- NULL
	}
	
	prcBootPass <- 0.5
	numBoot <- max(min(10/(1-ciLevs[2]+ciLevs[1]),1000),100)
	def <- coef(hfitl$allfits[[mi]])
	pred <- evalHillEqn(hfitl$conc,def)
	res <- hfitl$act-pred
	
	pinds <- hfitl$allfits[[mi]]$pinds
	fixed <- def
	fixed[pinds] <- NA
	bCoefs <- array(0,dim=c(0,4))
	for (i in 1:numBoot) {
		cact <- pred + sample(res,length(res),replace=TRUE)
		mod <- fitHillDrm(hfitl$conc,cact,fixed=fixed,def=def,llims=hfitl$mlims[,1],ulims=hfitl$mlims[,2])
		coefs <- def
		coefs[pinds] <- mod$par
		bCoefs <- rbind(bCoefs,coefs)
	}
	if (nrow(bCoefs)<numBoot*prcBootPass) {
		nhfitl$ciPass <- FALSE
	} else {
		qMat <- apply(bCoefs,2,quantile,probs=ciLevs)
		nhfitl <- c(nhfitl,list(ciPass=TRUE,ciLevs=ciLevs,ciMInd=mi,bCoefs=bCoefs,ciVec=as.vector(qMat)))
	}
	return(nhfitl)
}

hillConcCorrect <- function(conc,act,parv,sigr=1) {
	concvec <- exp(seq(min(log(conc[conc>0])),max(log(conc[conc>0])),length=1000))
	actvec <- evalHillEqn(concvec,parv)
	cconc <- conc
	for (i in 1:length(conc)) {
		if (conc[i]==0) { next }
		cconc[i] <- concvec[which.min((log10(concvec)-log10(conc[i]))^2+(1/sigr)^2*(actvec-act[i])^2)]
	}
	return(cconc)
}

evalHillEqn <- function(conc,parv) { return(evalHillEqn_int(conc,parv)) }
evalHillEqn_int <- function(conc,parv,fixed=c(NA,NA,NA,NA),calcderivs=FALSE) {
	Det <- 1+(conc/exp(parv[4]))^(-parv[1])
	Ei <- parv[2]+(parv[3]-parv[2])/Det
	if (!calcderivs) { return(Ei) }
	derivs <- array(0,dim=c(length(conc),4))
	if (is.na(fixed[1])) { derivs[,1] <- ((parv[3]-parv[2])/(Det^2))*((conc/exp(parv[4]))^(-parv[1]))*log(conc/exp(parv[4])) }
	if (is.na(fixed[2])) { derivs[,2] <- 1-1/Det }
	if (is.na(fixed[3])) { derivs[,3] <- 1/Det }
	if (is.na(fixed[4])) { derivs[,4] <- -parv[1]*((parv[3]-parv[2])/(Det^2))*((conc/exp(parv[4]))^(-parv[1])) }
	derivs[which(is.nan(derivs))] <- 0
	return(cbind(Ei,derivs))
}
invertHillEqn <- function(val,parv) {
	concout <- rep(0,times=length(val))
	if (sign(parv[3]-parv[2])*sign(parv[1])>0) {
		concout[val>=max(parv[2],parv[3])] <- Inf
		concout[val<=min(parv[2],parv[3])] <- 0
	} else {
		concout[val>=max(parv[2],parv[3])] <- 0
		concout[val<=min(parv[2],parv[3])] <- Inf
	}
	Eval <- (val-parv[2])/(parv[3]-val)
	cinds <- which(Eval>0)
	concout[cinds] <- exp(parv[4])*Eval[cinds]^(1/parv[1])
	return(concout)
}

fitHillDrm <- function(conc,act,fixed=c(NA,NA,NA,NA),llims=NULL,ulims=NULL,def=NULL) {
	npar <- length(which(is.na(fixed)))
	if (is.null(def)) {
		full_def <- dr_start(conc,act)
		def <- full_def[which(is.na(fixed))]
		full_def[which(!is.na(fixed))] <- fixed[which(!is.na(fixed))]
	} else if (length(def)==4) {
		full_def <- def
		def <- full_def[which(is.na(fixed))]
	} else if (length(def)==npar) {
		full_def <- fixed
		full_def[which(is.na(fixed))] <- def
	} else { stop("Number of starting values must match number of free parameters.") }
	
	fpars <- which(is.na(fixed))
	if (full_def[2]<full_def[3]) {
		full_llims <- c(0,2*full_def[2]-full_def[3],(full_def[2]+full_def[3])/2,full_def[4]-3)
		full_ulims <- c(0,(full_def[2]+full_def[3])/2,2*full_def[3]-full_def[2],full_def[4]+3)
	} else {
		full_llims <- c(0,(full_def[2]+full_def[3])/2,2*full_def[3]-full_def[2],full_def[4]-3)
		full_ulims <- c(0,2*full_def[2]-full_def[3],(full_def[2]+full_def[3])/2,full_def[4]+3)
	}
	if (full_def[1]<0) {
		full_llims[1] <- -10
		full_ulims[1] <- -0.1
	} else {
		full_llims[1] <- 0.1
		full_ulims[1] <- 10
	}
	if (!is.null(llims)) {
		if (length(llims)==4) { full_llims[!is.na(llims)] <- pmax(full_llims,llims)[!is.na(llims)] }
		else if (length(llims)==npar) { full_llims[fpars[!is.na(llims)]] <- pmax(full_llims,llims)[fpars[!is.na(llims)]] }
		else { stop("Number of lower limit values must match number of free parameters.") }
	}
	if (!is.null(ulims)) {
		if (length(ulims)==4) { full_ulims[!is.na(ulims)] <- pmin(full_ulims,ulims)[!is.na(ulims)] }
		else if (length(ulims)==npar) { full_ulims[fpars[!is.na(ulims)]] <- pmin(full_ulims,ulims)[fpars[!is.na(ulims)]] }
		else { stop("Number of upper limit values must match number of free parameters.") }
	}
	wllims <- full_llims[fpars]
	wulims <- full_ulims[fpars]
	if (length(def)==4) { def <- def[which(is.na(fixed))] }
	fixed[which(!is.na(fixed))] <- full_def[which(!is.na(fixed))]
	
	parvHillErr <- function(parv) {
		fullparv <- fixed
		fullparv[which(is.na(fixed))] <- parv
		ei <- evalHillEqn_int(conc,fullparv,fixed,calcderivs=FALSE)
		err <- act-ei
		return(sum(err^2))
	}
	parvHillDerivs <- function(parv) {
		fullparv <- fixed
		fullparv[which(is.na(fixed))] <- parv
		ederivs <- evalHillEqn_int(conc,fullparv,fixed,calcderivs=TRUE)
		oderivs <- rep(0,times=4)
		err <- act-ederivs[,1]
		for (i in 1:4) { oderivs[i] <- -2*sum(err*ederivs[,i+1]) }
		oderivs <- oderivs[which(is.na(fixed))]
		return(oderivs)
	}
	
	parscale <- c(1,0.01*abs(full_def[3]-full_def[2]),0.01*abs(full_def[3]-full_def[2]),1)
	control <- list(parscale=parscale[which(is.na(fixed))])
	nls <- optim(def,parvHillErr,parvHillDerivs,method="L-BFGS-B",lower=wllims,upper=wulims,control=control,hessian=FALSE)
	nls$mlims <- rbind(llims,ulims)
	return(nls)
}
dr_start <- function(conc,act) { 
    tact <- act[which(conc>0)]
	tconc <- conc[which(conc>0)]
	
	scaleInc <- 0.001
	yRange <- quantile(tact,c(0.05,0.95),names=FALSE)
	tact <- pmin(pmax(tact,yRange[1]),yRange[2])
	lenyRange <- scaleInc * diff(yRange)
	yMin <-  yRange[1] - lenyRange
	yMax <-  yRange[2] + lenyRange
    
	lhs <- log((tact - yMin)/(yMax - tact))
	rhs <- log(tconc)
	lmOut <- lm(lhs ~ rhs, data = data.frame(lhs, rhs)[which(!is.infinite(lhs)), ])
	hill <- coef(lmOut)[[2]]
	logEC50 <- coef(lmOut)[[1]] / -hill
	if (hill<0) { out <- c(-hill,yMax,yMin,logEC50) }
	else { out <- c(hill, yMin, yMax, logEC50) }
	return(out)
}
getHillAIC <- function(mod,conc,act) {
	coefs <- coef(mod)
	if (length(coefs)==4) { pred <- evalHillEqn(conc,coefs) }
	else if (length(coefs)==1) { pred <- rep(coefs[1],times=length(conc)) }
	else { stop("Model is not valid!") }
	df <- length(conc)
	ssErr <- sum((pred-act)^2)
	lLik <- -(df/2)*(log(2*pi)+log(ssErr)-log(df)+1)
	k <- length(mod$parv)
	aic <- 2*k - 2*lLik + 2*k*(k+1)/(length(conc)-k-1)
	return(aic)
}
modelSelect <- function(modAIC, modK, ik=2) {
	evRatioCut <- 15
	bmi <- which.min(ifelse(modK==ik,modAIC,NA))
	while(any(exp(0.5*(modAIC[bmi]-modAIC[which(modK>modK[bmi])]))>evRatioCut)) {
		bmi <- which.min(ifelse(modK==modK[bmi]+1,modAIC,NA))
	}
	return(bmi)
}
