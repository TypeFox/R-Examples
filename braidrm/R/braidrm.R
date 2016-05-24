summary.braidrm <- function(object, ...) {
	if (is.null(object$ciPass) || !object$ciPass) {
		TAB <- coef(object)
	} else {
		cv <- object$ciVec
		TAB <- cbind(CILow = cv[seq(1,length(cv),by=2)],Estimate = coef(object),CIHigh = cv[seq(2,length(cv),by=2)])
	}
	res <- list(call=object$call,coefficients=TAB)
	class(res) <- "summary.braidrm"
	return(res)
}
print.summary.braidrm <- function(x, ...) {
	cat("Call:\n")
	print(x$call)
	cat("\n")
	print(round(x$coefficients,digits=4))
}
print.braidrm <- function(x, ...) {
	cat("Call:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(x$fullpar)
}

braidrm <- function(model,data,getCIs=TRUE,fixed="kappa2",startparv=NULL,llims=NULL,ulims=NULL,...) UseMethod("braidrm")

braidrm.default <- function(model,data,getCIs=TRUE,fixed="kappa2",startparv=NULL,llims=NULL,ulims=NULL,...) {
	concs <- model
	act <- data
	if (ncol(concs)!=2) { stop("Parameter 'concs' must be an array with two columns.") }
	conc1 <- as.vector(concs[,1])
	conc2 <- as.vector(concs[,2])
	act <- as.vector(act)
	if (length(conc1)!=length(conc2) | length(conc1)!=length(act)) {
		stop("Parameters 'concs', and 'act' must contain the same number of elements.")
	}
	cdef <- startparv
	cfixed <- fixed
	for (iter in 1:5) {
		nls <- try(fitBRAIDrsm(conc1,conc2,act,fixed=cfixed,def=cdef,llims=llims,ulims=ulims),silent=TRUE)
		if (class(nls) != "try-error" && nls$convergence==0) { break }
		else if (iter==5) { stop(paste("Convergence error:",nls$message)) }
		cdef <- nls$fullpar
		cdef[1:5] <- log(cdef[1:5])
		cdef[1:7] <- cdef[1:7]+runif(7,min=-0.01,max=0.01)
		# cdef[8:10] <- mean(cdef[8:10])+runif(1,min=-0.01,max=0.01)
		cdef[1:5] <- exp(cdef[1:5])
		cfixed <- nls$fixed
		if ((is.character(fixed)||length(fixed)<10) && !is.na(cfixed[10])) { cfixed[10] <- cdef[10] }
		cdef <- cdef[which(is.na(cfixed))]
		cdef <- pmax(cdef,nls$mlims[1,])
		cdef <- pmin(cdef,nls$mlims[2,])
	}
	ostart <- nls$odef
	nfixed <- nls$fixed
	fullpar <- nls$fullpar
	nfixed[which(!is.na(nfixed))] <- fullpar[which(!is.na(nfixed))]
	ccoef <- fullpar[which(is.na(nfixed))]
	pred_act <- evalBRAIDrsm(conc1,conc2,fullpar)
	bfit <- list(conc1=conc1,conc2=conc2,act=act,fitted.values=pred_act,residuals=(pred_act-act),
					ostart=ostart,fixed=nfixed,mlims=nls$mlims,coefficients=ccoef,fullpar=fullpar,
					convergence=nls$convergence,message=nls$message,call=match.call())
	
	cfnames <- c("IDMA","IDMB","na","nb","delta","kappa","E0","EfA","EfB","EfAB")
	names(bfit$coefficients) <- cfnames[which(is.na(nfixed))]
	class(bfit) <- "braidrm"
	if (getCIs) { bfit <- getBRAIDbootstrap(bfit) }
	return(bfit)
}
braidrm.formula <- function(model,data,...) {
	mf <- model.frame(formula=model, data=data)
	concs <- model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- model.response(mf)
	bfit <- braidrm.default(concs,act,...)
	bfit$call <- match.call()
	return(bfit)
}
getBRAIDbootstrap <- function(bfit,ciLevs=c(0.025,0.975),numBoot=NULL) {
	if (class(bfit)!="braidrm") { stop("Input 'bfit' must be of class 'braidrm'.") }
	if (!is.null(bfit$ciPass)) {
		warning("Input already has bootstrapped coefficients (or has failed bootstrapping).")
		return(bfit)
	}
	prcBootPass <- 0.5
	if (is.null(numBoot)) { numBoot <- max(min(10/(1-ciLevs[2]+ciLevs[1]),1000),100) }
	bCoefs <- array(coef(bfit),dim=c(1,length(coef(bfit))))
	for (i in 1:numBoot) {
		bact <- fitted(bfit) + sample(resid(bfit),length(resid(bfit)),replace=TRUE)
		nls <- try(fitBRAIDrsm(bfit$conc1,bfit$conc2,bact,fixed=bfit$fixed,
							def=coef(bfit),llims=bfit$mlims[1,],ulims=bfit$mlims[2,]),silent=TRUE)
		if (class(nls)!="try-error" && nls$convergence==0 && !(TRUE %in% (is.infinite(nls$par) | is.nan(nls$par)))) {
			cpar <- nls$fullpar[which(is.na(bfit$fixed))]
			bCoefs <- rbind(bCoefs,array(cpar,dim=c(1,length(cpar))))
		}
	}
	nbfit <- bfit
	if (nrow(bCoefs)<numBoot*prcBootPass) {
		nbfit$ciPass <- FALSE
	} else {
		qMat <- apply(bCoefs,2,quantile,probs=ciLevs)
		nbfit <- c(nbfit,list(ciPass=TRUE,ciLevs=ciLevs,bCoefs=bCoefs,ciVec=as.vector(qMat)))
	}
	class(nbfit) <- "braidrm"
	return(nbfit)
}
calcBRAIDconfint <- function(bfit,parfunc,civals=NULL) {
	if (class(bfit)!="braidrm") { stop("Input 'bfit' must be of class 'braidrm'.") }
	if (class(parfunc)!="function") { stop("Input 'parfunc' must be a function of one variable.") }
	if (is.null(bfit$ciPass) || !bfit$ciPass) { stop("Input 'bfit' must have bootstrapped coefficients.") }
	if (is.null(civals)) { civals <- bfit$ciLevs }
	
	ostart <- bfit$ostart
	fixed <- bfit$fixed
	fullpar <- bfit$fullpar
	isIncr <- sign(fullpar[10]-fullpar[7])
	parv2fullpar <- function(parv) {
		tfullpar <- fullpar
		tfullpar[which(is.na(fixed))] <- parv
		if (is.na(fixed[10])) {
			if (is.na(fixed[8])||is.na(fixed[9])) {
				if (!is.na(fixed[8]) && ostart[8]==ostart[10]) { tfullpar[8] <- tfullpar[10] }
				else if (!is.na(fixed[9]) && ostart[9]==ostart[10]) { tfullpar[9] <- tfullpar[10] }
				else { tfullpar[10] <- tfullpar[10]+isIncr*max(isIncr*tfullpar[8:9]) }
			}
			else if (ostart[8]==ostart[10] && ostart[9]==ostart[10]) { tfullpar[8:9] <- tfullpar[10] }
		} else if ((is.na(fixed[8]) || is.na(fixed[9])) && ostart[10]==isIncr*max(isIncr*ostart[8:9])) {
			tfullpar[10] <- isIncr*max(isIncr*tfullpar[8:9])
		}
		return(tfullpar)
	}
	
	outval <- parfunc(fullpar)
	outmat <- outval
	if (length(outmat)==1) { outmat <- array(outmat,dim=c(1,1)) }
	for (b in 1:nrow(bfit$bCoefs)) { outmat <- cbind(outmat,parfunc(parv2fullpar(bfit$bCoefs[b,]))) }
	outci <- apply(outmat,1,quantile,probs=civals)
	fullout <- cbind(as.vector(outci[1,]),outval,as.vector(outci[2,]))
	return(fullout)
}

# Params 1 and 2: IDMA and IDMB
# Params 3 and 4: na and nb
# Param 5: delta
# Param 6: kappa
# Param 7: E0
# Params 8 and 9: EA and EB
# Param 10: EAB
evalBRAIDrsm <- function(DA,DB,parv) { return(evalBRAIDrsm_int(DA,DB,parv)) }
evalBRAIDrsm_int <- function(DA,DB,parv,fixed=c(NA,NA,NA,NA,1,NA,NA,NA,NA,NA),calcderivs=FALSE) {
	concA <- DA
	concB <- DB
	if (abs(parv[8]-parv[7])>abs(parv[10]-parv[7]) || abs(parv[9]-parv[7])>abs(parv[10]-parv[7])) {
		stop("Invalid parameter values.")
	}
	if (parv[8]==parv[10]) {
		Aden <- 1
		Amod <- (concA/parv[1])^(sqrt(parv[3]/parv[4])/parv[5])
	} else {
		Aden <- 1+(1-((parv[8]-parv[7])/(parv[10]-parv[7])))*((concA/parv[1])^parv[3])
		Amod <- ((parv[8]-parv[7])/(parv[10]-parv[7]))*((concA/parv[1])^parv[3])
		iinds <- which(is.infinite(Aden)|is.infinite(Amod))
		Aden[iinds] <- 1-(parv[8]-parv[7])/(parv[10]-parv[7])
		Amod[iinds] <- (parv[8]-parv[7])/(parv[10]-parv[7])
		Amod <- (Amod/Aden)^(1/(parv[5]*sqrt(parv[3]*parv[4])))
		Aden[iinds] <- Inf
	}
	if (parv[9]==parv[10]) {
		Bden <- 1
		Bmod <- (concB/parv[2])^(sqrt(parv[4]/parv[3])/parv[5])
	} else {
		Bden <- 1+(1-((parv[9]-parv[7])/(parv[10]-parv[7])))*((concB/parv[2])^parv[4])
		Bmod <- ((parv[9]-parv[7])/(parv[10]-parv[7]))*((concB/parv[2])^parv[4])
		iinds <- which(is.infinite(Bden)|is.infinite(Bmod))
		Bden[iinds] <- 1-(parv[9]-parv[7])/(parv[10]-parv[7])
		Bmod[iinds] <- (parv[9]-parv[7])/(parv[10]-parv[7])
		Bmod <- (Bmod/Bden)^(1/(parv[5]*sqrt(parv[3]*parv[4])))
		Bden[iinds] <- Inf
	}
	ABrad <- sqrt(Amod*Bmod)
	if (parv[6]!=0) { ABrad[which(Amod==0|Bmod==0)] <- 0
	} else { ABrad[which(Amod==0|Bmod==0|ABrad==Inf)] <- 0 }
	Exp <- Amod+Bmod+parv[6]*ABrad
	Exp[is.infinite(Amod)|is.infinite(Bmod)] <- Inf
	Den <- 1+Exp^(-parv[5]*sqrt(parv[3]*parv[4]))
	Ei <- parv[7]+(parv[10]-parv[7])/Den
	if (!calcderivs) { return(Ei) }
	Emod <- (parv[10]-parv[7])*(Exp^(-parv[5]*sqrt(parv[3]*parv[4])-1))/(Den^2)
	Emod[which(Exp==0)] <- 0
	Emodz <- which(Emod==0 | is.infinite(Emod))
	derivs <- array(0,dim=c(length(Ei),10))
	if (is.na(fixed[1])) { derivs[,1] <- -parv[3]*Emod*(Amod+parv[6]*ABrad/2)/(Aden*parv[1]) }
	if (is.na(fixed[2])) { derivs[,2] <- -parv[4]*Emod*(Bmod+parv[6]*ABrad/2)/(Bden*parv[2]) }
	if (is.na(fixed[3]) || is.na(fixed[4]) || is.na(fixed[5])) {
		AlogA <- Amod*log(Amod)
		AlogA[which(Amod==0)] <- 0
		BlogB <- Bmod*log(Bmod)
		BlogB[which(Bmod==0)] <- 0
		RlogR <- ABrad*log(ABrad)
		RlogR[which(ABrad==0)] <- 0
		ElogE <- Exp*log(Exp)
		ElogE[which(Exp==0)] <- 0
		derivs[,5] <- -sqrt(parv[3]*parv[4])*Emod*(AlogA+BlogB+parv[6]*RlogR-ElogE)
		if (is.na(fixed[3])) {
			derivs[,3] <- Emod*(Amod+parv[6]*ABrad/2)*log(concA/parv[1])/Aden
			derivs[which(concA==0),3] <- 0
			derivs[,3] <- derivs[,3]+parv[5]*derivs[,5]/(2*parv[3])
		}
		if (is.na(fixed[4])) {
			derivs[,4] <- Emod*(Bmod+parv[6]*ABrad/2)*log(concB/parv[2])/Bden
			derivs[which(concB==0),4] <- 0
			derivs[,4] <- derivs[,4]+parv[5]*derivs[,5]/(2*parv[4])
		}
	}
	if (is.na(fixed[6])) { derivs[,6] <- parv[5]*sqrt(parv[3]*parv[4])*Emod*ABrad }
	derivs[Emodz,1:6] <- 0
	if (is.na(fixed[7]) || is.na(fixed[8]) || is.na(fixed[9]) || is.na(fixed[10])) {
		if (is.na(fixed[10]) && !is.na(fixed[8]) && !is.na(fixed[9]) && parv[8]==parv[10] && parv[9]==parv[10]) {
			derivs[,10] <- 1/Den
		} else if (xor(is.na(fixed[8]),is.na(fixed[10])) && parv[8]==parv[10]) {
			derivs[,9] <- Emod*(Bmod+parv[6]*ABrad/2)*(1+(concB/parv[2])^parv[4])/(Bden*(parv[9]-parv[7]))
			derivs[Emodz,9] <- 0
			if (is.na(fixed[8])) { derivs[,8] <- 1/Den-derivs[,9]*(parv[9]-parv[7])/(parv[8]-parv[7]) }
			else { derivs[,10] <- 1/Den-derivs[,9]*(parv[9]-parv[7])/(parv[10]-parv[7]) }
		} else if (xor(is.na(fixed[9]),is.na(fixed[10])) && parv[9]==parv[10]) {
			derivs[,8] <- Emod*(Amod+parv[6]*ABrad/2)*(1+(concA/parv[1])^parv[3])/(Aden*(parv[8]-parv[7]))
			derivs[Emodz,8] <- 0
			if (is.na(fixed[9])) { derivs[,9] <- 1/Den-derivs[,8]*(parv[8]-parv[7])/(parv[9]-parv[7]) }
			else { derivs[,10] <- 1/Den-derivs[,8]*(parv[8]-parv[7])/(parv[10]-parv[7]) }
		} else {
			derivs[,8] <- Emod*(Amod+parv[6]*ABrad/2)*(1+(concA/parv[1])^parv[3])/(Aden*(parv[8]-parv[7]))
			derivs[,9] <- Emod*(Bmod+parv[6]*ABrad/2)*(1+(concB/parv[2])^parv[4])/(Bden*(parv[9]-parv[7]))
			derivs[Emodz,8:9] <- 0
			derivs[,10] <- 1/Den-derivs[,8]*(parv[8]-parv[7])/(parv[10]-parv[7])-derivs[,9]*(parv[9]-parv[7])/(parv[10]-parv[7])
		}
		if (is.na(fixed[7])) { derivs[,7] <- 1-derivs[,8]-derivs[,9]-derivs[,10] }
	}
	return(cbind(Ei,derivs))
}
invertBRAIDrsm <- function(val,DA=NULL,DB=NULL,parv) {
	if ((is.null(DA)&&is.null(DB))||(!is.null(DA)&&!is.null(DB))) {
		stop("Exactly one of the inputs 'DA' and 'DB' must be NULL (the default).")
	}
	conc1 <- DA
	conc2 <- DB
	if (is.null(conc1)) {
		if (length(val)!=length(conc2)) {
			if (length(val)%%length(conc2)==0) { conc2 <- rep(conc2,times=length(val)/length(conc2)) }
			else if (length(conc2)%%length(val)==0) { val <- rep(val,times=length(conc2)/length(val)) }
			else { stop("The length of the longer of inputs 'val' and 'conc2' is not a multiple of the length of the shorter.") }
		}
	} else {
		if (length(val)!=length(conc1)) {
			if (length(val)%%length(conc1)==0) { conc1 <- rep(conc1,times=length(val)/length(conc1)) }
			else if (length(conc1)%%length(val)==0) { val <- rep(val,times=length(conc1)/length(val)) }
			else { stop("The length of the longer of inputs 'val' and 'conc1' is not a multiple of the length of the shorter.") }
		}
	}
	concout <- rep(0,times=length(val))
	incr <- sign(parv[10]-parv[7])
	Eval <- (val-parv[7])/(parv[10]-val)
	concout[which(incr*(val-parv[7])<=0)] <- 0
	concout[which(incr*(val-parv[10])>=0)] <- Inf
	cinds <- which(Eval>0)
	Eval <- Eval[cinds]^(1/(parv[5]*sqrt(parv[3]*parv[4])))
	if (is.null(conc1)) {
		if (parv[9]==parv[10]) {
			Bmod <- (conc2[cinds]/parv[2])^(sqrt(parv[4]/parv[3])/parv[5])
		} else {
			Bden <- 1+(1-((parv[9]-parv[7])/(parv[10]-parv[7])))*((conc2[cinds]/parv[2])^parv[4])
			Bmod <- ((parv[9]-parv[7])/(parv[10]-parv[7]))*((conc2[cinds]/parv[2])^parv[4])
			iinds <- which(is.infinite(Bden)|is.infinite(Bmod))
			Bden[iinds] <- 1-(parv[9]-parv[7])/(parv[10]-parv[7])
			Bmod[iinds] <- (parv[9]-parv[7])/(parv[10]-parv[7])
			Bmod <- (Bmod/Bden)^(1/(parv[5]*sqrt(parv[3]*parv[4])))
		}
		if (parv[6]>=0) {
			concout[cinds[which(Bmod>=Eval)]] <- 0
			ccinds <- which(Bmod<Eval)
		} else {
			concout[cinds[which(Bmod>=Eval/(1-(parv[6]^2)/4))]] <- 0
			ccinds <- which(Bmod<Eval/(1-(parv[6]^2)/4))
		}
		Bmod <- Bmod[ccinds]
		Eval <- Eval[ccinds]
		cinds <- cinds[ccinds]
		Jval <- Eval+((parv[6]^2)/2-1)*Bmod
		if (parv[6]!=0) { Jval <- Jval-parv[6]*sqrt(((parv[6]^2)/4-1)*Bmod^2+Eval*Bmod) }
		Jval <- Jval^(parv[5]*sqrt(parv[3]*parv[4]))
		if (parv[8]!=parv[10]) {
			concout[cinds[which(Jval>=(parv[8]-parv[7])/(parv[10]-parv[8]))]] <- Inf
			cinds <- cinds[which(Jval<(parv[8]-parv[7])/(parv[10]-parv[8]))]
			Jval <- Jval[which(Jval<(parv[8]-parv[7])/(parv[10]-parv[8]))]
		}
		conc1c <- parv[1]*(Jval/((1+Jval)*(parv[8]-parv[7])/(parv[10]-parv[7])-Jval))^(1/parv[3])
		concout[cinds] <- conc1c
	} else {
		if (parv[8]==parv[10]) {
			Amod <- (conc1/parv[1])^(sqrt(parv[3]/parv[4])/parv[5])
		} else {
			Aden <- 1+(1-((parv[8]-parv[7])/(parv[10]-parv[7])))*((conc1/parv[1])^parv[3])
			Amod <- ((parv[8]-parv[7])/(parv[10]-parv[7]))*((conc1/parv[1])^parv[3])
			iinds <- which(is.infinite(Aden)|is.infinite(Amod))
			Aden[iinds] <- 1-(parv[8]-parv[7])/(parv[10]-parv[7])
			Amod[iinds] <- (parv[8]-parv[7])/(parv[10]-parv[7])
			Amod <- (Amod/Aden)^(1/(parv[5]*sqrt(parv[3]*parv[4])))
		}
		if (parv[6]>=0) {
			concout[cinds[which(Amod>=Eval)]] <- 0
			ccinds <- which(Amod<Eval)
		} else {
			concout[cinds[which(Amod>=Eval/(1-(parv[6]^2)/4))]] <- 0
			ccinds <- which(Amod<Eval/(1-(parv[6]^2)/4))
		}
		Amod <- Amod[ccinds]
		Eval <- Eval[ccinds]
		cinds <- cinds[ccinds]
		Jval <- Eval+((parv[6]^2)/2-1)*Amod
		if (parv[6]!=0) { Jval <- Jval-parv[6]*sqrt(((parv[6]^2)/4-1)*Amod^2+Eval*Amod) }
		Jval <- Jval^(parv[5]*sqrt(parv[3]*parv[4]))
		if (parv[9]!=parv[10]) {
			concout[cinds[which(Jval>=(parv[9]-parv[7])/(parv[10]-parv[9]))]] <- Inf
			cinds <- cinds[which(Jval<(parv[9]-parv[7])/(parv[10]-parv[9]))]
			Jval <- Jval[which(Jval<(parv[9]-parv[7])/(parv[10]-parv[9]))]
		}
		conc2c <- parv[2]*(Jval/((1+Jval)*(parv[9]-parv[7])/(parv[10]-parv[7])-Jval))^(1/parv[4])
		concout[cinds] <- conc2c
	}
	return(concout)
}

# (...,NA,NA,NA)                : All are varying, EAB >= EA,EB
# (...,NA,NA,~) EAB=rmax(EA,EB) : EA and EB varying, EAB varies with maximum
# (...,NA,NA,~) EAB>rmax(EA,EB) : EA and EB varying below fixed EAB
# (...,NA,~,NA) EAB=EB>=EA      : EA varying below EAB and EB varying together
# (...,NA,~,NA) EAB>EB;EAB>=EA  : EA varying, EB fixed EAB varying above EA and EB
# (...,~,~,NA) EAB=EA=EB        : All three varying together
# (...,~,~,NA) EAB=EA>EB        : EAB and EA varying together above EB
# (...,~,~,NA) EAB>EA & EAB>EB  : EAB varying above fixed EA and EB
# (...,NA,~,~) EAB=rmax(EA,EB)  : EA varying, EB fixed, EAB varies with maximum
# (...,NA,~,~) EAB>EA;EAB>EB    : EA varying below fixed EAB, EB fixed
# (...,~,~,~)                   : All three fixed
fitBRAIDrsm <- function(conc1,conc2,act,def=NULL,fixed=NULL,llims=NULL,ulims=NULL) {
	# Fill in 'fixed' and 'def'
	if (is.null(fixed)) { tfixed <- c(1,2,3,4,6,7,8,9) }
	else { tfixed <- fixed }
	if (length(tfixed)!=10 || length(which(is.na(tfixed)))==0) {
		if (is.character(tfixed)) {
			if (tfixed=="kappa1") { fixinds <- c(1,2,3,4,6,7,10) }
			else if (tfixed=="kappa2") { fixinds <- c(1,2,3,4,6,7,8,9) }
			else if (tfixed=="delta1") { fixinds <- c(1,2,3,4,5,7,10) }
			else if (tfixed=="delta2") { fixinds <- c(1,2,3,4,5,7,8,9) }
			else if (tfixed=="kappa3") { fixinds <- c(1,2,3,4,6,7,8,9,10) }
			else if (tfixed=="delta3") { fixinds <- c(1,2,3,4,5,7,8,9,10) }
			else if (tfixed=="ebraid") { fixinds <- 1:10 }
			else { stop(sprintf("Unknown model name '%s'.",tfixed)) }
		} else {
			if (length(which(is.na(tfixed)))>0) { stop("'NA' values are only permitted in raw 'fixed' vectors.") }
			fixinds <- sort(unique(round(tfixed)))
			if (min(fixinds)<1 || max(fixinds)>10) { stop("Indexed 'fixed' vectors may only contain values from 1 to 10.") }
		}
	} else { fixinds <- which(is.na(tfixed)) }
	npar <- length(fixinds)
	if (npar==0) { stop("At least one parameter must be allowed to vary.") }
	if (is.null(def)) {
		inds <- which(conc1+conc2>0)
		sdef <- dr_start(conc1[inds]+conc2[inds],act[inds])
		odef <- c(exp(sdef[4]),exp(sdef[4]),sdef[1],sdef[1],1,0,sdef[2],sdef[3],sdef[3],sdef[3])
		if (length(tfixed)!=10 || length(which(is.na(tfixed)))==0) {
			tfixed <- odef
			tfixed[fixinds] <- NA
		} else {
			odef[which(!is.na(tfixed))] <- tfixed[which(!is.na(tfixed))]
			if (odef[10]>odef[7] && !is.na(tfixed[10])) { odef[8:10] <- pmin(odef[8:10],tfixed[10]) }
			else if (!is.na(tfixed[10])) { odef[8:10] <- pmax(odef[8:10],tfixed[10]) }
		}
	} else {
		if (length(tfixed)!=10 || length(which(is.na(tfixed)))==0) {
			if (length(def)==10) {
				odef <- def
				tfixed <- odef
				tfixed[fixinds] <- NA
			} else { stop("When using index vectors or model names for 'fixed', parameter 'def' must specify all 10 values.") }
		} else {
			if (length(def)==10) {
				odef <- def
				tfixed <- odef
				tfixed[fixinds] <- NA
			} else if (length(def)==npar) {
				odef <- tfixed
				odef[fixinds] <- def
			} else { stop("Input 'def' must have as many values as free parameters or all parameter values (length 10).") }
		}
	}
	fixed <- tfixed
	isIncr <- sign(odef[10]-odef[7])
	# medAct <- isIncr*min(max(isIncr*median(act),isIncr*odef[7]),min(isIncr*odef[8:10]))
	medAct <- isIncr*(isIncr*odef[7]+min(isIncr*odef[8:10]))/2
	
	# Fill in lower limit values
	if (is.null(llims)) { llims <- rep(NA,times=10) }
	else {
		if (length(llims)<10 && length(llims)!=npar) { stop("Improper length for lower limit vector.") }
		else if (length(llims)<10) {
			tllims <- rep(NA,times=10)
			tllims[which(is.na(fixed))] <- llims
			llims <- tllims
		}
	}
	full_llims <- c(min(exp(-12),odef[1]/exp(3)),min(exp(-12),odef[2]/exp(3)),1/10,1/10,1/10,-1.999,0,0,0,0)
	if (isIncr>0) { full_llims[7:8] <- c(2*odef[7]-odef[10],medAct+0.001) }
	else if (isIncr<0) { full_llims[7:8] <- c(medAct+0.001,2*odef[10]-odef[7]) }
	else { stop("Default limits must cover a non-zero width range.") }
	full_llims[9:10] <- full_llims[8]
	full_llims[which(!is.na(llims))] <- pmax(full_llims[which(!is.na(llims))],llims[which(!is.na(llims))])
	full_llims[1] <- max(full_llims[1],exp(2*min(log(conc1[conc1>0]))-max(log(conc1[conc1>0]))))
	full_llims[2] <- max(full_llims[2],exp(2*min(log(conc2[conc2>0]))-max(log(conc2[conc2>0]))))
	odef[which(is.na(fixed))] <- pmax(odef,full_llims)[which(is.na(fixed))]
	
	# Fill in upper limit values
	if (is.null(ulims)) { ulims <- rep(NA,times=10) }
	else {
		if (length(ulims)<10 && length(ulims)!=npar) { stop("Improper length for upper limit vector.") }
		else if (length(ulims)<10) {
			tulims <- rep(NA,times=10)
			tulims[which(is.na(fixed))] <- ulims
			ulims <- tulims
		}
	}
	full_ulims <- c(max(1,odef[1]*exp(3)),max(1,odef[2]*exp(3)),10,10,10,100,0,0,0,0)
	if (isIncr>0) { full_ulims[7:8] <- c(medAct-0.001,2*odef[10]-odef[7]) }
	else { full_ulims[7:8] <- c(2*odef[7]-odef[10],medAct-0.001) }
	full_ulims[9:10] <- full_ulims[8]
	full_ulims[which(!is.na(ulims))] <- pmin(full_ulims[which(!is.na(ulims))],ulims[which(!is.na(ulims))])
	full_ulims[1] <- min(full_ulims[1],exp(2*max(log(conc1[conc1>0]))-min(log(conc1[conc1>0]))))
	full_ulims[2] <- min(full_ulims[2],exp(2*max(log(conc2[conc2>0]))-min(log(conc2[conc2>0]))))
	odef[which(is.na(fixed))] <- pmin(odef,full_ulims)[which(is.na(fixed))]
	
	# Prepare parameter values and bounds for optimization
	tdef <- odef
	tdef[1:5] <- log(tdef[1:5])
	full_llims[1:5] <- log(full_llims[1:5])
	full_ulims[1:5] <- log(full_ulims[1:5])
	if (is.na(fixed[10])) {
		if (is.na(fixed[8])||is.na(fixed[9])) {
			if (!is.na(fixed[8]) && odef[8]==odef[10]) {
				tdef[10] <- tdef[10]-tdef[9]
				if (isIncr>0) {
					full_llims[10] <- 0
					full_ulims[10] <- full_ulims[10]-tdef[9]
				} else {
					full_llims[10] <- full_llims[10]-tdef[7]
					full_ulims[10] <- 0
				}
			} else if (!is.na(fixed[9]) && odef[9]==odef[10]) {
				tdef[10] <- tdef[10]-tdef[8]
				if (isIncr>0) {
					full_llims[10] <- 0
					full_ulims[10] <- full_ulims[10]-tdef[8]
				} else {
					full_llims[10] <- full_llims[10]-tdef[8]
					full_ulims[10] <- 0
				}
			} else {
				tdef[10] <- tdef[10]-isIncr*max(isIncr*tdef[8:9])
				if (isIncr>0) {
					full_llims[10] <- 0
					full_ulims[10] <- full_ulims[10]-isIncr*max(isIncr*tdef[8:9])
				} else {
					full_llims[10] <- full_llims[10]-isIncr*max(isIncr*tdef[8:9])
					full_ulims[10] <- 0
				}
			}
		} else {
			if (odef[8]!=odef[10]) {
				if (isIncr>0) { full_llims[10] <- max(full_llims[10],odef[8]) }
				else { full_ulims[10] <- min(full_ulims[10],odef[8]) }
			}
			if (odef[9]!=odef[10]) {
				if (isIncr>0) { full_llims[10] <- max(full_llims[10],odef[9]) }
				else { full_ulims[10] <- min(full_ulims[10],odef[9]) }
			}
		}
	} else if (isIncr*odef[10]>max(isIncr*odef[8:9])) {
		if (isIncr>0) { full_ulims[8:9] <- pmin(full_ulims[8:9],odef[10]) }
		else { full_llims[8:9] <- pmax(full_llims[8:9],odef[10]) }
	}
	llims <- full_llims[which(is.na(fixed))]
	ulims <- full_ulims[which(is.na(fixed))]
	
	parv2fullpar <- function(parv) {
		fullpar <- tdef
		fullpar[which(is.na(fixed))] <- parv
		fullpar[1:5] <- exp(fullpar[1:5])
		if (is.na(fixed[10])) {
			if (is.na(fixed[8])||is.na(fixed[9])) {
				if (!is.na(fixed[8]) && odef[8]==odef[10]) { 
					fullpar[10] <- fullpar[10]+fullpar[9]
					fullpar[8] <- fullpar[10]
				} else if (!is.na(fixed[9]) && odef[9]==odef[10]) {
					fullpar[10] <- fullpar[10]+fullpar[8]
					fullpar[9] <- fullpar[10]
				} else { fullpar[10] <- fullpar[10]+isIncr*max(isIncr*fullpar[8:9]) }
			} else {
				if (odef[8]==odef[10]) { fullpar[8] <- fullpar[10] }
				if (odef[9]==odef[10]) { fullpar[9] <- fullpar[10] }
			}
		} else if (is.na(fixed[8]) && is.na(fixed[9])) {
			if (odef[10]==isIncr*max(isIncr*odef[8:9])) { fullpar[10] <- isIncr*max(isIncr*fullpar[8:9]) }
		} else if (is.na(fixed[8]) && ((odef[8]==odef[10]&&isIncr*fullpar[8]>=isIncr*fullpar[9])
							|| isIncr*fullpar[8]>isIncr*odef[10])) { fullpar[10] <- fullpar[8] }
		else if (is.na(fixed[9]) && ((odef[9]==odef[10]&&isIncr*fullpar[9]>=isIncr*fullpar[8])
							|| isIncr*fullpar[9]>isIncr*odef[10])) { fullpar[10] <- fullpar[9] }
		return(fullpar)
	}
	parvRSMErr <- function(parv) {
		fullpar <- parv2fullpar(parv)
		ei <- evalBRAIDrsm(conc1,conc2,fullpar)
		err <- ei-act
		return(sum(err^2))
	}
	parvRSMDerivs <- function(parv) {
		fullpar <- parv2fullpar(parv)
		ederivs <- evalBRAIDrsm_int(conc1,conc2,fullpar,fixed=fixed,calcderivs=TRUE)
		err <- ederivs[,1]-act
		for (i in 1:5) { ederivs[,i+1] <- fullpar[i]*ederivs[,i+1] }
		if (is.na(fixed[10]) && (is.na(fixed[8])||is.na(fixed[9]))) {
			if (!is.na(fixed[8]) && odef[8]==odef[10]) { ederivs[,11] <- ederivs[,11]-ederivs[,10] }
			else if (!is.na(fixed[9]) && odef[9]==odef[10]) { ederivs[,11] <- ederivs[,11]-ederivs[,9] }
			else {
				if (isIncr*fullpar[8]>=isIncr*fullpar[9]) { ederivs[,11] <- ederivs[,11]-ederivs[,9] }
				else { ederivs[,11] <- ederivs[,11]-ederivs[,10] }
			}
		}
		oderivs <- rep(0,times=10)
		for (i in 1:10) { oderivs[i] <- 2*sum(err*ederivs[,i+1]) }
		oderivs <- oderivs[which(is.na(fixed))]
		return(oderivs)
	}
	
	escale <- max(abs(odef[5]-odef[9]),abs(odef[6]-odef[9]),abs(odef[10]-odef[9]))
	parscale <- c(1,1,1,1,1,1,0.1*escale,0.1*escale,0.1*escale,0.1*escale)
	control <- list(maxit=500,parscale=parscale[which(is.na(fixed))])
	nls <- optim(tdef[which(is.na(fixed))],parvRSMErr,parvRSMDerivs,method="L-BFGS-B",
					lower=llims,upper=ulims,hessian=FALSE,control=control)
	nls$fullpar <- parv2fullpar(nls$par)
	nls$fixed <- fixed
	nls$odef <- odef
	full_llims <- parv2fullpar(llims)
	llims <- full_llims[which(is.na(fixed))]
	full_ulims <- parv2fullpar(ulims)
	ulims <- full_ulims[which(is.na(fixed))]
	nls$mlims <- rbind(llims,ulims)
	return(nls)
}
