findBestBRAID <- function(model,data,defaults,startparv=NULL,llims=NULL,ulims=NULL,itype=1,getCIs=TRUE,crossval=TRUE,...) UseMethod("findBestBRAID")

findBestBRAID.default <- function(model,data,defaults,startparv=NULL,llims=NULL,ulims=NULL,itype=1,getCIs=TRUE,crossval=TRUE,...) {
	concs <- model
	act <- data
	if (ncol(concs)<2) { stop("Argument 'concs' must have at least two columns.") }
	
	# Construct a basic set of default values (if needed), and constrain defaults to lie
	# between any specified bounds
	if (is.null(startparv)) {
		inds <- which(concs[,1]+concs[,2]>0)
		sdef <- dr_start(concs[inds,1]+concs[inds,2],act[inds])
		syndef <- c(exp(sdef[4]),exp(sdef[4]),sdef[1],sdef[1],1,0,sdef[2],sdef[3],sdef[3],sdef[3])
	} else { syndef <- startparv }
	if (!is.null(llims)) {
		syndef[!is.na(llims)] <- pmax(syndef,llims)[!is.na(llims)]
		if (!is.null(startparv)) { startparv[!is.na(llims)] <- pmax(startparv,llims)[!is.na(llims)] }
	}
	if (!is.null(ulims)) {
		syndef[!is.na(ulims)] <- pmin(syndef,ulims)[!is.na(ulims)]
		if (!is.null(startparv)) { startparv[!is.na(ulims)] <- pmin(startparv,ulims)[!is.na(ulims)] }
	}

	# Construct array of 'fixed' vectors representing all models to be tested
	synfix <- c(1,1,1,1,1,0,0,1,1,1)
	synfix[7] <- defaults[1]
	synfix[8:10] <- defaults[2]
	fixarr <- array(synfix,dim=c(10,10))
	if (itype==1) { fixarr[c(1:4,6),] <- NA
	} else { fixarr[1:5,] <- NA }
	fixarr[7,c(1,3:6)] <- NA
	fixarr[8,c(1:3,7)] <- NA
	fixarr[9,c(1,2,4,8)] <- NA
	if (!is.null(startparv)) {
		fixarr[c(8,9),c(5,9)] <- startparv[10]
		if (startparv[10]>startparv[7]) {
			fixarr[10,c(1,2)] <- max(startparv[8:9])
			fixarr[10,c(3,7)] <- max(startparv[8],defaults[2])
			fixarr[10,c(4,8)] <- max(startparv[9],defaults[2])
		} else {
			fixarr[10,c(1,2)] <- min(startparv[8:9])
			fixarr[10,c(3,7)] <- min(startparv[8],defaults[2])
			fixarr[10,c(4,8)] <- min(startparv[9],defaults[2])
		}
	}
	fixarr[10,c(5,9)] <- NA
	allmods <- TRUE
	if (itype==3 | itype==4) {
		tfixarr <- fixarr
		tfixarr[5,] <- 1
		tfixarr[6,] <- NA
		fixarr <- cbind(fixarr,tfixarr)
		if (itype==4) {
			tfixarr[5,] <- NA
			fixarr <- cbind(fixarr,tfixarr)
		}
	} else if (itype==5) {
		fixarr[5,] <- 1
	}
	
	if (crossval) {
		# Randomly assign all data points to one of four blocks for cross-validation
		ptclass <- rep(0,times=nrow(concs))
		ninds <- which(concs[,1]==0 & concs[,2]==0)
		if (length(ninds)!=0) { ptclass[ninds] <- sample(length(ninds)) }
		linds <- which(concs[,1]==0 & concs[,2]!=0)
		if (length(linds)!=0) { ptclass[linds] <- sample(length(linds))+max(ptclass) }
		pinds <- which(concs[,1]!=0 & concs[,2]==0)
		if (length(pinds)!=0) { ptclass[pinds] <- sample(length(pinds))+max(ptclass) }
		cinds <- which(concs[,1]!=0 & concs[,2]!=0)
		if (length(cinds)!=0) { ptclass[cinds] <- sample(length(cinds))+max(ptclass) }
		else if (itype<5) { stop("No positive combinations were found! Optimization is ill-defined except when 'itype' equals 5 (no interaction).") }
		ptclass <- ptclass%%4
	}
	
	modK <- modAIC <- rep(Inf,times=ncol(fixarr))
	modpars <- fixarr
	for (ii in 1:ncol(fixarr)) {
		modK[ii] <- length(which(is.na(fixarr[,ii])))
		
		# Fit best BRAID parameters for the given model
		if (!is.null(startparv)) { cstart <- startparv[is.na(fixarr[,ii])] }
		else { cstart <- NULL }
		if (!is.null(llims)) { cllims <- llims[is.na(fixarr[,ii])] }
		else { cllims <- NULL }
		if (!is.null(ulims)) { culims <- ulims[is.na(fixarr[,ii])] }
		else { culims <- NULL }
		bfit <- try(braidrm(concs,act,getCIs=FALSE,fixed=fixarr[,ii],startparv=cstart,llims=cllims,ulims=culims),silent=TRUE)
		if (class(bfit)=='braidrm' && !any((bfit$fullpar[8:10]-bfit$fullpar[7])*(defaults[2]-defaults[1])<0)) {
			# Store best fit parameters for later optimization
			modpars[which(is.na(fixarr[,ii])),ii] <- coef(bfit)
			if ((ii%%10) %in% c(1:4,7,8)) {
				if (syndef[10]>syndef[7]) { modpars[10,ii] <- max(modpars[8:9,ii]) }
				else { modpars[10,ii] <- min(modpars[8:9,ii]) }
			} else if ((ii%%10) %in% c(5,9)) { modpars[8:9,ii] <- modpars[10,ii] }
			cstart <- modpars[,ii]
			sumerr <- sum(resid(bfit)^2)
		} else {
			# If fit unsuccessful, fill parameter column with default values
			modpars[which(is.na(fixarr[,ii])),ii] <- syndef[which(is.na(fixarr[,ii]))]
			sumerr <- NA
		}
		
		if (crossval) {
			# Evaluate model fit for four cross-validation blocks
			sumerr <- 0
			for (jj in 0:3) {
				trn <- ptclass==jj
				tst <- !trn
				bfit <- try(braidrm(concs[trn,],act[trn],getCIs=FALSE,fixed=fixarr[,ii],
									startparv=cstart,llims=cllims,ulims=culims),silent=TRUE)
				if (class(bfit)=='braidrm') {
					sumerr <- sumerr+sum((evalBRAIDrsm(concs[tst,1],concs[tst,2],bfit$fullpar)-act[tst])^2)
				} else {
					sumerr <- NA
					break
				}
			}
		}
		
		# Calculate log-likelihood and AIC
		if (!is.na(sumerr)) {
			df <- nrow(concs)
			llik <- -(df/2)*(log(2*pi)+log(sumerr)-log(df)+1)
			modAIC[ii] <- 2*modK[ii]-2*llik+2*modK[ii]*(modK[ii]+1)/(df-modK[ii]-1)
		}
	}
	
	# Select best model using AIC, and refit, with confidence intervals if desired
	fixarr[!is.na(fixarr)] <- modpars[!is.na(fixarr)]
	bmi <- modelSelect(modAIC,modK,ik=min(modK))
	if (!is.null(llims)) { cllims <- llims[is.na(fixarr[,bmi])] }
	else { cllims <- NULL }
	if (!is.null(ulims)) { culims <- ulims[is.na(fixarr[,bmi])] }
	else { culims <- NULL }
	bfit <- braidrm(concs,act,getCIs=getCIs,fixed=fixarr[,bmi],startparv=modpars[,bmi],llims=cllims,ulims=culims)
	bfit$call <- match.call()
	return(bfit)
}
findBestBRAID.formula <- function(model,data,...) {
	mf <- model.frame(formula=model, data=data)
	concs <- model.matrix(attr(mf, "terms"), data=mf)
	tms <- attr(concs,"assign")
	for (i in seq(length(tms),1,by=-1)) {
		if (tms[i]==0) { concs <- concs[,-i] }
	}
	act <- model.response(mf)
	bfit <- findBestBRAID.default(concs,act,...)
	bfit$call <- match.call()
	return(bfit)
}

runBRAIDanalysis <- function(data,defaults,llims=NULL,ulims=NULL,itype=1,compounds=NULL,corrconc=FALSE,corrsigr=1) {
	if (corrconc && is.null(data$well)) { stop("Concentrations can only be corrected if wells are specified for all measurements.") }
	if (is.null(compounds)) { wdata <- data }
	else {
		wdata <- data[((data$compound1==compounds[1] | data$conc1==0) & (data$compound2==compounds[2] | data$conc2==0)),]
	}
	if (corrconc) {
		wdata$cconc1 <- wdata$conc1
		wdata$cconc2 <- wdata$conc2
	}
	
	ldata <- wdata[wdata$conc2==0,]
	if (nrow(ldata)>0) {
		if (!is.null(llims)) {
			cllims <- llims[c(3,7,8,1)]
			cllims[4] <- log(cllims[4])
		} else { cllims <- NULL }
		if (!is.null(ulims)) {
			culims <- ulims[c(3,7,8,1)]
			culims[4] <- log(culims[4])
		} else { culims <- NULL }
		
		hfit1 <- findBestHill(ldata$conc1,ldata$act,defaults=defaults,llims=cllims,ulims=culims)
		hfit1 <- getHillBootstrap(hfit1)
		if (corrconc) {
			tldata <- data.frame(well=unique(ldata$well[ldata$conc1>0]),conc1=0,act=0)
			for (i in 1:nrow(tldata)) {
				rel <- which(ldata$well==tldata$well[i] & ldata$conc1>0)
				if (length(unique(ldata$conc1[rel]))>1) {
					stop("For concentration corrections, matching wells must contain the same compound concentration.")
				}
				tldata$conc1[i] <- ldata$conc1[rel[1]]
				tldata$act[i] <- mean(ldata$act[rel])
			}
			tldata$cconc1 <- hillConcCorrect(tldata$conc1,tldata$act,coef(hfit1$allfits[[hfit1$bestModIdx]]),sigr=corrsigr)
			for (i in which(wdata$conc1!=0)) {
				rel <- which(tldata$well==wdata$well[i])
				if (length(rel)>0) { wdata$cconc1[i] <- tldata$cconc1[rel[1]] }
			}
		}
	} else { hfit1 <- NULL }
	
	
	
	pdata <- wdata[wdata$conc1==0,]
	if (nrow(ldata)>0) {
		if (!is.null(llims)) {
			cllims <- llims[c(4,7,9,2)]
			cllims[4] <- log(cllims[4])
		} else { cllims <- NULL }
		if (!is.null(ulims)) {
			culims <- ulims[c(4,7,9,2)]
			culims[4] <- log(culims[4])
		} else { culims <- NULL }
		hfit2 <- findBestHill(pdata$conc2,pdata$act,defaults=defaults,llims=cllims,ulims=culims)
		hfit2 <- getHillBootstrap(hfit2)
		if (corrconc) {
			tpdata <- data.frame(well=unique(pdata$well[pdata$conc2>0]),conc2=0,act=0)
			for (i in 1:nrow(tpdata)) {
				rel <- which(pdata$well==tpdata$well[i] & pdata$conc2>0)
				if (length(unique(pdata$conc2[rel]))>1) {
					stop("For concentration corrections, matching wells must contain the same compound concentration.")
				}
				tpdata$conc2[i] <- pdata$conc2[rel[1]]
				tpdata$act[i] <- mean(pdata$act[rel])
			}
			tpdata$cconc2 <- hillConcCorrect(tpdata$conc2,tpdata$act,coef(hfit2$allfits[[hfit2$bestModIdx]]),sigr=corrsigr)
			for (i in which(wdata$conc2!=0)) {
				rel <- which(tpdata$well==wdata$well[i])
				if (length(rel)>0) { wdata$cconc2[i] <- tpdata$cconc2[rel[1]] }
			}
		}
	} else { hfit2 <- NULL }
	
	if (!is.null(hfit1) && !is.null(hfit2)) {
		par1 <- coef(hfit1$allfits[[hfit1$bestModIdx]])
		par2 <- coef(hfit2$allfits[[hfit2$bestModIdx]])
		cstart <- c(exp(par1[4]),exp(par2[4]),par1[1],par2[1],1,0,0,par1[3],par2[3],0)
		if (defaults[1]<defaults[2]) {
			cstart[7] <- min(par1[2],par2[2])
			cstart[10] <- max(par1[3],par2[3])
		} else {
			cstart[7] <- max(par1[2],par2[2])
			cstart[10] <- min(par1[3],par2[3])
		}
	} else { cstart <- NULL }
	if (!is.null(llims) && !is.null(cstart)) { cstart[!is.na(llims)] <- pmax(cstart,llims)[!is.na(llims)] }
	if (!is.null(ulims) && !is.null(cstart)) { cstart[!is.na(ulims)] <- pmin(cstart,ulims)[!is.na(ulims)] }
	if (corrconc) { bfit <- findBestBRAID(act~cconc1+cconc2,wdata,defaults=defaults,startparv=cstart,llims=llims,ulims=ulims,itype=itype) }
	else { bfit <- findBestBRAID(act~conc1+conc2,wdata,defaults=defaults,startparv=cstart,llims=llims,ulims=ulims,itype=itype) }
	
	brdAnalysis <- list(concs=cbind(wdata$conc1,wdata$conc2),act=wdata$act,
							corrconc=corrconc,corrsigr=corrsigr,braidFit=bfit,hfit1=hfit1,hfit2=hfit2)
	return(brdAnalysis)
}
