BlockPredict <- function(object, predpointsID)
{
	if(length(grep("binomial",tolower(object$args$family))) > 0)
		stop("No block prediction for binomial family")
	if(length(grep("poisson",tolower(object$args$family))) > 0)
		stop("No block prediction for Poisson family")
	datao <- object$ssn.object@obspoints@SSNPoints[[1]]@point.data
	ocoord <- object$ssn.object@obspoints@SSNPoints[[1]]@point.coords
	theta <- object$estimates$theta
	ind <- object$sampinfo$ind.obs
	nobs <- length(ind)
	cm <- object$args$CorModels
	useTailDownWeight <-object$args$useTailDownWeight
	REs <- object$sampinfo$REs
	distord <- order(as.integer(as.character(datao[,"netID"])),
				datao[,"pid"])

	a.mat <- NULL
	b.mat <- NULL
	net.zero <- NULL
	w.matrix <- NULL
	dist.hydro <- NULL
	xcoord <- NULL
	ycoord <- NULL
	xyobs <- NULL
	xypred <- NULL
	rnames <- NULL
	##get the predicted data data.frame and coordinates from glmssn object
	for(i in 1:length(object$ssn.object@predpoints@SSNPoints))
		if(object$ssn.object@predpoints@ID[i] == predpointsID){
			datap <- object$ssn.object@predpoints@SSNPoints[[i]]@point.data
			pcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@point.coords
			netpcoord <- object$ssn.object@predpoints@
				SSNPoints[[i]]@network.point.coords
		}

	npred <- length(datap[,1])
	distordp <- order(as.integer(as.character(datap[,"netID"])),
 		datap[,"pid"])
	datap <- datap[distordp,]
	pcoord <- pcoord[distordp,]
	netpcoord <- netpcoord[distordp,]
	rnames <- NULL
	cnames <- NULL
	if(any(c(grep("tailup",cm),grep("taildown",cm)))) {
		nets <- levels(object$ssn.object@data$netID)
		net.count <- length(nets)
		dist.junc.a <- matrix(0, nrow = nobs, ncol = npred)
		dist.junc.b <- matrix(0, nrow = npred, ncol = nobs)
		net.zero <-  matrix(0, nrow = nobs, ncol = npred)
		nSoFari <- 0
		nSoFarj <- 0
		nIDs <- sort(unique(c(as.integer(as.character(datao[,"netID"])),
		as.integer(as.character(datap[,"netID"])))))
 		##loop through the networks to build covariance matrix
		for(i in 1:length(nIDs)) {
			ind.obs <- object$ssn.object@obspoints@
				SSNPoints[[1]]@point.data$netID == as.numeric(nIDs[i])
			site.no <- nrow(object$ssn.object@obspoints@
				SSNPoints[[1]]@point.data[ind.obs,])
			ind.pred <- datap$netID == as.numeric(i)
			pred.no <- nrow(datap[ind.pred,])
			if(site.no*pred.no != 0) {
				workspace.name.a <- paste("dist.net", nIDs[i], 
					".a.RData", sep = "")
				workspace.name.b <- paste("dist.net", nIDs[i], 
					".b.RData", sep = "")
				path.a <- file.path(object$ssn.object@path, 
					"distance", predpointsID, workspace.name.a)
				path.b <- file.path(object$ssn.object@path, 
					"distance", predpointsID, workspace.name.b)
				## distance matrix a
				file_handle = file(path.a, open="rb")
				distmata <- unserialize(file_handle)
				close(file_handle)
				## distance matrix b
				file_handle = file(path.b, open="rb")
				distmatb <- unserialize(file_handle)
				close(file_handle)
				## pid's for observed data in matrices above
				ordoi <- order(as.numeric(rownames(distmata)))
				ordpi <- order(as.numeric(rownames(distmatb)))
				ni <- length(ordoi)
				nj <- length(ordpi)
				distmata <- distmata[ordoi, ordpi]
				distmatb <- distmatb[ordpi, ordoi]
				if(any(rownames(distmata)!=colnames(distmatb)))
					stop("rownames of distmata do not match colnames of distmatb")
				if(any(colnames(distmata)!=rownames(distmatb)))
					stop("colnames of distmata do not match rownames of distmatb")
				# rownames are are adopted from A
				dist.junc.a[(nSoFari + 1):(nSoFari + ni),
					(nSoFarj + 1):(nSoFarj + nj)] <- distmata
				dist.junc.b[(nSoFarj + 1):(nSoFarj + nj),
					(nSoFari + 1):(nSoFari + ni)] <- distmatb
 				net.zero[(nSoFari + 1):(nSoFari + ni),
					(nSoFarj + 1):(nSoFarj + nj)] <- 1
			} else {
				ni <- sum(as.integer(as.character(datao[,"netID"])) == nIDs[i])
				nj <- sum(as.integer(as.character(datap[,"netID"])) == nIDs[i])
			}
			nSoFari <- nSoFari + ni
			nSoFarj <- nSoFarj + nj
		}
		rownames(dist.junc.a) <- rownames(datao)[distord]
		colnames(dist.junc.a) <- rownames(datap)[distordp]
		## creat A matrix (longest distance to junction of two points)
		a.mat <- pmax(dist.junc.a,t(dist.junc.b))
		## creat B matrix (shorted distance to junction of two points)
		b.mat <- pmin(dist.junc.a,t(dist.junc.b))
		## get hydrologic distance
		dist.hydro <- as.matrix(dist.junc.a + t(dist.junc.b))

		## create indicator matrix of flow-connected
		if(any(grep("tailup",cm))) {
			flow.con.mat <- 1 - (b.mat > 0)*1
			## weight matrix based on additive function
			addfunccol <- object$args$addfunccol
			w.matrix <- sqrt(pmin(outer(datao[distord,addfunccol],
				rep(1, times = npred)),
				t(outer(datap[,addfunccol],rep(1, times = nobs)))) /
				pmax(outer(datao[distord,addfunccol],rep(1, times = npred)),
				t(outer(datap[,addfunccol],rep(1, times = nobs))))) *
				flow.con.mat*net.zero
			if(any(rownames(w.matrix)!=rownames(a.mat)))
				stop("rownames of w.matrix do not match rownames of a.mat")
			if(any(colnames(w.matrix)!=colnames(a.mat)))
				stop("colnames of w.matrix do not match colnames of a.mat")
		}
	}
	# necessary Euclidean distance data
	xyobs <- ocoord
	x.samp <- ocoord[distord,1,drop = F]
	y.samp <- ocoord[distord,2,drop = F]
	xypred <- pcoord
	x.pred <- pcoord[distordp,1,drop = F]
	y.pred <- pcoord[distordp,2,drop = F]
	if(any(rownames(x.samp)!=rownames(a.mat)))
	stop("rownames of x.samp do not match rownames of a.mat")
	if(any(rownames(x.pred)!=colnames(a.mat)))
	stop("rownames of x.pred do not match colnames of a.mat")
	REPs <- NULL
	if(!is.null(REs)) {
		REnames <- names(REs)
		for(ii in 1:length(REnames)) if(any(is.na(datap[,REnames[ii]])))
			stop("Cannot having missing values when creating random effects")
		REOs <- list()
		REPs <- list()
		ObsSimDF <- datao
		PredSimDF <- datap
		## model matrix for a RE factor
		for(ii in 1:length(REnames)){ 
			#we'll add "o" to observed levels and "p" to prediction
			# levels so create all possible levels
			plevels <- unique(c(levels(PredSimDF[,REnames[[ii]]]),
			paste("o",levels(ObsSimDF[,REnames[[ii]]]),sep = ""),
			paste("p",levels(PredSimDF[,REnames[[ii]]]),sep = "")))
			# sites with prediction levels same as observation levels
			pino <- PredSimDF[,REnames[[ii]]] %in% ObsSimDF[,REnames[[ii]]]
			#add "o" to observed levels
			ObsSimDF[,REnames[[ii]]] <- paste("o", 
				ObsSimDF[,REnames[[ii]]], sep = "")
			ObsSimDF[,REnames[[ii]]] <- as.factor(as.character(
				ObsSimDF[,REnames[[ii]]]))
			#add all possible levels to prediction data frame
			levels(PredSimDF[,REnames[[ii]]]) <- plevels
			# add "o" to prediction sites with observation levels
			if(any(pino)) PredSimDF[pino,REnames[[ii]]] <- paste("o", 
				PredSimDF[pino,REnames[[ii]]], sep = "")
			# add "p" to all predicition sites without observation levels
			if(any(!pino)) PredSimDF[!pino,REnames[[ii]]] <- paste("p", 
				PredSimDF[!pino,REnames[[ii]]], sep = "")
			PredSimDF[,REnames[[ii]]] <- as.factor(as.character(
				PredSimDF[,REnames[[ii]]]))
			# now get down to just levels with "o" & "p" added
			blevels <- unique(c(levels(ObsSimDF[,REnames[[ii]]]),
				levels(PredSimDF[,REnames[[ii]]])))
			ObsSimDF[,REnames[[ii]]] <- factor(ObsSimDF[,REnames[[ii]]],
				levels = blevels, ordered = FALSE)
			PredSimDF[,REnames[[ii]]] <- factor(PredSimDF[,REnames[[ii]]],
				levels = blevels, ordered = FALSE)
			# now ordering of factors in Z matrices should be compatible
			# with obs x obs Z matrices
			REOs[[ii]] <- model.matrix(~ObsSimDF[distord,
				REnames[[ii]]] - 1)
			REPs[[ii]] <- model.matrix(~PredSimDF[,
				REnames[[ii]]] - 1)
			rownames(REOs[[ii]]) <- datao[distord,"pid"]
			rownames(REPs[[ii]]) <- datap[,"pid"]
			if(any(rownames(REOs[[ii]])!=rownames(a.mat)))
				stop("rownames RE for obs do not match rownames of a.mat")
			if(any(rownames(REPs[[ii]])!=colnames(a.mat)))
				stop("rownames RE for preds do not match colnames of a.mat")
		}
		## corresponding block matrix 
		for(ii in 1:length(REnames)) REPs[[ii]] <- REOs[[ii]] %*% t(REPs[[ii]])
	}
	Vpred <- makeCovMat(theta = theta, dist.hydro = dist.hydro,
		a.mat = a.mat, b.mat = b.mat, w.matrix = w.matrix,
		net.zero = net.zero, x.row = x.samp, y.row = y.samp,
		x.col = x.pred, y.col = y.pred, useTailDownWeight = useTailDownWeight,
		CorModels = cm, use.nugget = FALSE,
		use.anisotropy = object$args$use.anisotropy, REs = REPs)
	Vpred <- Vpred[ind,]

		#---------------------------------------------------------------------
		#               PREDICTION BY PREDICTION COVARIANCE MATRIX
		#---------------------------------------------------------------------
		a.mat <- NULL
		b.mat <- NULL
		net.zero <- NULL
		w.matrix <- NULL
		dist.hydro <- NULL
		xcoordp <- NULL
		ycoordp <- NULL
		dist.Euclid <- NULL
		if(length(grep("tail",cm)) > 0){
			nets <- levels(datap$netID)
			net.count <- length(nets)		  
			dist.junc <- matrix(0, nrow = npred, ncol = npred)
			net.zero <-  matrix(0, nrow = npred, ncol = npred)
			nsofar <- 0
			nIDs <- as.integer(as.character(sort(unique(datap[,"netID"]))))
			for(net.num in nIDs){
				ind.pred <- datap$netID == as.numeric(net.num)
				pred.no <- nrow(datap[ind.pred,])
				if(pred.no != 0) {
					workspace.name <- paste("dist.net", net.num, 
						".RData", sep = "")
					path <- file.path(object$ssn.object@path, 
                    	"distance", predpointsID, workspace.name)
					if(!file.exists(path)) {
						stop("Unable to locate required distance matrix")
					}
					file_handle <- file(path, open="rb")
					distmat <- unserialize(file_handle)
					close(file_handle)
					ni <- length(distmat[1,])
					ordpi <- order(as.numeric(rownames(distmat)))
					dist.junc[(nsofar + 1):(nsofar + ni),(nsofar + 1):
						(nsofar + ni)] <- distmat[ordpi, ordpi]
					net.zero[(nsofar + 1):(nsofar + ni),(nsofar + 1):
						(nsofar + ni)] <- 1
					nsofar <- nsofar + ni
				}
			}
			# maximum distance to common junction between two sites
			a.mat <- pmax(dist.junc,t(dist.junc))
			# minimum distance to common junction between two sites
			b.mat <- pmin(dist.junc,t(dist.junc))
			# hydrological distance
			dist.hydro <- as.matrix(dist.junc + t(dist.junc))*net.zero
			if(length(grep("tailup",cm)) > 0){
				flow.con.mat <- 1 - (b.mat > 0)*1			
				w.matrix <- sqrt(pmin(outer(datap[,addfunccol],
					rep(1, times = npred)),
					t(outer(datap[,addfunccol],rep(1, times = npred)))) /
					pmax(outer(datap[,addfunccol],rep(1, times = npred)),
					t(outer(datap[,addfunccol],rep(1, times = npred))))) *
					  flow.con.mat*net.zero
			}		  
		}
		xcoordp <- pcoord[, 1, drop = F]
		ycoordp <- pcoord[, 2, drop = F]
		REPPs <- NULL	  
		if(!is.null(REs)) {
			REPPs <- list()
			for(ii in 1:length(REnames)){ 
				## model matrix for a RE factor
				REPPs[[ii]] <- model.matrix(~PredSimDF[,REnames[ii]] - 1)
				## corresponding block matrix 
				REPPs[[ii]] <- REPPs[[ii]] %*% t(REPPs[[ii]])
			}
		}
		Vpredpred <- makeCovMat(theta = theta, dist.hydro, 
			a.mat, b.mat, w.matrix, net.zero, xcoordp, ycoordp, 
			xcoordp, ycoordp, cm, use.nugget = object$args$use.nugget, 
			use.anisotropy = object$args$use.anisotropy,
			REPPs, useTailDownWeight = useTailDownWeight)


    ## get a list of response and covariate names
	response.col <- object$args$zcol
	mod.names <- as.character(attr(terms(object$args$formula,
		data = object$ssn.object@obspoints@SSNPoints[[1]]@point.data),"variables"))
        ## get the number of names ( + 1, as the first is always "list")
	nc.tmp <- length(mod.names)
        ## if there are any covariates ...
	ind.allcov <- rep(TRUE, times = length(datap[,1]))
	if(nc.tmp > 2) {
            ## create a FALSE for a record with missing values of the covariates
            for(i in 3:nc.tmp) ind.allcov <- ind.allcov & !is.na(datap[,mod.names[i]])
	}
        ## prediction sample size without missing covariates
	np.allcov <- sum(ind.allcov)
        ## add response variable as -1 in datap
	datap[,response.col] <- NA
	## add prediction standard errors as NAs in datap
	datap[,paste(response.col,".predSE", sep = "")] <- NA
        ## remove records that had any missing values for any of the covariates
	datap1 <- datap[ind.allcov,]
        ## add response variable as -1 in datap
	datap1[,response.col] <- -1
	## add prediction standard errors as NAs in datap
	datap1[,paste(response.col,".predSE", sep = "")] <- NA

	formula <- object$args$formula
	## create design matrix for prediction data set
	mf <- model.frame(formula, data = datap1)
	mt <- attr(mf, "terms")
	Xpred <- model.matrix(mt, mf, contrasts)
	Xpred <- as.matrix(Xpred[,object$sampinfo$cutX1toX2])

	Vi <- object$estimates$Vi
	covb <- object$estimates$covb
	Xobs <- object$sampinfo$X
	z <- object$sampinfo$z
	n <- object$sampinfo$obs.sample.size
	p <- object$sampinfo$rankX

	Vp <- Vpred[,ind.allcov]
    Vpp <- Vpredpred[ind.allcov, ind.allcov]

	XXSiXi <- Xobs %*% covb
	XSi <- t(Xobs) %*% Vi
	XB <- apply(Xpred,2,mean)
	cB <- apply(Vp,1,mean)
	CBB <- mean(Vpp)
	r1 <- XB - XSi %*% cB
	m <- covb %*% r1
	tlam <- t(cB + XXSiXi %*% r1) %*% Vi
	outpt <- data.frame(
		BlockPredEst = tlam %*% z,
		BlockPredSE = sqrt(CBB - tlam %*% cB + t(m) %*% XB)
	)
	outpt
}


