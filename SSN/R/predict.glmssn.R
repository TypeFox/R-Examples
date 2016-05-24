predict.glmssn <- function(object, predpointsID, ...)
{
	if(length(object$estimates$Warnlog) > 0 &&
		length(grep("Algorithm diverging",object$estimates$Warnlog)) > 0)
		stop("No predictions for diverging algorithm")
	datao <- object$ssn.object@obspoints@SSNPoints[[1]]@point.data
	ocoord <- object$ssn.object@obspoints@SSNPoints[[1]]@point.coords
	theta <- object$estimates$theta
	ind <- object$sampinfo$ind.obs
	nobs <- length(ind)
	CorModels <- object$args$CorModels
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
	if(predpointsID == "_MissingObs_")
	{
		if(length(grep("tail",CorModels)) > 0 | useTailDownWeight == TRUE){
			if(length(grep("taildown",CorModels)) > 1)
				stop("Cannot have more than 1 tailup model")
			n.all <- length(datao[,1])
			dist.junc <- matrix(0, nrow = nobs, ncol = nobs)
			net.zero <-  matrix(0, nrow = nobs, ncol = nobs)
			nsofar <- 0
			nIDs <- sort(as.integer(as.character(unique(datao[,"netID"]))))
			for(i in nIDs){
				workspace.name <- paste("dist.net", i, ".RData", sep="")
				path <- file.path(object$ssn.object@path,
				"distance", "obs", workspace.name)
				file_handle <- file(path, "rb")
				distmat <- unserialize(file_handle)
				close(file_handle)

				oi <- order(as.numeric(rownames(distmat)))
				distmat <- distmat[oi, oi, drop = F]
				rnames <- c(rnames,rownames(distmat))
				ni <- length(distmat[1,])
				dist.junc[(nsofar + 1):(nsofar + ni),
					(nsofar + 1):(nsofar + ni)] <- distmat
				net.zero[(nsofar + 1):(nsofar + ni),
					(nsofar + 1):(nsofar + ni)] <- 1
				nsofar <- nsofar + ni
			}
			rownames(dist.junc) <- rnames
			# maximum distance to common junction between two sites
			a.mat <- pmax(dist.junc,t(dist.junc))
			# minimum distance to common junction between two sites
			b.mat <- pmin(dist.junc,t(dist.junc))
			# hydrological distance
			dist.hydro <- as.matrix(dist.junc + t(dist.junc))*net.zero
                        flow.con.mat <- 1 - (b.mat > 0)*1
			addfunccol <- object$args$addfunccol
			w.matrix <- sqrt(pmin(outer(datao[distord,addfunccol],
				rep(1, times = nobs)),
				t(outer(datao[distord,addfunccol],rep(1, times = nobs)))) /
				pmax(outer(datao[distord,addfunccol],rep(1, times = nobs)),
				t(outer(datao[distord,addfunccol],rep(1, times = nobs))))) *
				flow.con.mat*net.zero

			if(length(grep("tailup",CorModels)) > 0){
				if(length(grep("tailup",CorModels)) > 1)
					stop("Cannot have more than 1 tailup model")
				flow.con.mat <- 1 - (b.mat > 0)*1
			}
		}
		xcoord <- ocoord[distord,1,drop = F]
		ycoord <- ocoord[distord,2,drop = F]

		REind <- which(names(datao) %in% CorModels)
		if(length(REind)) {
			REs <- list()
			REnames <- sort(names(datao)[REind])
			## model matrix for a RE factor
      		for(ii in 1:length(REind)) REs[[REnames[ii]]] <-
				model.matrix(~datao[distord,REnames[ii]] - 1)
				rownames(REs[[REnames[ii]]]) <- datao[distord,"pid"]
			## corresponding block matrix
			for(ii in 1:length(REind)) REs[[ii]] <- REs[[ii]] %*% t(REs[[ii]])
		}

		Vpred <- makeCovMat(theta = theta, dist.hydro = dist.hydro,
			a.mat = a.mat, b.mat = b.mat, w.matrix = w.matrix,
			net.zero = net.zero, x.row = xcoord, y.row = ycoord,
			x.col = xcoord, y.col = ycoord,
			CorModels, useTailDownWeight = useTailDownWeight,
			use.nugget = object$args$use.nugget,
			use.anisotropy = FALSE, REs)
		for(i in 1:length(object$ssn.object@predpoints@SSNPoints))
			if(object$ssn.object@predpoints@ID[i] == "_MissingObs_"){
				datap <- object$ssn.object@predpoints@SSNPoints[[i]]@point.data
				distordp <- order(as.integer(as.character(datap[,"netID"])),
					datap[,"pid"])
				datap <- datap[distordp, , drop = F]
				pcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@
					point.coords
				pcoord <- pcoord[distordp, , drop = F]
				netpcoord <- object$ssn.object@predpoints@
					SSNPoints[[i]]@network.point.coords
				netpcoord <- netpcoord[distordp, , drop = F]
			}
		Vpred <- Vpred[!(datao[distord,"pid"] %in% datap[,"pid"]),
			datao[distord,"pid"] %in% datap[,"pid"], drop = F]
	}
	if(predpointsID != "_MissingObs_")
	{

		#get the predicted data data.frame and coordinates from glmssn object
		for(i in 1:length(object$ssn.object@predpoints@SSNPoints))
			if(object$ssn.object@predpoints@ID[i] == predpointsID){
				datap <- object$ssn.object@predpoints@SSNPoints[[i]]@point.data
				pcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@point.coords
				netpcoord <- object$ssn.object@predpoints@SSNPoints[[i]]@network.point.coords
			}

		npred <- length(datap[,1])
		distordp <- order(as.integer(as.character(datap[,"netID"])),
				datap[,"pid"])
		datap <- datap[distordp, , drop = F]
		pcoord <- pcoord[distordp, , drop = F]
		netpcoord <- netpcoord[distordp, , drop = F]
		rnames <- NULL
		cnames <- NULL
		if(length(grep("tail",CorModels)) > 0){
			if(length(grep("taildown",CorModels)) > 1)
				stop("Cannot have more than 1 tailup model")
			dist.junc.a <- matrix(0, nrow = nobs, ncol = npred)
			dist.junc.b <- matrix(0, nrow = npred, ncol = nobs)
			net.zero <-  matrix(0, nrow = nobs, ncol = npred)
			nSoFari <- 0
			nSoFarj <- 0
			nIDs <- sort(unique(c(as.integer(as.character(datao[,"netID"])),
				as.integer(as.character(datap[,"netID"])))))
			#loop through the networks to build covariance matrix
			for(i in nIDs) {
				ind.obs <- as.numeric(as.character(datao$netID)) ==
					as.numeric(i)
				site.no <- nrow(datao[ind.obs,])
				ind.pred <- as.numeric(as.character(datap$netID)) ==
					as.numeric(i)
				pred.no <- nrow(datap[ind.pred,])
				if(site.no*pred.no != 0) {
					workspace.name.a <- paste("dist.net", i,
						".a.RData", sep = "")
					workspace.name.b <- paste("dist.net", i,
						".b.RData", sep = "")
					path.a <- file.path(object$ssn.object@path,
						"distance", predpointsID, workspace.name.a)
					path.b <- file.path(object$ssn.object@path,
						"distance", predpointsID, workspace.name.b)
					# distance matrix a
					file_handle = file(path.a, open="rb")
					distmata <- unserialize(file_handle)
					close(file_handle)
					# distance matrix b
					file_handle = file(path.b, open="rb")
					distmatb <- unserialize(file_handle)
					close(file_handle)
    				# pid's for observed data in matrices above
					ordoi <- order(as.numeric(rownames(distmata)))
					ordpi <- order(as.numeric(rownames(distmatb)))
					ni <- length(ordoi)
					nj <- length(ordpi)
					distmata <- distmata[ordoi, ordpi, drop = F]
					distmatb <- distmatb[ordpi, ordoi, drop = F]
					if(any(rownames(distmata)!=colnames(distmatb)))
					stop("rownames of distmata do not match colnames of distmatb")
					if(any(colnames(distmata)!=rownames(distmatb)))
					stop("colnames of distmata do not match rownames of distmatb")
					dist.junc.a[(nSoFari + 1):(nSoFari + ni),
						(nSoFarj + 1):(nSoFarj + nj)] <- distmata
					dist.junc.b[(nSoFarj + 1):(nSoFarj + nj),
						(nSoFari + 1):(nSoFari + ni)] <- distmatb
					net.zero[(nSoFari + 1):(nSoFari + ni),
						(nSoFarj + 1):(nSoFarj + nj)] <- 1
				}
				else
				{
					ni <- site.no
					nj <- pred.no
				}
				nSoFari <- nSoFari + ni
				nSoFarj <- nSoFarj + nj
			}
			# rownames are for observed data, colnames are prediction data
			rnames <- datao[distord,"pid"]
			cnames <- datap[,"pid"]
			rownames(dist.junc.a) <- rnames
			colnames(dist.junc.a) <- cnames
			# creat A matrix (longest distance to junction of two points)
    		a.mat <- pmax(dist.junc.a,t(dist.junc.b))
    		# creat B matrix (shorted distance to junction of two points)
    		b.mat <- pmin(dist.junc.a,t(dist.junc.b))
    		# get hydrologic distance
    		dist.hydro <- as.matrix(dist.junc.a + t(dist.junc.b))

			# create indicator matrix of flow-connected
			if(length(grep("tailup",CorModels)) > 0){
				if(length(grep("tailup",CorModels)) > 1)
					stop("Cannot have more than 1 tailup model")
				flow.con.mat <- 1 - (b.mat > 0)*1
    			# weight matrix based on additive function
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
		# create Euclidean matrix among observed data
		xyobs <- ocoord
		x.samp <- ocoord[distord,1,drop = F]
		y.samp <- ocoord[distord,2,drop = F]
		xypred <- pcoord
		x.pred <- pcoord[, 1, drop = F]
		y.pred <- pcoord[, 2, drop = F]
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
			for(ii in 1:length(REnames)) REPs[[ii]] <-
				REOs[[ii]] %*% t(REPs[[ii]])
		}
		Vpred <- makeCovMat(theta = theta, dist.hydro = dist.hydro,
				a.mat = a.mat, b.mat = b.mat, w.matrix = w.matrix,
				net.zero = net.zero, x.row = x.samp, y.row = y.samp,
				x.col = x.pred, y.col = y.pred,
				CorModels = CorModels,
				useTailDownWeight = useTailDownWeight,
				use.nugget = FALSE,
				use.anisotropy = object$args$use.anisotropy, REs = REPs)
		Vpred <- Vpred[ind, , drop = F]


	}

	# get a list of response and covariate names
	response.col <- object$args$zcol
	mod.names <- as.character(attr(terms(object$args$formula,
		data = object$ssn.object@obspoints@SSNPoints[[1]]@point.data),"variables"))
	# get the number of names ( + 1, as the first is always "list")
	nc.tmp <- length(mod.names)
	# if there are any covariates ...
	ind.allcov <- rep(TRUE, times = length(datap[,1]))
	if(nc.tmp > 2) {
		# create a FALSE for a record with missing values of the covariates
		for(i in 3:nc.tmp) ind.allcov <- ind.allcov & !is.na(datap[,mod.names[i]])
	}
	# prediction sample size without missing covariates
	np.allcov <- sum(ind.allcov)
	# add response variable as -1 in datap
	datap[,response.col] <- NA
	# add prediction standard errors as NAs in datap
	datap[,paste(response.col,".predSE", sep = "")] <- NA
	# remove records that had any missing values for any of the covariates
	datap1 <- datap[ind.allcov,]
	# add response variable as -1 in datap
	datap1[,response.col] <- -1
	# add prediction standard errors as NAs in datap
	datap1[,paste(response.col,".predSE", sep = "")] <- NA

	formula <- object$args$formula
	# create design matrix for prediction data set
	mf <- model.frame(formula, data = datap1)
	mt <- attr(mf, "terms")
	Xpred <- model.matrix(mt, mf, contrasts)
	Xpred <- Xpred[,object$sampinfo$cutX1toX2, drop = F]

	# get the sum of partial sills
	sumparsil <- sum(theta[attr(theta,"type") == "parsill"])
  	Vi <- object$estimates$Vi
  	covb <- object$estimates$covb
  	Xobs <- object$sampinfo$X
  	z <- object$sampinfo$z
  	n <- object$sampinfo$obs.sample.size
  	p <- object$sampinfo$rankX

	parsilvec <- rep(sumparsil, times = length(Vpred[1,]))

	if(object$args$family == "poisson") {
		beta.hat <- object$estimates$betahat
		eta.hatp <- Xpred[ind.allcov,] %*% beta.hat
		eta.hato <- Xobs %*% beta.hat
		#diagonal elements of Delta~^{-1} of my manuscript
		Del.ip <- as.vector(1/exp(eta.hatp))
		Del.io <- as.vector(1/exp(eta.hato))
		#diagonal elements of A^(1/2) of my manuscript
		A.5p <- as.vector(sqrt(exp(eta.hatp)))
		A.5o <- as.vector(sqrt(exp(eta.hato)))
		Vpred <- t((Del.ip*A.5p)*t(Del.io*A.5o*Vpred))
		parsilvec <- sumparsil*(A.5p*Del.ip)^2
	}
	if(object$args$family == "binomial") {
		beta.hat <- object$estimates$betahat
		eta.hatp <- Xpred[ind.allcov,] %*% beta.hat
		eta.hato <- Xobs %*% beta.hat
		#diagonal elements of Delta~^{-1} of my manuscript
		Del.ip <- as.vector((1 + exp(eta.hatp))^2/exp(eta.hatp))
		Del.io <- as.vector((1 + exp(eta.hato))^2/exp(eta.hato))
		#diagonal elements of A^(1/2) of my manuscript
		A.5p <- as.vector(sqrt(exp(eta.hatp)/(1 +
			exp(eta.hatp))^2))
		A.5o <- as.vector(sqrt(exp(eta.hato)/(1 +
			exp(eta.hato))^2/object$sampinfo$trialsvec))
		Vpred <- t((Del.ip*A.5p)*t(Del.io*A.5o*Vpred))
		parsilvec <- sumparsil*(A.5p*Del.ip)^2
	}

  	M <- rbind(Vpred, t(Xpred), parsilvec)
  	XXSiXi <- Xobs %*% covb
  	XSi <- t(Xobs) %*% Vi
  	pred.out <- t(apply(M, 2, UK4Apply, covb = covb,
    	XXSiXi = XXSiXi, XSi = XSi, Vi = Vi, z = z, n = n, p = p))
  	datap1[,response.col] <- pred.out[,1]
  	datap1[,paste(response.col,".predSE", sep = "")] <- pred.out[,2]
  	datap[ind.allcov,] <- datap1
  	#put the predictions in the predicted data data.frame and return the SSN object
  	for(i in 1:length((object$ssn.object@predpoints@SSNPoints)))
    	if(object$ssn.object@predpoints@ID[i] == predpointsID){
      		object$ssn.object@predpoints@SSNPoints[[i]]@point.data <-
				datap[order(distordp),]
    	}
  	object$args$predpointsID <- predpointsID
  	class(object) <- "glmssn.predict"

  	object

}

