glmssn1 <- function(formula, ssn.object,
	family = "Gaussian",
	CorModels = c("Exponential.tailup", "Exponential.taildown",
		"Exponential.Euclid"),
	use.nugget = TRUE,
	use.anisotropy = FALSE,
	addfunccol = NULL,
	trialscol = NULL,
	EstMeth = "REML",
	useTailDownWeight = FALSE,
	trans.power = NULL,
	trans.shift = 0,
	control = list(max.range.factor = 4,
		trunc.pseudo = NULL,
		maxiter.pseudo = 20,
		beta.converge = 1e-5)
)
{

        Warnlog <- NULL
	data <- ssn.object@obspoints@SSNPoints[[1]]@point.data
	data <- cbind(data, ssn.object@obspoints@SSNPoints[[1]]@point.coords)
	xcol <- "coords.x1"
	ycol <- "coords.x2"
	family <- tolower(family)
  if(!"max.range.factor" %in% names(control)) control$max.range.factor <- 4
  if(!"maxiter.pseudo" %in% names(control)) control$maxiter.pseudo <- 20

	nIDs <- sort(as.integer(as.character(unique(data[,"netID"]))))
	dist.junc <- matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
	net.zero <-  matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
	nsofar <- 0
	distord <- order(data[,"netID"],data[,"pid"])
	names(distord) <- rownames(data)[distord]
	for(i in nIDs) {
		workspace.name <- paste("dist.net", i, ".RData", sep = "")
		path <- file.path(ssn.object@path, "distance", "obs",
			workspace.name)
		if(!file.exists(path)) {
			stop("Unable to locate required distance matrix")
		}
		file_handle <- file(path, open="rb")
		distmat <- unserialize(file_handle)
		ordpi <- order(as.numeric(rownames(distmat)))
		close(file_handle)
		ni <- length(distmat[1,])
		dist.junc[(nsofar + 1):(nsofar + ni),(nsofar + 1):
			(nsofar + ni)] <- distmat[ordpi, ordpi, drop = F]
		net.zero[(nsofar + 1):(nsofar + ni),(nsofar + 1):(nsofar + ni)] <- 1
		nsofar <- nsofar + ni
	}

	## Check all arguments have correct form
	Err <- arg.error.check.multi(CorModels = CorModels,
		use.nugget = use.nugget,
		use.anisotropy = use.anisotropy,
		addfunccol = addfunccol,
		family = family,
		EstMeth = EstMeth,
		ssn = ssn.object)

	if(Err$Err == 1) return(print(Err$message))

	mf <- match.call(expand.dots = FALSE)
	dataXY.out <- dataXY(formula, data,
		family = family, trialscol = trialscol,
		trans.power = trans.power,
		trans.shift = trans.shift,
		CorModels = CorModels,
		distord = distord)

	REs <- dataXY.out$REs
        REmodelmatrices <- dataXY.out$REmodelmatrices
	n.all <- dataXY.out$sampsizes$n.all
	a.mat <- NULL
	b.mat <- NULL
	a.mat.data <- NULL
	b.mat.data <- NULL
	dist.hydro <- NULL
	dist.hydro.data <- NULL
	w.matrix.data <- NULL
	ind <- dataXY.out$indvecs$ind.allxy

	## create any necessary matrices from distance and flow matrices
	if(!is.null(dist.junc) ) {
		## maximum distance to common junction between two sites
		a.mat <- pmax(dist.junc,t(dist.junc))
		a.mat.data <- a.mat[ind,ind]
		## minimum distance to common junction between two sites
		b.mat <- pmin(dist.junc,t(dist.junc))
		b.mat.data <- b.mat[ind,ind]
		## hydrological distance
		dist.hydro <- as.matrix(dist.junc + t(dist.junc))
		## subset stream distance to observed locations only
		dist.hydro.data <- dist.hydro[ind, ind]
	}
	if(length(grep("tailup",CorModels)) | useTailDownWeight == TRUE) {
		if(missing(addfunccol) || is.null(addfunccol) ||
			length(addfunccol) == 0 || !(addfunccol %in% colnames(data)))
	    	stop("The specified value for addfunccol was invalid")
		flow.con.mat <- 1 - (b.mat > 0)*1

		w.matrix <- sqrt(pmin(outer(data[distord,addfunccol],
			rep(1, times = n.all)),
			t(outer(data[distord,addfunccol],rep(1, times = n.all)))) /
			pmax(outer(data[distord,addfunccol],rep(1, times = n.all)),
			t(outer(data[distord,addfunccol],rep(1, times = n.all)))))*
			flow.con.mat*net.zero
		w.matrix.data <- w.matrix[ind, ind]
	}
	net.zero.data <- net.zero[ind,ind]

	xcoord <- data[distord,xcol]
	ycoord <- data[distord,ycol]
	xcoord.data <- xcoord[ind]
	ycoord.data <- ycoord[ind]

	z <- dataXY.out$respvecs$z
	X2 <- dataXY.out$Xmats$X2
	n.allxy <- dataXY.out$sampsizes$n.allxy

	trialsvec <- NULL
	## Initial parameter estimates, fixed effects, pseudo-data
	if(family == "binomial") {
		if(is.null(trialscol)) trialsvec <- rep(1, times = length(z))
		if(!is.null(trialscol)) trialsvec <-
			dataXY.out$datasets$data2[,trialscol]
		beta.hat <- glm(formula, data,
			family = "binomial")$coefficients
		beta.current <- beta.hat
		eta.hat <- X2 %*% beta.hat
		##diagonal elements of Delta~^{-1} of my manuscript
		Del.i <- as.vector((1 + exp(eta.hat))^2/exp(eta.hat))
		##diagonal elements of A^(1/2) of my manuscript
		A.5 <- as.vector(sqrt(exp(eta.hat)/(1 +
			exp(eta.hat))^2/trialsvec))
		##Binomial pseudo data
		zt <- Del.i*(z - exp(eta.hat)/(1 + exp(eta.hat))) + eta.hat
		if(!is.null(control$trunc.pseudo)) {
			zt[zt > control$trunc.pseudo] <- control$trunc.pseudo
			zt[zt < -control$trunc.pseudo] <- -control$trunc.pseudo
			if(max(abs(zt)) == control$trunc.pseudo) Warnlog <- c(Warnlog,
				paste("Psuedo-data truncated to +/-",control$trunc.pseudo))
		}
		## stop if zt gets large; causes numerical instabililty when exponentiated
		if(max(abs(zt)) > 80 & is.null(control$trunc.pseudo))
			return(list("Bad initial pseudo-data",zt))
	}
	if(family == "poisson") {
		beta.hat <- glm(formula,
			data, family = "poisson")$coefficients
		beta.current <- beta.hat
		eta.hat <- X2 %*% beta.hat
		##diagonal elements of Delta~^{-1} of my manuscript
		Del.i <- as.vector(1/exp(eta.hat))
		##diagonal elements of A^(1/2) of my manuscript
		A.5 <- as.vector(sqrt(exp(eta.hat)))
		##Poisson pseudo data
		zt <- Del.i*(z - exp(eta.hat)) + eta.hat
		if(!is.null(control$trunc.pseudo)) {
			zt[zt > control$trunc.pseudo] <- control$trunc.pseudo
			zt[zt < -control$trunc.pseudo] <- -control$trunc.pseudo
			if(max(abs(zt)) == control$trunc.pseudo) Warnlog <- c(Warnlog,
				paste("Psuedo-data truncated to +/-",control$trunc.pseudo))
		}
		## stop if zt gets large; causes numerical instabililty when exponentiated
		if(max(abs(zt)) > 80 & is.null(control$trunc.pseudo))
			return(list("Bad initial pseudo-data",zt))
	}
	if(family == "gaussian"){
		A.5 <-  NULL
		Del.i <-  NULL
		zt <- z
	}

	##We're going to pass this environment in to the m2LL.stream
	##funtion and it will track all evaluations of the -loglik surface
	loglik.environment <- environment()
	assign("RESULT",NULL,loglik.environment)

	#set maximum range parameters
	maxrang <- NULL
	mrf <- control$max.range.factor
	if(length(grep("tailup",CorModels)) > 0)
		maxrang <- c(maxrang, NA, mrf*max(dist.hydro))
	if(length(grep("taildown",CorModels)) > 0)
		maxrang <- c(maxrang, NA, mrf*max(dist.hydro))
	if(length(grep("Euclid",CorModels)) > 0)
		maxrang <- c(maxrang, NA, mrf*max(distGeo(xcoord.data,
			ycoord.data, xcoord.data, ycoord.data, 1)))
	if(length(REs)) maxrang <- c(maxrang, rep(NA, times = length(REs)))
	if(use.nugget == TRUE) maxrang <- c(maxrang, NA)

	# ----------- START LOOPING HERE ---------------------------

	# create an indicator to stop looping
	stoploop <- 0
	# keep track of the number of iterations
	iter <- 0
	# keep track of number of inner iterations for beta
	inner.iter2 <- NULL
	# begin looping
	while(stoploop == 0) {
		## initial parameter estimates
		if(iter == 0) {
			theta <- theta.ini(z = zt, X = X2,
				CorModels= CorModels,
					use.nugget = use.nugget, use.anisotropy = use.anisotropy,
					dist.hydro.data = dist.hydro.data, x.dat = xcoord.data,
					y.dat = ycoord.data, REs = REs)
					attributes(theta) ## scale, type, terms
					TH.scale <- attributes(theta)$scale
					TH.type <- attributes(theta)$type
					TH.terms <- attributes(theta)$terms
		}

		if(length(theta > 0)) {
		## optimizing for covariance parameter estimates using ML or REML
			if(length(theta) ==1) {
				lowerb <- log(.1*exp(theta))
				upperb <- log(10*exp(theta))
				parmest1.out <- optimize(m2LL.stream, interval =
					c(lowerb,upperb), m2LLdata = zt, X = X2,
					dist.hydro = dist.hydro.data, weight = w.matrix.data,
					net.zero = net.zero.data,
					a.mat = a.mat.data, b.mat = b.mat.data,
					x.dat = xcoord.data, y.dat = ycoord.data,
					Del.i = Del.i, A.5 = A.5,
					CorModels = CorModels, scale = TH.scale, 
					useTailDownWeight = useTailDownWeight,
					use.nugget = use.nugget, use.anisotropy = use.anisotropy,
					EstMeth = EstMeth, loglik.environment=loglik.environment,
					REs = REs, maxrang = maxrang)
				lowerb <- log(.5*exp(parmest1.out$minimum))
				upperb <- log(2*exp(parmest1.out$minimum))
				parmest2.out <- optimize(m2LL.stream, interval =
					c(lowerb,upperb), m2LLdata = zt, X = X2,
				dist.hydro = dist.hydro.data, weight = w.matrix.data,
					net.zero = net.zero.data,
					a.mat = a.mat.data, b.mat = b.mat.data,
					x.dat = xcoord.data, y.dat = ycoord.data,
					Del.i = Del.i, A.5 = A.5,
					CorModels = CorModels, scale = TH.scale, 
					useTailDownWeight = useTailDownWeight,
					use.nugget = use.nugget, use.anisotropy = use.anisotropy,
					EstMeth = EstMeth, loglik.environment=loglik.environment,
					REs = REs, maxrang = maxrang)
				parmest.out <- parmest2.out
				theta <- parmest2.out$minimum
				m2LL <- parmest2.out$objective
			}

			if(length(theta) > 1) {
				if(iter == 0) {
				## try Nelder-Mead
					parmest1.out <- optim(theta, m2LL.stream, m2LLdata = zt,
						X = X2, dist.hydro = dist.hydro.data, weight =
						w.matrix.data, net.zero = net.zero.data,
						a.mat = a.mat.data, b.mat = b.mat.data,
						x.dat = xcoord.data, y.dat = ycoord.data,
						Del.i = Del.i, A.5 = A.5,
						CorModels = CorModels, useTailDownWeight = useTailDownWeight,
						use.nugget = use.nugget, use.anisotropy = use.anisotropy,
						EstMeth = EstMeth,
						method = "BFGS",loglik.environment=
						loglik.environment, REs = REs, scale = TH.scale,
						maxrang = maxrang)
				}
				## try BFGS
				parmest2.out <- optim(theta, m2LL.stream, m2LLdata = zt, X = X2,
					dist.hydro = dist.hydro.data, weight = w.matrix.data,
					net.zero = net.zero.data,
					a.mat = a.mat.data, b.mat = b.mat.data,
					x.dat = xcoord.data, y.dat = ycoord.data,
					Del.i = Del.i, A.5 = A.5,
					CorModels = CorModels, useTailDownWeight = useTailDownWeight,
					use.nugget = use.nugget, use.anisotropy = use.anisotropy,
					EstMeth = EstMeth, method = "Nelder-Mead",loglik.environment=
					loglik.environment,hessian=TRUE, scale = TH.scale,
					REs = REs, maxrang = maxrang)
				if(iter == 0 & parmest1.out$value <
					parmest2.out$value) parmest.out <- parmest1.out
				else parmest.out <- parmest2.out
				theta <- parmest.out$par
				m2LL <- parmest.out$value
			}

				## go back to original scale for covariance parameters
			parmest <- untrans.theta(theta = theta, scale = TH.scale)
		}

		if(is.null(theta)) parmest <- 1
		if(length(CorModels) == 0) V <- diag(rep(parmest, times = n.allxy))
		else V <- makeCovMat(parmest, dist.hydro = dist.hydro.data,
			w.matrix = w.matrix.data, net.zero = net.zero.data,
			a.mat = a.mat.data, b.mat = b.mat.data, 
			x.row = xcoord.data, y.row = ycoord.data,
			x.col = xcoord.data, y.col = ycoord.data,
			CorModels = CorModels, useTailDownWeight = useTailDownWeight,
			use.nugget = use.nugget, use.anisotropy = use.anisotropy, REs)

		if(!is.null(Del.i)) {
			## V <- diag(Del.i) %*% diag(A.5) %*% V %*% diag(A.5) %*% diag(Del.i)
			V <- Del.i*A.5*t((Del.i*A.5) * V)
		}
		##	if(min(svd(V)$d) < 0 ) {browser()}
		qrV <- try(qr(V), silent = T)
		##	if(class(qrV) != "qr") {browser()}
		ViX <- try(solve(qrV,X2), silent = T)
		if(class(ViX) == "try-error")
			return("Algorith diverged -- estimation stopped")
		covbi <- t(X2) %*% ViX
		covb <- solve(covbi)
		if(family != "gaussian") beta.old <- beta.hat
		beta.hat <- covb %*% t(ViX) %*% zt
		if(family == "gaussian") beta.old <- beta.hat

		if(family == "binomial") {
			beta.current <- beta.hat
			eta.hat <- X2 %*% beta.hat
			##diagonal elements of Delta~^{-1} of my manuscript
			Del.i <- as.vector((1 + exp(eta.hat))^2/exp(eta.hat))
			##diagonal elements of A^(1/2) of my manuscript
			A.5 <- as.vector(sqrt(exp(eta.hat)/(1 + exp(eta.hat))^2/trialsvec))
			##Binomial pseudo data
			zt <- Del.i*(z - exp(eta.hat)/(1 + exp(eta.hat))) + eta.hat
			if(!is.null(control$trunc.pseudo)) {
				zt[zt > control$trunc.pseudo] <- control$trunc.pseudo
				zt[zt < -control$trunc.pseudo] <- -control$trunc.pseudo
				if(max(abs(zt)) == control$trunc.pseudo) Warnlog <- c(Warnlog,
					paste("Psuedo-data truncated to +/-",control$trunc.pseudo))
			}
			##			browser()
			## stop if zt gets large; causes numerical instabililty when exponentiated
			if(max(abs(zt)) > 80 & is.null(control$trunc.pseudo)) {
				Warnlog <- c(Warnlog,
					"Algorithm diverging, thus terminated -- results questionable")
				stoploop <- 1
			}
		}
		if(family == "poisson") {
			beta.current <- beta.hat
			eta.hat <- X2 %*% beta.hat
			##diagonal elements of Delta~^{-1} of my manuscript
			Del.i <- as.vector(1/exp(eta.hat))
			##diagonal elements of A^(1/2) of my manuscript
			A.5 <- as.vector(sqrt(exp(eta.hat)))
			##Poisson pseudo data
			zt <- Del.i*(z - exp(eta.hat)) + eta.hat
			if(!is.null(control$trunc.pseudo)) {
				zt[zt > control$trunc.pseudo] <- control$trunc.pseudo
				zt[zt < -control$trunc.pseudo] <- -control$trunc.pseudo
				if(max(abs(zt)) == control$trunc.pseudo) Warnlog <- c(Warnlog,
					paste("Psuedo-data truncated to +/-",control$trunc.pseudo))
			}
			##			browser()
			## stop if zt gets large; causes numerical instabililty when exponentiated
			if(max(abs(zt)) > 80 & is.null(control$trunc.pseudo)) {
				Warnlog <- c(Warnlog,
					"Algorithm diverging, thus terminated -- results questionable")
				stoploop <- 1
			}
		}
		if(family == "gaussian"){
			A.5 <-  NULL
			Del.i <-  NULL
			zt <- z
		}
		##convergence criteria on the fixed effect parameters
		non0ind <- beta.old != 0
		if(all(abs((beta.hat[non0ind] - beta.old[non0ind])/
			beta.old[non0ind]) < control$beta.converge)) stoploop <- 1
		if (iter > control$maxiter.pseudo) {
			Warnlog <- c(Warnlog,
			paste("More than", control$maxiter.pseudo,
				"iterations, algorithm terminated -- results questionable"))
			stoploop <- 1
		}

		iter <- iter + 1

	}

	# ----------- DONE LOOPING HERE ---------------------------

	# inverse covariance matrix between observed locations
	Vi <- solve(qrV)

	ViX <- solve(qrV,X2)
	covbi <- t(X2) %*% ViX
	covb <- solve(covbi)
	bhat.se <- sqrt(diag(covb))
	b.hat <- covb %*% t(X2) %*% Vi %*% zt

    bhat.se <- sqrt(diag(covb))
    p <-  dataXY.out$Xmats$p
	if(use.nugget == TRUE)  nugget <- parmest[length(parmest)]

	is.na(data[,dataXY.out$respvecs$response.col])
	nobs <- length(ssn.object@obspoints@SSNPoints[[1]]@point.data[,1])
	# if any missing data, create prediction set in the SSN object
	if(dataXY.out$sampsizes$n.allcov > n.allxy) {
		TempSSNPoints <- ssn.object@obspoints@SSNPoints
		TempPredPoints <- ssn.object@predpoints
		listlen <- length(TempPredPoints@SSNPoints)
		TempPredPoints@SSNPoints[[listlen+1]] <- TempSSNPoints[[1]]
	  	TempPredPoints@SSNPoints[[listlen+1]]@network.point.coords <-
			TempPredPoints@SSNPoints[[listlen+1]]@network.point.coords[
			is.na(data[,dataXY.out$respvecs$response.col]), ,drop = F]
		TempPredPoints@SSNPoints[[listlen+1]]@point.coords <-
			junk <- TempPredPoints@SSNPoints[[listlen+1]]@point.coords[
			is.na(data[,dataXY.out$respvecs$response.col]), ,drop = F]
		TempPredPoints@SSNPoints[[listlen+1]]@point.data <-
			TempPredPoints@SSNPoints[[listlen+1]]@point.data[
			is.na(data[,dataXY.out$respvecs$response.col]), ,drop = F]
		TempPredPoints@ID <- c(TempPredPoints@ID, "_MissingObs_")
		ssn.object@predpoints <- TempPredPoints
	}

	if(is.null(theta)) {
		nugget <- 1
		parmest <- 1
		attr(parmest,"scale") <- "natural"
		attr(parmest, "type") <- "parsill"
		attr(parmest, "terms") <- "nugget(held at 1)"
		m2LL <- NULL
		parmest.out <- NULL
	} else {
    	attr(parmest,"scale") <- TH.scale
    	attr(parmest,"type") <- TH.type
    	attr(parmest,"terms") <- TH.terms
	}
	if(use.nugget == FALSE) nugget <- 0
	outpt <- list(
		args = list(
		formula = formula,
		zcol = dataXY.out$respvecs$response.col,
		family = family,
		CorModels = CorModels,
		useTailDownWeight = useTailDownWeight,
		use.nugget = use.nugget,
		use.anisotropy = use.anisotropy,
		addfunccol = addfunccol,
		trialscol = trialscol,
	    EstMeth = EstMeth,
	    trans.power = trans.power,
		trans.shift = trans.shift,
		algorithm = "orig"
	),
	ssn.object = ssn.object,
	sampinfo = list(
		ind.obs = ind,
		ind.RespNA = dataXY.out$indvecs$ind.RespNA,
		sample.size = nobs,
		obs.sample.size = n.allxy,
		missing.sample.size = nobs - n.allxy,
		rankX = p,
		z = zt,
		trialsvec = trialsvec,
		X = X2,
		effnames = dataXY.out$Xmats$effnames,
		setzero = dataXY.out$indvecs$setzero,
		setNA = dataXY.out$indvecs$setNA,
		setNA2 = dataXY.out$indvecs$setNA2,
		cutX1toX2 = dataXY.out$indvecs$cutX1toX2,
		REs = REs,
                REmodelmatrices = REmodelmatrices
	),
	estimates = list(
		theta = parmest,
 		nugget = nugget,
		V = V,
		Vi = Vi,
		betahat = b.hat,
		covb = covb,
		covbi = covbi,
		m2LL = m2LL,
		Warnlog = Warnlog
	),
	loglik.surface=get("RESULT",loglik.environment),
		optimOutput=parmest.out
	)
	class(outpt) <- "glmssn"
	outpt
}

