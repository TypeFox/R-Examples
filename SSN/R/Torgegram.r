Torgegram <- function(object, ResponseName,
	maxlag = NULL, nlag = 6, inc = 0, nlagcutoff = 15, 
	EmpVarMeth = "MethMoment")
{
	if(class(object) == "influenceSSN") object <- object$ssn.object
	data <- object@obspoints@SSNPoints[[1]]@point.data
	data <- cbind(data, object@obspoints@SSNPoints[[1]]@point.coords)
	xcol <- "coords.x1"
	ycol <- "coords.x2"
	distord <- order(data[,"netID"],data[,"pid"])
	names(distord) <- rownames(data)[distord]

	nIDs <- sort(as.integer(as.character(unique(data[,"netID"]))))
	dist.junc <- matrix(0, nrow = length(data[,1]), ncol = length(data[,1]))
	Dif2s <- NULL
	sqrtDifs <- NULL
	Dists <- NULL
	FCons <- NULL
	Covp <- NULL
	mnz <- mean(data[,ResponseName], na.rm = TRUE)
	nsofar <- 0
	for(i in nIDs) {
		workspace.name <- paste("dist.net", i, ".RData", sep = "")
		path <- file.path(object@path, "distance", "obs", 
			workspace.name)
		if(!file.exists(path)) {
			stop("Unable to locate required distance matrix")
		}
		file_handle <- file(path, open="rb")
		distmat <- unserialize(file_handle)
		close(file_handle)
		ordi <- order(as.numeric(rownames(distmat)))
		distmat <- distmat[ordi, ordi, drop = F]
		Ds <- distmat + t(distmat)
		FCs <- 1 - (pmin(distmat,t(distmat)) > 0)*1
		zs <- data[data$pid %in% rownames(distmat), 
			ResponseName, drop = F]
		attr(zs,"pid") <- data[data$pid %in% rownames(distmat),"pid"]
		zs <- zs[order(attr(zs,"pid")),ResponseName]
		ind <- !is.na(zs)
		ni <- sum(ind)
		zs <- zs[ind, drop = F]
		Ds <- Ds[ind,ind, drop = F]
		FCs <- FCs[ind,ind, drop = F]
		rs <- zs - mnz
		CP <- rs%o%rs
		diff <- abs(zs%o%rep(1, times = ni) - 
			rep(1, times = ni)%o%zs)
		diff2 <- diff^2
		sqrtd <- sqrt(diff)
		Dif2s <- c(Dif2s, as.vector(diff2[col(diff2) < row(diff2)]))
		sqrtDifs <- c(sqrtDifs, as.vector(sqrtd[col(sqrtd) < row(sqrtd)]))
		Dists <- c(Dists, as.vector(Ds[col(Ds) < row(Ds)]))
		FCons <- c(FCons, as.vector(FCs[col(FCs) < row(FCs)]))
		Covp <- c(Covp, as.vector(CP[col(CP) < row(CP)]))
	}
	if(is.null(maxlag)) maxlag <- max(Dists)
	indmax <- Dists <= maxlag
	Dists <- Dists[indmax]
	Dif2s <- Dif2s[indmax]
	sqrtDifs <- sqrtDifs[indmax]
	FCons <- FCons[indmax]
	Covp <- Covp[indmax]

	if(inc <= 0) inc <- maxlag/nlag

	store.results <- matrix(data = NA, ncol = 6, nrow = nlag,
		dimnames = list(1:nlag, c("distance.connect", "gam.connect", "np.connect",
		"distance.unconnect", "gam.unconnect","np.unconnect")))
	for ( i in 1:nlag){
		ind1 <- (Dists >= (i-1)*inc) & (Dists < i*inc) & (FCons == 1)
		ind2 <- (Dists >= (i-1)*inc) & (Dists < i*inc) & (FCons == 0)
		nclass1 <- sum(ind1)
		nclass2 <- sum(ind2)
		if(EmpVarMeth == "MethMoment") {
			cv1 <- mean(Dif2s[ind1])
			cv2 <- mean(Dif2s[ind2])
		}
		if(EmpVarMeth == "RobustMean") {
			cv1 <- ((mean(sqrtDifs[ind1]))^4)/(.457 + .494/sum(ind1))
			cv2 <- ((mean(sqrtDifs[ind2]))^4)/(.457 + .494/sum(ind2))
		}
		if(EmpVarMeth == "RobustMedian") {
			cv1 <- (median(sqrtDifs[ind1]))^4/.457
			cv2 <- (median(sqrtDifs[ind2]))^4/.457
		}
		if(EmpVarMeth == "Covariance" | EmpVarMeth == "CovMean") {
			cv1 <- mean(Covp[ind1])
			cv2 <- mean(Covp[ind2])
		}
		mean.dis1 <- mean(Dists[ind1])
		mean.dis2 <- mean(Dists[ind2])
		if(nclass1 > 0 | nclass2 > 0) store.results[i,] <-
				c(mean.dis1, cv1, nclass1, mean.dis2, cv2, nclass2)
	}
	store.results[,"gam.connect"] <- store.results[,"gam.connect"]/2
	store.results[,"gam.unconnect"] <- store.results[,"gam.unconnect"]/2
	indc <- store.results[,"np.connect"] >= nlagcutoff 
	store.results[!indc,1:3] <- NA	
	indu <- store.results[,"np.unconnect"] >= nlagcutoff
	store.results[!indu,4:6] <- NA	

	ev <- as.data.frame(store.results)

	ev[is.nan(ev[,"distance.unconnect"]),"distance.unconnect"] <- NA
	ev[is.nan(ev[,"gam.unconnect"]),"gam.unconnect"] <- NA
	ev[is.nan(ev[,"distance.connect"]),"distance.connect"] <- NA
	ev[is.nan(ev[,"gam.connect"]),"gam.connect"] <- NA

	ev <- as.data.frame(ev)
	ev <- as.list(ev)
	for(i in 1:length(ev))
		ev[[i]] <- ev[[i]][!is.na(ev[[i]])]
	ev$call <- as.list(match.call()[-1])
	class(ev) <- "Torgegram"
	ev
}


