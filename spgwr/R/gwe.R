# Copyright 2001-2004 Roger Bivand and Danlin Yu
# 

#gw.dists <- function(dp, pt, lonlat=FALSE) {
#	n <- nrow(dp)
#	dists <- numeric(n)
#	res <- .C("gw_dists", as.double(dp[,1]), as.double(dp[,2]),
#		as.double(pt[1]), as.double(pt[2]), as.integer(n),
#		as.double(dists), as.integer(lonlat), PACKAGE="spgwr")[[6]]
#	res
#}

gw.adapt <- function(dp, fp, quant, longlat=FALSE) {
	n1 <- nrow(dp)
	n2 <- nrow(fp)
        storage.mode(dp) <- "double"
        storage.mode(fp) <- "double"
	dists <- numeric(n1)
	bw <- numeric(n2)
	if (quant > 1) {
		factor <- quant
		quant <- 1
		res <- .C("gw_adapt", dp[,1], dp[,2],
			fp[,1], fp[,2], as.integer(n1),
			as.integer(n2), as.double(bw), as.double(quant),
			dists, as.integer(longlat), 
			PACKAGE="spgwr")[[7]]
		res <- res*factor
	} else {
		n.ideal <- n1*quant
		n.lower <- floor(n.ideal)
		n.higher <- n.lower+1
		q1 <- n.lower/n1
		q2 <- n.higher/n1
		res1 <- .C("gw_adapt", dp[,1], dp[,2],
			fp[,1], fp[,2], as.integer(n1),
			as.integer(n2), as.double(bw), as.double(q1),
			dists, as.integer(longlat), 
			PACKAGE="spgwr")[[7]]
		res2 <- .C("gw_adapt", dp[,1], dp[,2],
			fp[,1], fp[,2], as.integer(n1),
			as.integer(n2), as.double(bw), as.double(q2),
			dists, as.integer(longlat), 
			PACKAGE="spgwr")[[7]]

		res <- (n.ideal - n.lower)*res2 + (n.higher - n.ideal)*res1
	}
	res
}


gw.cov <- function(x, vars, fp, adapt=NULL, bw, gweight=gwr.bisquare, 
		cor=TRUE, var.term=FALSE, longlat=NULL) {
	p4s <- as.character(NA)
	Polys <- NULL
	fp.missing <- missing(fp)
	if (is(x, "SpatialPolygonsDataFrame")) {
		Polys <- as(x, "SpatialPolygons")
		gridded <- gridded(x)
		dp <- coordinates(x)
		p4s <- proj4string(x)
		data <- as(x, "data.frame")
                if (is.null(longlat) || !is.logical(longlat)) {
	            if (!is.na(is.projected(x)) && !is.projected(x)) {
                        longlat <- TRUE
                    } else {
                        longlat <- FALSE
                    }
                }
	} else if (is(x, "SpatialPointsDataFrame")) {
		gridded <- gridded(x)
		dp <- coordinates(x)
		p4s <- proj4string(x)
		data <- as(x, "data.frame")
                if (is.null(longlat) || !is.logical(longlat)) {
	            if (!is.na(is.projected(data)) && !is.projected(data)) {
                        longlat <- TRUE
                    } else {
                        longlat <- FALSE
                    }
                }
	} else stop("x must be a Spatial Polygons or Points DataFrame")
        if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
        if (is.integer(vars)) vars <- names(x)[vars]
        stopifnot(is.character(vars))
	x <- as.matrix(data[, vars])
	if (any(is.na(x))) stop("x contains NAs")
	nc <- ncol(x)
	if (is.null(colnames(x))) colnames(x) <- paste("V", 1:nc, sep="")
	cn <- colnames(x)
	n1 <- nrow(dp)
	if (n1 != nrow(x)) stop ("Differing row numbers between x and dp")

	# prepare bandwidths for data points
	if (is.null(adapt)) {
		if (!missing(bw)) bw0 <- rep(bw, n1)
		else stop("Bandwidth must be given for non-adaptive weights")
	}
	else bw0 <- gw.adapt(dp=dp, fp=dp, quant=adapt, longlat=longlat)
	dm <- numeric(nc)
	rss <- numeric(nc)
	trhat <- 0
	for (i in 1:n1) { # establish residuals for data points and 
			# calculate hat matrix trace
		dxs <- spDistsN1(dp, dp[i,], longlat=longlat)
		if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
		wts <- gweight(dist2=dxs^2, bw0[i])
		if (any(wts < 0 | is.na(wts))) {
			print(dxs)
			print(wts)
        		stop(paste("Invalid weights for i:", i)) }
		if (sum(wts) == 0) {
			print(dxs)
			print(wts)
        		stop(paste("Invalid weights for i:", i)) }
		for (j in 1:nc) {
			dm[j] <- weighted.mean(x[,j], wts)
			rss[j] <- rss[j] + (x[i,j] - dm[j])^2
		}
		trhat <- trhat + wts[i]/sum(wts)
	}
	dm <- sqrt(rss/(n1 - trhat))
	# global adjusted residual standard error

	if (missing(fp)) {
#		fp.given <- FALSE
		fp <- dp
		colnames(fp) <- colnames(dp)
	} else {
#		fp.given <- TRUE
                stopifnot(is(fp, "Spatial"))
		gridded <- gridded(fp)
	}
        fp <- coordinates(fp)

	# prepare bandwidths for fitting/estimation points
        
	n2 <- nrow(fp)
	if (is.null(adapt)) {
		if (!missing(bw)) bw <- rep(bw, n2)
		else stop("Bandwidth must be given for non-adaptive weights")
	}
	else bw <- gw.adapt(dp=dp, fp=fp, quant=adapt, longlat=longlat)
	
	ut <- ((nc^2 -nc)/2)
	res <- matrix(NA, nrow=n2, ncol=((4*nc)+(ifelse(cor, 2*ut, ut))))
	rnm <- paste("mean", cn, sep=".")
	rnm <- c(rnm, paste("sd", cn, sep="."))
	rnm <- c(rnm, paste("sem", cn, sep="."))
	rnm <- c(rnm, paste("diff", cn, sep="."))
	if (nc > 1) {
		corn <- outer(cn, cn, paste, sep=",")
		rnm <- c(rnm, 
			paste("cov(", corn[upper.tri(corn)], ")", sep=""))
		if (cor) rnm <- c(rnm, 
			paste("cor(", corn[upper.tri(corn)], ")", sep=""))
	}
	
	colnames(res) <- rnm
	swts <- numeric(n2)
	swts2 <- numeric(n2)

	for (i in 1:n2) {
		dxs <- spDistsN1(dp, fp[i,], longlat=longlat)
		if (any(!is.finite(dxs)))
			dxs[which(!is.finite(dxs))] <- 0
#		if (!is.finite(dxs[i])) dxs[i] <- 0
		wts <- gweight(dxs^2, bw[i])
		if (any(wts < 0 | is.na(wts))) {
			print(dxs)
			print(wts)
        		stop(paste("Invalid weights for i:", i))}
		swts[i] <- sum(wts)
		if (swts[i] == 0) {
			print(dxs)
			print(wts)
        		stop(paste("Invalid weights for i:", i))}
		swts2[i] <- sum((wts/swts[i])^2)
		res1 <- cov.wt(as.matrix(x), wts, cor=cor)
		res[i, 1:nc] <- res1$center
		if (var.term) {
			cov <- res1$cov * sqrt(1 - swts2[i])
		} else {
			cov <- res1$cov
		}
		sd <- sqrt(diag(cov))
		res[i, (nc+1):(2*nc)] <- sd
		if (nc > 1) {
			res[i, ((4*nc)+1):((4*nc)+ut)] <- 
				res1$cov[upper.tri(res1$cov)]
			if (cor) {
			        sdinv <- diag(1/sd, nrow(cov))
        			corr <- sdinv %*% cov %*% sdinv
				res[i, ((4*nc)+ut+1):((4*nc)+2*ut)] <- 
					corr[upper.tri(corr)]
			}
		}
	}
	gxbar <- apply(x, 2, mean)
	means <- grep("mean", colnames(res))
	for (i in means) res[,((nc*2)+i)] <- dm[i] * sqrt(swts2)
	for (i in means) 
		res[,((nc*3)+i)] <- (gxbar[i] - res[,i]) / res[,((nc*2)+i)]
	SDF <- SpatialPointsDataFrame(coords=fp, data=data.frame(res, 
		fp), proj4string=CRS(p4s))
	if (gridded) gridded(SDF) <- TRUE
	else if (!is.null(Polys) && fp.missing) {
		df <- data.frame(SDF@data)
		rownames(df) <- sapply(slot(Polys, "polygons"),
                    function(i) slot(i, "ID"))
		SDF <- SpatialPolygonsDataFrame(Sr=Polys, data=df)
	}
	res <- list(SDF=SDF, bandwidth=bw, adapt=adapt,  
		gweight=deparse(substitute(gweight)))
	class(res) <- c("gw.cov", "matrix")
	attr(res, "swts") <- swts
	attr(res, "sdswts2") <- swts2
	attr(res, "g.se") <- dm
	res
}

