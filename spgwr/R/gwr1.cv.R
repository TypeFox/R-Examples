# Copyright 2006-8 Roger Bivand
# 

ggwr.sel <- function(formula, data = list(), coords, 
	adapt=FALSE, gweight=gwr.Gauss, family=gaussian, verbose=TRUE, 
	longlat=NULL, RMSE=FALSE, tol=.Machine$double.eps^0.25) {
	if (!is.logical(adapt)) stop("adapt must be logical")
	if (is(data, "Spatial")) {
		if (!missing(coords))
		    warning("data is Spatial* object, ignoring coords argument")
		coords <- coordinates(data)
                if (is.null(longlat) || !is.logical(longlat)) {
	            if (!is.na(is.projected(data)) && !is.projected(data)) {
                        longlat <- TRUE
                    } else {
                        longlat <- FALSE
                    }
                }
		data <- as(data, "data.frame")
	}
        if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
	if (missing(coords))
		stop("Observation coordinates have to be given")
    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data"), names(mf), 0)
    	mf <- mf[c(1, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())

	mt <- attr(mf, "terms")
	y <- model.extract(mf, "response") 

	if (!adapt) {
		bbox <- cbind(range(coords[,1]), range(coords[,2]))
		difmin <- spDistsN1(bbox, bbox[2,], longlat)[1]
		if (any(!is.finite(difmin)))
			difmin[which(!is.finite(difmin))] <- 0
		beta1 <- difmin/1000
		beta2 <- difmin
		opt <- optimize(ggwr.cv.f, lower=beta1, upper=beta2, 
			maximum=FALSE, formula=formula, data=data, 
			family=family, coords=coords, y=y,
			gweight=gweight, verbose=verbose, longlat=longlat, 
			RMSE=RMSE, tol=tol)
		bdwt <- opt$minimum
		res <- bdwt
	} else {
		beta1 <- 0
		beta2 <- 1
		opt <- optimize(ggwr.cv.adapt.f, lower=beta1, 
			upper=beta2, maximum=FALSE, formula=formula, data=data, 
			family=family, coords=coords, y=y,
			gweight=gweight, verbose=verbose, longlat=longlat, 
			RMSE=RMSE, tol=tol)
		q <- opt$minimum
		res <- q
	}
	res
}


ggwr.cv.f <- function(bandwidth, formula, data, family, coords, y,
	gweight, verbose=TRUE, longlat=FALSE, RMSE=FALSE) {
    n <- nrow(coords)
    cv <- numeric(n)
    options(show.error.messages = FALSE)
    for (i in 1:n) {
	dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
	if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
	w.i <- gweight(dxs^2, bandwidth)
        w.i[i] <- 0
	if (any(w.i < 0 | is.na(w.i)))
       		stop(paste("Invalid weights for i:", i))
	datai <- data.frame(data, w.i=w.i)
        lm.i <- try(glm(formula=formula, data=datai, family=family, 
		weights=w.i))
        if(!inherits(lm.i, "try-error")) {
            cv[i] <- y[i] - predict(lm.i, newdata=data[i,,drop=FALSE],
		type="response")
        }
    }
    score <- sum(t(cv) %*% cv)
    if (RMSE) score <- sqrt(score/n)
#    score <- sqrt(sum(t(cv) %*% cv)/n)
    options(show.error.messages = TRUE)
    if (verbose) cat("Bandwidth:", bandwidth, "CV score:", score, "\n")
    score
}


ggwr.cv.adapt.f <- function(q, formula, data, family, coords, y, gweight, 
	verbose=TRUE, longlat=FALSE, RMSE=FALSE) {
    n <- nrow(coords)
    cv <- numeric(n)
    bw <- gw.adapt(dp=coords, fp=coords, quant=q, longlat=longlat)
    options(show.error.messages = FALSE)
    for (i in 1:n) {
	dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
	if (!is.finite(dxs[i])) dxs[i] <- .Machine$double.xmax/2
	w.i <- gweight(dxs^2, bw[i])
        w.i[i] <- 0
	if (any(w.i < 0 | is.na(w.i)))
       		stop(paste("Invalid weights for i:", i))
	datai <- data.frame(data, w.i=w.i)
        lm.i <- try(glm(formula=formula, data=datai, family=family, 
		weights=w.i))
        if(!inherits(lm.i, "try-error")) {
            cv[i] <- y[i] - predict(lm.i, newdata=data[i,,drop=FALSE],
		type="response")
        }
    }
    score <- sum(t(cv) %*% cv)
    if (RMSE) score <- sqrt(score/n)
#    score <- sqrt(sum(t(cv) %*% cv)/n)
    options(show.error.messages = TRUE)
    if (verbose) cat("Adaptive q:", q, "CV score:", score, "\n")
    score
}

