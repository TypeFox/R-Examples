# Copyright 2006-8 Roger Bivand
# 

ggwr <- function(formula, data = list(), coords, bandwidth, 
	gweight=gwr.Gauss, adapt=NULL, fit.points, family=gaussian,
	longlat=NULL, type=c("working", "deviance", "pearson", "response")) {
	type <- match.arg(type)
	resid_name <- paste(type, "resids", sep="_")
	this.call <- match.call()
	p4s <- as.character(NA)
	Polys <- NULL
	if (is(data, "SpatialPolygonsDataFrame")) 
		Polys <- as(data, "SpatialPolygons")
	if (is(data, "Spatial")) {
		if (!missing(coords))
		    warning("data is Spatial* object, ignoring coords argument")
		coords <- coordinates(data)
		p4s <- proj4string(data)
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
	if (is.null(colnames(coords))) 
		colnames(coords) <- c("coord.x", "coord.y")

    	if (is.character(family)) 
            family <- get(family, mode = "function", envir = parent.frame())
    	if (is.function(family)) family <- family()
    	if (is.null(family$family)) {
            print(family)
            stop("'family' not recognized")
    	}
     	if (missing(data)) data <- environment(formula)
    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data"), names(mf), 0)
    	mf <- mf[c(1, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())

	mt <- attr(mf, "terms")
	y <- model.extract(mf, "response") # extract X and Y
	x <- model.matrix(mt, mf)

	offset <- model.offset(mf)
	if (!is.null(offset) && length(offset) != NROW(y))
	    stop("number of offsets should equal number of observations")
	if (is.null(offset)) offset <- rep(0, length(c(y)))

	glm_fit <- glm.fit(x=x, y=y, offset=offset, family=family)


	if (missing(fit.points)) {
		fp.given <- FALSE
		fit.points <- coords
		colnames(fit.points) <- colnames(coords)
	} else fp.given <- TRUE
	griddedObj <- FALSE
	if (is(fit.points, "Spatial")) {
		Polys <- NULL
		if (is(fit.points, "SpatialPolygonsDataFrame")) {
			Polys <- Polygons(fit.points)
			fit.points <- coordinates(fit.points)
		} else {
			griddedObj <- gridded(fit.points)
			fit.points <- coordinates(fit.points)
		}
	}

	n <- NROW(fit.points)
	rownames(fit.points) <- NULL
	if (is.null(colnames(fit.points))) colnames(fit.points) <- c("x", "y")
	m <- NCOL(x)
	if (NROW(x) != NROW(coords))
		stop("Input data and coordinates have different dimensions")
	if (is.null(adapt)) {
		if (!missing(bandwidth)) {
			bw <- bandwidth
			bandwidth <- rep(bandwidth, n)
		} else stop("Bandwidth must be given for non-adaptive weights")
	} else {
		bandwidth <- gw.adapt(dp=coords, fp=fit.points, quant=adapt,
			longlat=longlat)
		bw <- bandwidth
	}
	if (any(bandwidth < 0)) stop("Invalid bandwidth")
	gwr.b <- matrix(nrow=n, ncol=m)
#	assign(resid_name, numeric(n))
	v_resids <- numeric(n)
	colnames(gwr.b) <- colnames(x)
	lhat <- NA
	sum.w <- numeric(n)
	dispersion <- numeric(n)
	for (i in 1:n) {
		dxs <- spDistsN1(coords, fit.points[i,], longlat=longlat)
		if (any(!is.finite(dxs)))
			dxs[which(!is.finite(dxs))] <- .Machine$double.xmax/2
		w.i <- gweight(dxs^2, bandwidth[i])
		if (any(w.i < 0 | is.na(w.i)))
        		stop(paste("Invalid weights for i:", i))
		lm.i <- glm.fit(y=y, x=x, weights=w.i, offset=offset,
			family=family)
		sum.w[i] <- sum(w.i)
		gwr.b[i,] <- coefficients(lm.i)
		if (!fp.given) v_resids[i] <- residuals.glm(lm.i, type=type)[i]
		else is.na(v_resids[i]) <- TRUE
    		df.r <- lm.i$df.residual
        	if (lm.i$family$family %in% c("poisson", "binomial")) 
            		dispersion[i] <- 1
        	else {
			if (df.r > 0) {
            			dispersion[i] <- sum((lm.i$weights * 
			    	    lm.i$residuals^2)[lm.i$weights >  0])/df.r
        		} else {
            			dispersion[i] <- NaN
			}
        	}
	}
	df <- data.frame(sum.w=sum.w, gwr.b, dispersion=dispersion)
	df[[resid_name]] <- v_resids

	SDF <- SpatialPointsDataFrame(coords=fit.points, 
			data=df, proj4string=CRS(p4s))

	if (griddedObj) {
		gridded(SDF) <- TRUE
	} else {
		if (!is.null(Polys)) {
			df <- data.frame(SDF@data)
			rownames(df) <- sapply(slot(Polys, "polygons"),
                            function(i) slot(i, "ID"))
			SDF <- SpatialPolygonsDataFrame(Sr=Polys, data=df)
		}
	}
	z <- list(SDF=SDF, lhat=lhat, lm=glm_fit, results=NULL, 
		bandwidth=bw, adapt=adapt, hatmatrix=FALSE, 
		gweight=deparse(substitute(gweight)), fp.given=fp.given,
                this.call=this.call)
	class(z) <- "gwr"
	invisible(z)
}


