CCmaps <- function(obj, zcol=NULL, cvar=NULL, cvar.names=NULL,
	..., names.attr, scales = list(draw = FALSE), 
	xlab = NULL, ylab = NULL, aspect = mapasp(obj, xlim, ylim), 
	sp.layout = NULL, xlim = bbox(obj)[1,], ylim = bbox(obj)[2,]) {
	stopifnot(is(obj, "SpatialPolygonsDataFrame")) 
	stopifnot(!is.null(zcol), !is.null(cvar))
	n <- length(slot(obj, "polygons"))
	stopifnot(length(zcol) == 1L)
	ncc <- length(cvar)
	stopifnot(ncc <= 2, ncc > 0)
	if (is.null(cvar.names)) cvar.names <- names(cvar)
	nlcc <- integer(ncc)
	lcc <- vector(mode="list", length=ncc)
#	fcc <- logical(nlcc)
	fcc <- logical(ncc)
	for (i in 1:ncc) {
	    ccc <- class(cvar[[i]])
	    stopifnot(ccc %in% c("factor", "shingle"))
	    fcc[i] <- ccc == "factor"
	    stopifnot(length(cvar[[i]]) == n)
	    nlcc[i] <- nlevels(cvar[[i]])
	    lcc[[i]] <- levels(cvar[[i]])
	}
	obj <- obj[zcol]
	zcol <- names(obj)
	Outside <- function(x, y, z) (x < y | x > z)
	if (ncc == 1) {
		if (fcc[1]) {
	            for (j in 1:nlcc[1]) {
		        vn <- paste(cvar.names[1], lcc[[1]][j], sep="_")
		        io <- as.character(cvar[[1]]) != lcc[[1]][j]
		        obj[[vn]] <- obj[[zcol]]
		        is.na(obj[[vn]]) <- io
	            }
   	        } else {
		    ilcc <- do.call("rbind", lcc[[1]])
	            for (j in 1:nlcc[1]) {
		        vn <- paste(cvar.names[1], j, sep="_")
			io <- Outside(cvar[[1]], ilcc[j,1], ilcc[j,2]) 
		        obj[[vn]] <- obj[[zcol]]
		        is.na(obj[[vn]]) <- io
		    }
	    }
	    nms <- names(obj)
	    nms <- nms[-(match(zcol, nms))]
	    if (fcc[1]) {
	        print(spplot(obj, zcol=nms, ..., scales = scales, 
		    xlab = xlab, ylab = ylab, aspect = aspect, 
		    sp.layout = sp.layout, xlim = xlim, ylim = ylim, 
		    strip=strip.custom(which.given=1, 
		    factor.levels=lcc[[1]], par.strip.text=list(cex=0.8), 
		    bg="grey95")))
	    } else {
	        print(spplot(obj, zcol=nms, ..., scales = scales, 
		    xlab = xlab, ylab = ylab, aspect = aspect, 
		    sp.layout = sp.layout, xlim = xlim, ylim = ylim, 
		    strip=strip.custom(which.given=1, 
		    shingle.intervals=as.matrix(lcc[[1]]), 
		    var.name=cvar.names[1], par.strip.text=list(cex=0.8), 
		    bg="grey95", fg="grey75")))
	    }
	} else {
	    if (all(fcc)) {
	        for (i in 1:nlcc[1]) {
		    for (j in 1:nlcc[2]) {
		        vn <- paste(cvar.names[1], lcc[[1]][i], cvar.names[2], 
			    lcc[[2]][j], sep="_")
		        obj[[vn]] <- obj[[zcol]]
			ioi <- as.character(cvar[[1]]) != lcc[[1]][i] 
			ioj <- as.character(cvar[[2]]) != lcc[[2]][j]
			io <- ioi | ioj
		        is.na(obj[[vn]]) <- io
		    }
		}
	        nms <- names(obj)
	        nms <- nms[-(match(zcol, nms))]
		lcc1 <- lcc[[1]]
		xlcc <- NULL
		for (i in 1:nlcc[1]) {
		    xlcc <- c(xlcc, rep(lcc1[i], nlcc[2]))
		}
		lcc2 <- lcc[[2]]
		xlcc2 <- rep(lcc2, nlcc[1])
	        print(spplot(obj, zcol=nms, ..., scales = scales, 
		    xlab = xlab, ylab = ylab, aspect = aspect, 
		    sp.layout = sp.layout, xlim = xlim, ylim = ylim, 
		    strip=strip.custom(which.given=1,
		    factor.levels=xlcc, 
		    par.strip.text=list(cex=0.8), bg="grey95"), 
		    strip.left=strip.custom(which.given=1, 
		    factor.levels=xlcc2, 
		    par.strip.text=list(cex=0.8), bg="grey95")))
	    } else if (any(fcc)) {
	      if (fcc[1]) {
		jlcc <- do.call("rbind", lcc[[2]])
	        for (i in 1:nlcc[1]) {
		    for (j in 1:nlcc[2]) {
		        vn <- paste(cvar.names[1], lcc[[1]][i], cvar.names[2], 
			    j, sep="_")
		        obj[[vn]] <- obj[[zcol]]
			ioi <- as.character(cvar[[1]]) != lcc[[1]][i] 
			ioj <- Outside(cvar[[2]], jlcc[j,1], jlcc[j,2])
			io <- ioi | ioj
		        is.na(obj[[vn]]) <- io
		    }
		}
	        nms <- names(obj)
	        nms <- nms[-(match(zcol, nms))]
		lcc1 <- lcc[[1]]
		xlcc <- NULL
		for (i in 1:nlcc[1]) {
		    xlcc <- c(xlcc, rep(lcc1[i], nlcc[2]))
		}
		lcc2 <- matrix(unlist(lcc[[2]]), ncol=2, byrow=TRUE)
		xlcc2 <- matrix(rep(t(lcc2), nlcc[1]), byrow=TRUE, ncol=2)
	        print(spplot(obj, zcol=nms, ..., scales = scales, 
		    xlab = xlab, ylab = ylab, aspect = aspect, 
		    sp.layout = sp.layout, xlim = xlim, ylim = ylim, 
		    strip=strip.custom(which.given=1,
		    factor.levels=xlcc, 
		    par.strip.text=list(cex=0.8), bg="grey95"), 
		    strip.left=strip.custom(which.given=1, 
		    shingle.intervals=xlcc2, var.name=cvar.names[2], 
		    par.strip.text=list(cex=0.8), bg="grey95", fg="grey75")))
	      } else {
		ilcc <- do.call("rbind", lcc[[1]])
	        for (i in 1:nlcc[1]) {
		    for (j in 1:nlcc[2]) {
		        vn <- paste(cvar.names[1], i, cvar.names[2], 
			    lcc[[2]][j], sep="_")
		        obj[[vn]] <- obj[[zcol]]
			ioi <- Outside(cvar[[1]], ilcc[i,1], ilcc[i,2])
			ioj <- as.character(cvar[[2]]) != lcc[[2]][j]
			io <- ioi | ioj
		        is.na(obj[[vn]]) <- io
		    }
		}
	        nms <- names(obj)
	        nms <- nms[-(match(zcol, nms))]
		lcc1 <- matrix(unlist(lcc[[1]]), ncol=2, byrow=TRUE)
		xlcc <- matrix(ncol=2)
		for (i in 1:nlcc[1]) {
		    xlcc <- rbind(xlcc, matrix(rep(lcc1[i,], nlcc[2]), 
			ncol=2, byrow=TRUE))
		}
		xlcc <- xlcc[-1,]
		lcc2 <- lcc[[2]]
		xlcc2 <- rep(lcc2, nlcc[1])
	        print(spplot(obj, zcol=nms, ..., scales = scales, 
		    xlab = xlab, ylab = ylab, aspect = aspect, 
		    sp.layout = sp.layout, xlim = xlim, ylim = ylim, 
		    strip=strip.custom(which.given=1, 
		    shingle.intervals=xlcc, var.name=cvar.names[1], 
		    par.strip.text=list(cex=0.8), bg="grey95", fg="grey75"), 
		    strip.left=strip.custom(which.given=1,
		    factor.levels=xlcc2, 
		    par.strip.text=list(cex=0.8), bg="grey95")))
	      }
	    } else {
		ilcc <- do.call("rbind", lcc[[1]])
		jlcc <- do.call("rbind", lcc[[2]])
	        for (i in 1:nlcc[1]) {
		    for (j in 1:nlcc[2]) {
		        vn <- paste(cvar.names[1], i, cvar.names[2], j, sep="_")
		        obj[[vn]] <- obj[[zcol]]
			ioi <- Outside(cvar[[1]], ilcc[i,1], ilcc[i,2]) 
			ioj <- Outside(cvar[[2]], jlcc[j,1], jlcc[j,2])
			io <- ioi | ioj
		        is.na(obj[[vn]]) <- io
		    }
		}
	        nms <- names(obj)
	        nms <- nms[-(match(zcol, nms))]
		lcc1 <- matrix(unlist(lcc[[1]]), ncol=2, byrow=TRUE)
		xlcc <- matrix(ncol=2)
		for (i in 1:nlcc[1]) {
		    xlcc <- rbind(xlcc, matrix(rep(lcc1[i,], nlcc[2]), 
			ncol=2, byrow=TRUE))
		}
		xlcc <- xlcc[-1,]
		lcc2 <- matrix(unlist(lcc[[2]]), ncol=2, byrow=TRUE)
		xlcc2 <- matrix(rep(t(lcc2), nlcc[1]), byrow=TRUE, ncol=2)
	        print(spplot(obj, zcol=nms, ..., scales = scales, 
		    xlab = xlab, ylab = ylab, aspect = aspect, 
		    sp.layout = sp.layout, xlim = xlim, ylim = ylim, 
		    strip=strip.custom(which.given=1,
		    shingle.intervals=xlcc, var.name=cvar.names[1], 
		    par.strip.text=list(cex=0.8), bg="grey95", fg="grey75"), 
		    strip.left=strip.custom(which.given=1, 
		    shingle.intervals=xlcc2, var.name=cvar.names[2], 
		    par.strip.text=list(cex=0.8), bg="grey95", fg="grey75")))
	    }
	}
	invisible(obj)
}

