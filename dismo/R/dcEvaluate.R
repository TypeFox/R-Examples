


dcEvaluate <- function(p, a, reference, lonlat=TRUE, binsize=15, predp, preda, model, predictors, fun=predict) {
	
	if (missing(predp)) {
		p <- stats::na.omit(p)
		a <- stats::na.omit(a)
	} else {
		i <- is.na(p)
		if (any(i)) {
			p <- p[!i,]
			predp <- predp[!i,]
		}
		i <- is.na(a)
		if (any(i)) {
			a <- a[!i,]
			preda <- preda[!i,]
		}			
	}
	reference <- stats::na.omit(reference)
	
	dp <- apply(pointDistance(p, reference, longlat=lonlat), 1, min) / 1000
	da <- apply(pointDistance(a, reference, longlat=lonlat), 1, min) / 1000
	#if (is.null(dist)) {
	#	dist <- c(0, as.vector(quantile(dp,  probs=0:10/10)))
	#} else {
	#	if (length(dist) == 1) {
	#		n <- max(1, max(da) / dist)
	#		dist <- 0:n * dist
	#	}
	#}
	
	n <- round(length(dp) / binsize)
	dist <- c(0, as.vector(quantile(dp,  probs=0:n/n)))
	
	if (missing(predp)) {
		pv <- extract(predictors, p)
		pa <- extract(predictors, a)
		predp <- fun(model, pv)
		preda <- fun(model, pa)
	}
	e <- list()
	pwd <- TRUE
	for (d in 1:(length(dist)-1)) {
		
		if (pwd) {
			i <- which(dp > dist[d])
			ab <- pwdSample(p[i, ], a, reference, lonlat=lonlat, warn=TRUE) 
			i <- i[!is.na(ab)]
			j <- stats::na.omit(ab)
			abss <- preda[j]
			pres <- predp[i]
		} else {
			abss <- preda[da > dist[d] & da <= dist[d+1]]
			pres <- predp[dp > dist[d] & dp <= dist[d+1]]
		}
		if (NROW(pres) > 1 & NROW(abss) > 1) {
			e[[d]] <- evaluate(pres, abss)
		} else {
			e[[d]] <- NULL
		}
	}
	dist <- (dist[-1] + dist[-length(dist)]) / 2
	names(e) <- dist[2:(length(e)+1)]
	e
}
