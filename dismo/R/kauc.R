
.kAUC <- function(pxy, axy, ref, lonlat=TRUE, pv, av, model, null, groups=25, mins=10) {
	stopifnot(nrow(pxy)==nrow(pv))
	stopifnot(nrow(axy)==nrow(av))
	
	u <- ssb(pxy, axy, ref, lonlat=lonlat, avg=FALSE)
	up <- u[[1]]
	ua <- u[[2]]

	ret <- FALSE
	if (length(up) < 2* mins) { warning('not enough p samples') ; ret=TRUE}
	if (length(ua) < 2*mins) { warning('not enough a samples') ; ret=TRUE}
	if (ret) {
		res <- rep(NA, length(model))
		return(cbind(n=NA, distance=NA, res))
	}
	
	n <- min(length(up)/mins, groups)
	q <- quantile(up, 1:n / n)
	
	up <- cbind(up, 1)
	ua <- cbind(ua, 1)
	
	for (i in 1:length(q)) {
		up[up[,1] > q[i], 2] <- i+1
		ua[ua[,1] > q[i], 2] <- i+1
	}
	
	em1 <- vector(length=length(q), mode='list')
	em <- list()
	
	if (!is.list(model)) {
		model <- list(model)
	}
	lm <- length(model)
	
	for (k in 1:lm) {
		em[[k]] <- em1
	}
	en <- vector(length=length(q), mode='list')
	sb <- vector(length=length(q), mode='list')
	
	for (i in 1:length(q)) {
		npt <- pxy[up[,2]==i, ]
		nat <- axy[ua[,2]==i, ]

		if (nrow(npt) > 5 & nrow(nat) > 5) {
			#s <- pwdSample(npt[,1:2], nat[,1:2], ref, lonlat=lonlat, warn=FALSE)
			s <- pwdSample(npt[,1:2], nat[,1:2], ref, lonlat=lonlat)
			npt <- npt[!is.na(s), ]
			nat <- nat[stats::na.omit(s), ]
	
			pt <- pv[up[,2]==i, ]
			at <- av[ua[,2]==i, ]
			pt <- pt[!is.na(s), ]
			at <- at[s, ]

			if (nrow(pt) > 3 & nrow(at) > 3) {
				sb[[i]] <- ssb(npt, nat, ref, lonlat=lonlat, avg=TRUE)
				for (k in 1:lm) {
					em[[k]][[i]] <- evaluate(pt, at, model[[k]])
				}
				en[[i]] <- evaluate(npt, nat, null)
			} 
		}
	}
	
	nauc <- sapply(en, function(x) if (is.null(x)) NA  else x@auc)
	nn <- sapply(en, function(x) if (is.null(x)) NA  else x@np)

	res <- NULL
	for (i in 1:lm) {
		mauc <- sapply(em[[i]], function(x) if (is.null(x)) NA  else x@auc)
		cauc <- mauc - max(0.5, nauc) + 0.5
		res <- cbind(res, cauc)
	}
	colnames(res) <- paste('model', 1:lm, sep='_')
	sb <- sapply(sb, function(x) if (is.null(x)) NA  else mean(x))
	
	return(cbind(n=nn, distance=sb, res))
}

 
