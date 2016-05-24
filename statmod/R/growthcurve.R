meanT <- function(y1,y2) {
#  Mean t-statistic difference between two groups of growth curves
#  Columns are time points, rows are individuals
#  Gordon Smyth
#  14 Feb 2003

	if(is.null(dim(y1)) || is.null(dim(y2))) return(NA)
	y1 <- as.matrix(y1)
	y2 <- as.matrix(y2)
	if(ncol(y1) != ncol(y2)) stop("Number of time points must match")
	m1 <- colMeans(y1,na.rm=TRUE)
	m2 <- colMeans(y2,na.rm=TRUE)
	v1 <- apply(y1,2,var,na.rm=TRUE)
	v2 <- apply(y2,2,var,na.rm=TRUE)
	n1 <- apply(!is.na(y1),2,sum)
	n2 <- apply(!is.na(y2),2,sum)
	s <- ( (n1-1)*v1 + (n2-1)*v2 ) / (n1+n2-2)
	t.stat <- (m1-m2) / sqrt(s*(1/n1+1/n2))
	weighted.mean(t.stat,w=(n1+n2-2)/(n1+n2),na.rm=TRUE)
}

compareTwoGrowthCurves <- function(group,y,nsim=100,fun=meanT) {
#  Permutation test between two groups of growth curves
#  Columns are time points, rows are individuals
#  Gordon Smyth
#  14 Feb 2003

	group <- as.vector(group)
	g <- unique(group)
	if(length(g) != 2) stop("Must be exactly 2 groups")
	stat.obs <- fun(y[group==g[1],,drop=FALSE], y[group==g[2],,drop=FALSE])
	asbig <- 0
	for (i in 1:nsim) {
		pgroup <- sample(group)
		stat <- fun(y[pgroup==g[1],,drop=FALSE], y[pgroup==g[2],,drop=FALSE])
		if(abs(stat) >= abs(stat.obs)) asbig <- asbig+1
	}
	list(stat=stat.obs, p.value=asbig/nsim) 
}

compareGrowthCurves <- function(group,y,levels=NULL,nsim=100,fun=meanT,times=NULL,verbose=TRUE,adjust="holm") {
#  All pairwise permutation tests between groups of growth curves
#  Columns of y are time points, rows are individuals
#  Gordon Smyth
#  14 Feb 2003.  Last modified 17 Nov 2003.

	group <- as.character(group)
	if(is.null(levels)) {
		tab <- table(group)
		tab <- tab[tab >= 2]
		lev <- names(tab)
	} else
		lev <- as.character(levels)
	nlev <- length(lev)
	if(nlev < 2) stop("Less than 2 groups to compare")
	if(is.null(dim(y))) stop("y must be matrix-like")
	y <- as.matrix(y)
	if(!is.null(times)) y <- y[,times,drop=FALSE]

	g1 <- g2 <- rep("",nlev*(nlev-1)/2)
	stat <- pvalue <- rep(0,nlev*(nlev-1)/2)
	pair <- 0
	for (i in 1:(nlev-1)) {
		for (j in (i+1):nlev) {
			if(verbose) cat(lev[i],lev[j])
			pair <- pair+1
			sel <- group %in% c(lev[i],lev[j])
			out <- compareTwoGrowthCurves(group[sel],y[sel,,drop=FALSE],nsim=nsim,fun=fun)
			if(verbose) cat("\ ",round(out$stat,2),"\n")
			g1[pair] <- lev[i]
			g2[pair] <- lev[j]
			stat[pair] <- out$stat
			pvalue[pair] <- out$p.value
		}
	}
	tab <- data.frame(Group1=g1,Group2=g2,Stat=stat,P.Value=pvalue)
	tab$adj.P.Value <- p.adjust(pvalue,method=adjust)
	tab
}

plotGrowthCurves <- function(group,y,levels=sort(unique(group)),times=NULL,col=NULL,...) {
#  Plot growth curves with colors for groups
#  Columns of y are time points, rows are individuals
#  Gordon Smyth
#  30 May 2006.  Last modified 8 July 2006.

	group <- as.character(group)
	if(!is.null(levels)) levels <- as.character(levels)
	nlev <- length(levels)
	if(nlev < 2) stop("Less than 2 groups to compare")
	if(is.null(dim(y))) stop("y must be matrix-like")
	y <- as.matrix(y)
	if(!is.null(times)) y <- y[,times,drop=FALSE]

	if(is.null(col)) col <- 1:nlev
	group.col <- col[match(group,levels)]
	plot(col(y),y,type="n",xlab="Time",ylab="Response",...)
	x <- 1:ncol(y)
	for (i in 1:nrow(y)) {
		lines(x,y[i,],col=group.col[i])
	}
	yr <- range(y,na.rm=TRUE)
	legend(1,yr[2]-diff(yr)/40,legend=levels,col=col,lty=1)
	invisible()
}
