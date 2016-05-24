## Just for memory sake...
##samr.const.quantitative.response <- 0
##samr.const.twoclass.unpaired.response <- 1
##samr.const.survival.response <- 2
##samr.const.multiclass.response <- 3
##samr.const.oneclass.response <- 4
##samr.const.twoclass.paired.response <- 5
##samr.const.twoclass.unpaired.timecourse.response <- 6
##samr.const.oneclass.timecourse.response <- 7
##samr.const.twoclass.paired.timecourse.response <- 8


##SAM R associated functions
## individual functions for each response type

ttest.func <- function(x, y, s0 = 0, sd = NULL) {
	n1 <- sum(y == 1)
	n2 <- sum(y == 2)
	p <- nrow(x)
	m1 <- rowMeans(x[, y == 1, drop = F])
	m2 <- rowMeans(x[, y == 2, drop = F])
	if (is.null(sd)) {
		sd <- sqrt(((n2 - 1) * varr(x[, y == 2], meanx = m2) + 
			(n1 - 1) * varr(x[, y == 1], meanx = m1)) * (1/n1 + 
			1/n2)/(n1 + n2 - 2))
	}
	numer <- m2 - m1
	dif.obs <- (numer)/(sd + s0)
	return(list(tt = dif.obs, numer = numer, sd = sd))
}

wilcoxon.func <- function(x, y, s0 = 0) {
	n1 <- sum(y == 1)
	n2 <- sum(y == 2)
	p = nrow(x)
	r2 = rowSums(t(apply(x, 1, rank))[, y == 2, drop = F])
	numer = r2 - (n2/2) * (n2 + 1) - (n1 * n2)/2
	sd = sqrt(n1 * n2 * (n1 + n2 + 1)/12)
	tt = (numer)/(sd + s0)
	return(list(tt = tt, numer = numer, sd = rep(sd, p)))
}

onesample.ttest.func <- function(x, y, s0 = 0, sd = NULL) {
	n <- length(y)
	x <- x * matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
	m <- rowMeans(x)
	if (is.null(sd)) {
		sd <- sqrt(varr(x, meanx = m)/n)
	}
	dif.obs <- m/(sd + s0)
	return(list(tt = dif.obs, numer = m, sd = sd))
}

patterndiscovery.func = function(x, s0 = 0, eigengene.number = 1) {
	a = mysvd(x, n.components = eigengene.number)
	v = a$v[, eigengene.number]
	# here we try to guess the most interpretable orientation
	#   for the eigengene
	om = abs(a$u[, eigengene.number]) > quantile(abs(a$u[, eigengene.number]), 
		0.95)
	if (median(a$u[om, eigengene.number]) < 0) {
		v = -1 * v
	}
	aa = quantitative.func(x, v, s0 = s0)
	eigengene = cbind(1:nrow(a$v), v)
	dimnames(eigengene) = list(NULL, c("sample number", "value"))
	return(list(tt = aa$tt, numer = aa$numer, sd = aa$sd, eigengene = eigengene))
}

paired.ttest.func <- function(x, y, s0 = 0, sd = NULL) {
	nc <- ncol(x)/2
	o <- 1:nc
	o1 <- rep(0, ncol(x)/2)
	o2 <- o1
	for (j in 1:nc) {
		o1[j] <- (1:ncol(x))[y == -o[j]]
	}
	for (j in 1:nc) {
		o2[j] <- (1:ncol(x))[y == o[j]]
	}
	d <- x[, o2, drop = F] - x[, o1, drop = F]
	su <- x[, o2, drop = F] + x[, o1, drop = F]
	if (is.matrix(d)) {
		m <- rowMeans(d)
	}
	if (!is.matrix(d)) {
		m <- mean(d)
	}
	if (is.null(sd)) {
		if (is.matrix(d)) {
			sd <- sqrt(varr(d, meanx = m)/nc)
		}
		if (!is.matrix(d)) {
			sd <- sqrt(var(d)/nc)
		}
	}
	dif.obs <- m/(sd + s0)
	return(list(tt = dif.obs, numer = m, sd = sd))
}

cox.func <- function(x, y, censoring.status, s0 = 0) {
	# find the index matrix
	Dn <- sum(censoring.status == 1)
	Dset <- c(1:ncol(x))[censoring.status == 1]  # the set of observed
	ind <- matrix(0, ncol(x), Dn)
	# get the matrix
	for (i in 1:Dn) {
		ind[y > y[Dset[i]] - 1e-08, i] <- 1/sum(y > y[Dset[i]] - 
			1e-08)
	}
	ind.sums <- rowSums(ind)
	x.ind <- x %*% ind
	# get the derivatives
	numer <- x %*% (censoring.status - ind.sums)
	sd <- sqrt((x * x) %*% ind.sums - rowSums(x.ind * x.ind))
	tt <- numer/(sd + s0)
	return(list(tt = tt, numer = numer, sd = sd))
}

multiclass.func <- function(x, y, s0 = 0) {
	##assumes y is coded 1,2...
	nn <- table(y)
	m <- matrix(0, nrow = nrow(x), ncol = length(nn))
	v <- m
	for (j in 1:length(nn)) {
		m[, j] <- rowMeans(x[, y == j])
		v[, j] <- (nn[j] - 1) * varr(x[, y == j], meanx = m[, 
			j])
	}
	mbar <- rowMeans(x)
	mm <- m - matrix(mbar, nrow = length(mbar), ncol = length(nn))
	fac <- (sum(nn)/prod(nn))
	scor <- sqrt(fac * (apply(matrix(nn, nrow = nrow(m), ncol = ncol(m), 
		byrow = TRUE) * mm * mm, 1, sum)))
	sd <- sqrt(rowSums(v) * (1/sum(nn - 1)) * sum(1/nn))
	tt <- scor/(sd + s0)
	mm.stand = t(scale(t(mm), center = FALSE, scale = sd))
	return(list(tt = tt, numer = scor, sd = sd, stand.contrasts = mm.stand))
}

#quantitative.func <- function(x,y,s0=0){
#  yy <- y-mean(y)
#  temp <- x%*%yy
#mx=rowMeans(x)
#sxx <-rowSums( (x-mx%*%t(rep(1,ncol(x))))^2 )
#
#  scor <- temp/sxx
#  b0hat <- mean(y)-scor*mx
# yhat <-
#   matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+x*matrix(scor,nrow=nrow(x),ncol=ncol(x))
# ty <-
#   matrix(y,nrow=nrow(yhat),ncol=ncol(yhat),byrow=TRUE)
#  sigma <- sqrt(rowSums((ty-yhat)^2)/(ncol(yhat)-2))
#  sd <- sigma/sqrt(sxx)
#  tt <- scor/(sd+s0)
#  return(list(tt=tt, numer=scor, sd=sd))
#
#}

quantitative.func <- function(x, y, s0 = 0) {
	# regression of x on y
	my = mean(y)
	yy <- y - my
	temp <- x %*% yy
	mx = rowMeans(x)
	syy = sum(yy^2)
	scor <- temp/syy
	b0hat <- mx - scor * my
	ym = matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = T)
	xhat <- matrix(b0hat, nrow = nrow(x), ncol = ncol(x)) + ym * 
		matrix(scor, nrow = nrow(x), ncol = ncol(x))
	sigma <- sqrt(rowSums((x - xhat)^2)/(ncol(xhat) - 2))
	sd <- sigma/sqrt(syy)
	tt <- scor/(sd + s0)
	return(list(tt = tt, numer = scor, sd = sd))
}

timearea.func <- function(x, y, s0 = 0) {
	n <- ncol(x)
	xx <- 0.5 * (x[, 2:n] + x[, 1:(n - 1)]) * matrix(diff(y), 
		nrow = nrow(x), ncol = n - 1, byrow = T)
	numer <- rowMeans(xx)
	sd <- sqrt(varr(xx, meanx = numer)/n)
	tt <- numer/sqrt(sd + s0)
	return(list(tt = tt, numer = numer, sd = sd))
}

#########################################
detec.slab <- function(samr.obj, del, min.foldchange) {
	## find genes above and below the slab of half-width del
	# this calculation is tricky- for consistency, the slab
	#   condition picks
	# all genes that are beyond the first departure from the
	#   slab
	# then the fold change condition is applied (if applicable)
	n <- length(samr.obj$tt)
	tt <- samr.obj$tt
	evo <- samr.obj$evo
	numer <- samr.obj$tt * (samr.obj$sd + samr.obj$s0)
	tag <- order(tt)
	pup <- NULL
	foldchange.cond.up = rep(T, length(evo))
	foldchange.cond.lo = rep(T, length(evo))
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		foldchange.cond.up = samr.obj$foldchange >= min.foldchange
		foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
	}
	o1 <- (1:n)[(tt[tag] - evo > del) & evo > 0]
	if (length(o1) > 0) {
		o1 <- o1[1]
		o11 <- o1:n
		o111 <- rep(F, n)
		o111[tag][o11] <- T
		pup <- (1:n)[o111 & foldchange.cond.up]
	}
	plow <- NULL
	o2 <- (1:n)[(evo - tt[tag] > del) & evo < 0]
	if (length(o2) > 0) {
		o2 <- o2[length(o2)]
		o22 <- 1:o2
		o222 <- rep(F, n)
		o222[tag][o22] <- T
		plow <- (1:n)[o222 & foldchange.cond.lo]
	}
	return(list(plow = plow, pup = pup))
}

sumlengths <- function(aa) {
	length(aa$pl) + length(aa$pu)
}

## Jun added starts
samr.compute.delta.table <- function(samr.obj, min.foldchange = 0, 
	dels = NULL, nvals = 50) {
	res <- NULL
	if (samr.obj$assay.type == "array") {
		res <- samr.compute.delta.table.array(samr.obj, min.foldchange, 
			dels, nvals)
	}
	else if (samr.obj$assay.type == "seq") {
		res <- samr.compute.delta.table.seq(samr.obj, min.foldchange, 
			dels)
	}
	return(res)
}
## Jun added ends

## Jun added the first row below, and commented the row
#   after it
samr.compute.delta.table.array <- function(samr.obj, 
	min.foldchange = 0, dels = NULL, nvals = 50) {
	#samr.compute.delta.table <- function(samr.obj,
	#   min.foldchange=0, dels=NULL, nvals=50) {
	# computes delta table, starting with samr object 'a', for
	#   nvals values of delta
	lmax = sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
	if (is.null(dels)) {
		dels = (seq(0, lmax, length = nvals)^2)
	}
	col = matrix(1, nrow = length(samr.obj$evo), ncol = nvals)
	ttstar0 <- samr.obj$ttstar0
	tt <- samr.obj$tt
	n <- samr.obj$n
	evo <- samr.obj$evo
	nsim <- ncol(ttstar0)
	res1 <- NULL
	foldchange.cond.up = matrix(T, nrow = nrow(samr.obj$ttstar), 
		ncol = ncol(samr.obj$ttstar))
	foldchange.cond.lo = matrix(T, nrow = nrow(samr.obj$ttstar), 
		ncol = ncol(samr.obj$ttstar))
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		foldchange.cond.up = samr.obj$foldchange.star >= min.foldchange
		foldchange.cond.lo = samr.obj$foldchange.star <= 1/min.foldchange
	}
	cutup = rep(NA, length(dels))
	cutlow = rep(NA, length(dels))
	g2 = rep(NA, length(dels))
	errup = matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
	errlow = matrix(NA, ncol = length(dels), nrow = ncol(samr.obj$ttstar0))
	cat("", fill = T)
	cat("Computing delta table", fill = T)
	for (ii in 1:length(dels)) {
		cat(ii, fill = TRUE)
		ttt <- detec.slab(samr.obj, dels[ii], min.foldchange)
		cutup[ii] <- 1e+10
		if (length(ttt$pup > 0)) {
			cutup[ii] <- min(samr.obj$tt[ttt$pup])
		}
		cutlow[ii] <- -1e+10
		if (length(ttt$plow) > 0) {
			cutlow[ii] <- max(samr.obj$tt[ttt$plow])
		}
		g2[ii] = sumlengths(ttt)
		errup[, ii] = colSums(samr.obj$ttstar0 > cutup[ii] & 
			foldchange.cond.up)
		errlow[, ii] = colSums(samr.obj$ttstar0 < cutlow[ii] & 
			foldchange.cond.lo)
	}
	s <- sqrt(apply(errup, 2, var)/nsim + apply(errlow, 2, var)/nsim)
	gmed <- apply(errup + errlow, 2, median)
	g90 = apply(errup + errlow, 2, quantile, 0.9)
	res1 <- cbind(samr.obj$pi0 * gmed, samr.obj$pi0 * g90, g2, 
		samr.obj$pi0 * gmed/g2, samr.obj$pi0 * g90/g2, cutlow, 
		cutup)
	res1 <- cbind(dels, res1)
	# remove rows with #called=0
	#om=res1[,4]==0
	#res1=res1[!om,,drop=F]
	# remove duplicate rows with same # of genes called
	#omm=!duplicated(res1[,4])
	#res1=res1[omm,,drop=F]
	dimnames(res1) <- list(NULL, c("delta", "# med false pos", 
		"90th perc false pos", "# called", "median FDR", "90th perc FDR", 
		"cutlo", "cuthi"))
	return(res1)
}

detec.horiz <- function(samr.obj, cutlow, cutup, min.foldchange) {
	## find genes above or below horizontal cutpoints
	dobs <- samr.obj$tt
	n <- length(dobs)
	foldchange.cond.up = rep(T, n)
	foldchange.cond.lo = rep(T, n)
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		foldchange.cond.up = samr.obj$foldchange >= min.foldchange
		foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
	}
	pup <- (1:n)[dobs > cutup & foldchange.cond.up]
	plow <- (1:n)[dobs < cutlow & foldchange.cond.lo]
	return(list(plow = plow, pup = pup))
}

samr.plot <- function(samr.obj, del = NULL, min.foldchange = 0) {
	## make observed-expected plot
	## takes foldchange into account too
	if (is.null(del)) {
		del = sqrt(max(abs(sort(samr.obj$tt) - samr.obj$evo)))
	}
	LARGE = 1e+10
	b <- detec.slab(samr.obj, del, min.foldchange)
	bb <- c(b$pup, b$plow)
	b1 = LARGE
	b0 = -LARGE
	if (!is.null(b$pup)) {
		b1 <- min(samr.obj$tt[b$pup])
	}
	if (!is.null(b$plow)) {
		b0 <- max(samr.obj$tt[b$plow])
	}
	c1 <- (1:samr.obj$n)[sort(samr.obj$tt) >= b1]
	c0 <- (1:samr.obj$n)[sort(samr.obj$tt) <= b0]
	c2 <- c(c0, c1)
	foldchange.cond.up = rep(T, length(samr.obj$evo))
	foldchange.cond.lo = rep(T, length(samr.obj$evo))
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		foldchange.cond.up = samr.obj$foldchange >= min.foldchange
		foldchange.cond.lo = samr.obj$foldchange <= 1/min.foldchange
	}
	col = rep(1, length(samr.obj$evo))
	col[b$plow] = 3
	col[b$pup] = 2
	if (!is.null(samr.obj$foldchange[1]) & (min.foldchange > 
		0)) {
		col[!foldchange.cond.lo & !foldchange.cond.up] = 1
	}
	col.ordered = col[order(samr.obj$tt)]
	ylims <- range(samr.obj$tt)
	xlims <- range(samr.obj$evo)
	plot(samr.obj$evo, sort(samr.obj$tt), xlab = "expected score", 
		ylab = "observed score", ylim = ylims, xlim = xlims, 
		type = "n")
	points(samr.obj$evo, sort(samr.obj$tt), col = col.ordered)
	abline(0, 1)
	abline(del, 1, lty = 2)
	abline(-del, 1, lty = 2)
}

localfdr <- function(samr.obj, min.foldchange, perc = 0.01, 
	df = 10) {
	## estimates compute.localfdr at score 'd', using SAM
	#   object 'samr.obj'
	## 'd' can be a vector of d scores
	## returns estimate of symmetric fdr  as a percentage
	# this version uses a 1% symmetric window, and does not
	#   estimate fdr in
	# windows  having fewer than 100 genes
	## to use: first run samr and then pass the resulting fit
	#   object to
	## localfdr
	## NOTE: at most 20 of the perms are used to estimate the
	#   fdr (for speed sake)
	# I try two window shapes: symmetric and an assymetric one
	# currently I use the symmetric window to estimate the
	#   compute.localfdr
	ngenes = length(samr.obj$tt)
	mingenes = 50
	# perc is increased, in order to get at least mingenes in a
	#   window
	perc = max(perc, mingenes/length(samr.obj$tt))
	nperms.to.use = min(20, ncol(samr.obj$ttstar))
	nperms = ncol(samr.obj$ttstar)
	d = seq(sort(samr.obj$tt)[1], sort(samr.obj$tt)[ngenes], 
		length = 100)
	ndscore <- length(d)
	dvector <- rep(NA, ndscore)
	ind.foldchange = rep(T, length(samr.obj$tt))
	if (!is.null(samr.obj$foldchange[1]) & min.foldchange > 0) {
		ind.foldchange = (samr.obj$foldchange >= min.foldchange) | 
			(samr.obj$foldchange <= min.foldchange)
	}
	fdr.temp = function(temp, dlow, dup, pi0, ind.foldchange) {
		return(sum(pi0 * (temp >= dlow & temp <= dup & ind.foldchange)))
	}
	for (i in 1:ndscore) {
		pi0 <- samr.obj$pi0
		r <- sum(samr.obj$tt < d[i])
		r22 <- round(max(r - length(samr.obj$tt) * perc/2, 1))
		dlow.sym <- sort(samr.obj$tt)[r22]
		#      if(d[i]<0)
		#       {
		#         r2 <- max(r-length(samr.obj$tt)*perc/2, 1)
		# r22= min(r+length(samr.obj$tt)*perc/2,
		#   length(samr.obj$tt))
		#
		#          dlow <- sort(samr.obj$tt)[r2]
		#          dup=sort(samr.obj$tt)[r22]
		#       }
		r22 <- min(r + length(samr.obj$tt) * perc/2, length(samr.obj$tt))
		dup.sym <- sort(samr.obj$tt)[r22]
		#     if(d[i]>0)
		#      {
		# r2 <- min(r+length(samr.obj$tt)*perc/2,
		#   length(samr.obj$tt))
		#        r22 <- max(r-length(samr.obj$tt)*perc/2, 1)
		#        dup <- sort(samr.obj$tt)[r2]
		#        dlow <- sort(samr.obj$tt)[r22]
		#
		#       }
		# o <- samr.obj$tt>=dlow & samr.obj$tt<= dup &
		#   ind.foldchange
		oo <- samr.obj$tt >= dlow.sym & samr.obj$tt <= dup.sym & 
			ind.foldchange
		nsim <- ncol(samr.obj$ttstar)
		fdr <- rep(NA, nsim)
		fdr2 <- fdr
		if (!is.null(samr.obj$foldchange[1]) & min.foldchange > 
			0) {
			temp = as.vector(samr.obj$foldchange.star[, 1:nperms.to.use])
			ind.foldchange = (temp >= min.foldchange) | (temp <= 
				min.foldchange)
		}
		temp = samr.obj$ttstar0[, sample(1:nperms, size = nperms.to.use)]
		# fdr <-median(apply(temp,2,fdr.temp,dlow, dup, pi0,
		#   ind.foldchange))
		fdr.sym <- median(apply(temp, 2, fdr.temp, dlow.sym, 
			dup.sym, pi0, ind.foldchange))
		#      fdr <- 100*fdr/sum(o)
		fdr.sym <- 100 * fdr.sym/sum(oo)
		dlow.sym <- dlow.sym
		dup.sym <- dup.sym
		dvector[i] <- fdr.sym
	}
	om = !is.na(dvector) & (dvector != Inf)
	aa = smooth.spline(d[om], dvector[om], df = df)
	return(list(smooth.object = aa, perc = perc, df = df))
}

predictlocalfdr = function(smooth.object, d) {
	yhat = predict(smooth.object, d)$y
	yhat = pmin(yhat, 100)
	yhat = pmax(yhat, 0)
	return(yhat)
}

samr.compute.siggenes.table = function(samr.obj, del, 
	data, delta.table, min.foldchange = 0, all.genes = FALSE, 
	compute.localfdr = FALSE)
{
	## computes significant genes table, starting with samr
	#   object 'a' and 'delta.table'
	##  for a  **single** value del
	## if all.genes is true, all genes are printed (and value
	#   of del is ignored)
	if (is.null(data$geneid))
	{
		data$geneid = paste("g", 1:nrow(data$x), sep = "")
	}
	if (is.null(data$genenames))
	{
		data$genenames = paste("g", 1:nrow(data$x), sep = "")
	}
	if (!all.genes)
	{
		sig = detec.slab(samr.obj, del, min.foldchange)
	}
	if (all.genes)
	{
		p = length(samr.obj$tt)
		pup = (1:p)[samr.obj$tt >= 0]
		plo = (1:p)[samr.obj$tt < 0]
		sig = list(pup = pup, plo = plo)
	}
	if (compute.localfdr)
	{
		aa = localfdr(samr.obj, min.foldchange)
		if (length(sig$pup) > 0)
		{
			fdr.up = predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$pup])
		}
		if (length(sig$plo) > 0)
		{
			fdr.lo = predictlocalfdr(aa$smooth.object, samr.obj$tt[sig$plo])
		}
	}
	qvalues = NULL
	if (length(sig$pup) > 0 | length(sig$plo) > 0)
	{
		qvalues = qvalue.func(samr.obj, sig, delta.table)
	}
	res.up = NULL
	res.lo = NULL
	done = FALSE
	
	# two class unpaired or paired  (foldchange is reported)
	if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
		samr.obj$resp.type == samr.const.twoclass.paired.response))
	{
		if (!is.null(sig$pup))
		{
			res.up = cbind(sig$pup + 1, data$genenames[sig$pup], 
				data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup], 
				samr.obj$sd[sig$pup], samr.obj$foldchange[sig$pup], 
				qvalues$qvalue.up)
			if (compute.localfdr)
			{
				res.up = cbind(res.up, fdr.up)
			}
			temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
				"Score(d)", "Numerator(r)", "Denominator(s+s0)", 
				"Fold Change", "q-value(%)"))
			if (compute.localfdr)
			{
				temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
			}
			dimnames(res.up) = temp.names
		}
		if (!is.null(sig$plo))
		{
			res.lo = cbind(sig$plo + 1, data$genenames[sig$plo], 
				data$geneid[sig$plo], samr.obj$tt[sig$plo], samr.obj$numer[sig$plo], 
				samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo], 
				qvalues$qvalue.lo)
			if (compute.localfdr)
			{
				res.lo = cbind(res.lo, fdr.lo)
			}
			temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
				"Score(d)", "Numerator(r)", "Denominator(s+s0)", 
				"Fold Change", "q-value(%)"))
			if (compute.localfdr)
			{
				temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
			}
			dimnames(res.lo) = temp.names
		}
		done = TRUE
	}
	
	# multiclass
	if (samr.obj$resp.type == samr.const.multiclass.response)
	{
		if (!is.null(sig$pup))
		{
			res.up = cbind(sig$pup + 1, data$genenames[sig$pup], 
			data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup], 
			samr.obj$sd[sig$pup], samr.obj$stand.contrasts[sig$pup, ], qvalues$qvalue.up)
	
			if (compute.localfdr)
			{
				res.up = cbind(res.up, fdr.up)
			}
			
			collabs.contrast = paste("contrast-", as.character(1:ncol(samr.obj$stand.contrasts)), 
				sep = "")
			temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
			"Score(d)", "Numerator(r)", "Denominator(s+s0)", 
			collabs.contrast, "q-value(%)"))
			
			if (compute.localfdr)
			{
				temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
			}
			dimnames(res.up) = temp.names
		}
		res.lo = NULL
		done = TRUE
	}
	
	#all other cases
	if (!done)
	{
		if (!is.null(sig$pup))
		{
			res.up = cbind(sig$pup + 1, data$genenames[sig$pup], 
				data$geneid[sig$pup], samr.obj$tt[sig$pup], samr.obj$numer[sig$pup], 
				samr.obj$sd[sig$pup], samr.obj$foldchange[sig$pup], 
				qvalues$qvalue.up)
			if (compute.localfdr)
			{
				res.up = cbind(res.up, fdr.up)
			}
			temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
				"Score(d)", "Numerator(r)", "Denominator(s+s0)", 
				"q-value(%)"))
			if (compute.localfdr)
			{
				temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
			}
			dimnames(res.up) = temp.names
		}
		if (!is.null(sig$plo))
		{
			res.lo = cbind(sig$plo + 1, data$genenames[sig$plo], 
				data$geneid[sig$plo], samr.obj$tt[sig$plo], samr.obj$numer[sig$plo], 
				samr.obj$sd[sig$plo], samr.obj$foldchange[sig$plo], 
				qvalues$qvalue.lo)
			if (compute.localfdr)
			{
				res.lo = cbind(res.lo, fdr.lo)
			}
			temp.names = list(NULL, c("Row", "Gene ID", "Gene Name", 
				"Score(d)", "Numerator(r)", "Denominator(s+s0)", 
				"q-value(%)"))
			if (compute.localfdr)
			{
				temp.names[[2]] = c(temp.names[[2]], "localfdr(%)")
			}
			dimnames(res.lo) = temp.names
		}
		done = TRUE
	}
	if (!is.null(res.up))
	{
		o1 = order(-samr.obj$tt[sig$pup])
		res.up = res.up[o1, , drop = F]
	}
	if (!is.null(res.lo))
	{
		o2 = order(samr.obj$tt[sig$plo])
		res.lo = res.lo[o2, , drop = F]
	}
	color.ind.for.multi = NULL
	if (samr.obj$resp.type == samr.const.multiclass.response & !is.null(sig$pup))
	{
		color.ind.for.multi = 1 * (samr.obj$stand.contrasts[sig$pup, 
			] > samr.obj$stand.contrasts.95[2]) + (-1) * (samr.obj$stand.contrasts[sig$pup, 
			] < samr.obj$stand.contrasts.95[1])
	}
	ngenes.up = nrow(res.up)
	if (is.null(ngenes.up))
	{
		ngenes.up = 0
	}
	ngenes.lo = nrow(res.lo)
	if (is.null(ngenes.lo))
	{
		ngenes.lo = 0
	}
	return(list(genes.up = res.up, genes.lo = res.lo, color.ind.for.multi = color.ind.for.multi, 
		ngenes.up = ngenes.up, ngenes.lo = ngenes.lo))
}

qvalue.func = function(samr.obj, sig, delta.table) {
	# returns q-value as a percentage (out of 100)
	LARGE = 1e+10
	qvalue.up = rep(NA, length(sig$pup))
	o1 = sig$pup
	cutup = delta.table[, 8]
	FDR = delta.table[, 5]
	ii = 0
	for (i in o1) {
		o = abs(cutup - samr.obj$tt[i])
		o[is.na(o)] = LARGE
		oo = (1:length(o))[o == min(o)]
		oo = oo[length(oo)]
		ii = ii + 1
		qvalue.up[ii] = FDR[oo]
	}
	qvalue.lo = rep(NA, length(sig$plo))
	o2 = sig$plo
	cutlo = delta.table[, 7]
	ii = 0
	for (i in o2) {
		o = abs(cutlo - samr.obj$tt[i])
		o[is.na(o)] = LARGE
		oo = (1:length(o))[o == min(o)]
		oo = oo[length(oo)]
		ii = ii + 1
		qvalue.lo[ii] = FDR[oo]
	}
	# any qvalues that are missing, are set to 1 (the highest
	#   value)
	qvalue.lo[is.na(qvalue.lo)] = 1
	qvalue.up[is.na(qvalue.up)] = 1
	# ensure that each qvalue vector is monotone non-increasing
	o1 = order(samr.obj$tt[sig$plo])
	qv1 = qvalue.lo[o1]
	qv11 = qv1
	if (length(qv1) > 1) {
		for (i in 2:length(qv1)) {
			if (qv11[i] < qv11[i - 1]) {
				qv11[i] = qv11[i - 1]
			}
		}
		qv111 = qv11
		qv111[o1] = qv11
	}
	else {
		qv111 = qv1
	}
	o2 = order(samr.obj$tt[sig$pup])
	qv2 = qvalue.up[o2]
	qv22 = qv2
	if (length(qv2) > 1) {
		for (i in 2:length(qv2)) {
			if (qv22[i] > qv22[i - 1]) {
				qv22[i] = qv22[i - 1]
			}
		}
		qv222 = qv22
		qv222[o2] = qv22
	}
	else {
		qv222 = qv2
	}
	return(list(qvalue.lo = 100 * qv111, qvalue.up = 100 * qv222))
}

foldchange.twoclass = function(x, y, logged2) {
	#  if(logged2){x=2^x}
	m1 <- rowMeans(x[, y == 1, drop = F])
	m2 <- rowMeans(x[, y == 2, drop = F])
	if (!logged2) {
		fc = m2/m1
	}
	if (logged2) {
		fc = 2^{
			m2 - m1
		}
	}
	return(fc)
}

foldchange.paired = function(x, y, logged2) {
	#  if(logged2){x=2^x}
	nc <- ncol(x)/2
	o <- 1:nc
	o1 <- rep(0, ncol(x)/2)
	o2 <- o1
	for (j in 1:nc) {
		o1[j] <- (1:ncol(x))[y == -o[j]]
	}
	for (j in 1:nc) {
		o2[j] <- (1:ncol(x))[y == o[j]]
	}
	if (!logged2) {
		d <- x[, o2, drop = F]/x[, o1, drop = F]
	}
	if (logged2) {
		d <- x[, o2, drop = F] - x[, o1, drop = F]
	}
	if (!logged2) {
		fc <- rowMeans(d)
	}
	if (logged2) {
		fc <- 2^rowMeans(d)
	}
	return(fc)
}

est.s0 <- function(tt, sd, s0.perc = seq(0, 1, by = 0.05)) {
	## estimate s0 (exchangeability) factor for denominator.
	## returns the actual estimate s0 (not a percentile)
	br = unique(quantile(sd, seq(0, 1, len = 101)))
	nbr = length(br)
	a <- cut(sd, br, labels = F)
	a[is.na(a)] <- 1
	cv.sd <- rep(0, length(s0.perc))
	for (j in 1:length(s0.perc)) {
		w <- quantile(sd, s0.perc[j])
		w[j == 1] <- 0
		tt2 <- tt * sd/(sd + w)
		tt2[tt2 == Inf] = NA
		sds <- rep(0, nbr - 1)
		for (i in 1:(nbr - 1)) {
			sds[i] <- mad(tt2[a == i], na.rm = TRUE)
		}
		cv.sd[j] <- sqrt(var(sds))/mean(sds)
	}
	o = (1:length(s0.perc))[cv.sd == min(cv.sd)]
	# we don;t allow taking s0.hat to be 0th percentile when
	#   min sd is 0
	s0.hat = quantile(sd[sd != 0], s0.perc[o])
	return(list(s0.perc = s0.perc, cv.sd = cv.sd, s0.hat = s0.hat))
}

samr.missrate <- function(samr.obj, del, delta.table, 
	quant = NULL) {
	# returns miss rate as a percentage
	if (is.null(quant)) {
		if (samr.obj$resp.type != samr.const.multiclass.response) {
			quant = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.75, 0.8, 
				0.85, 0.9, 0.95, 1)
		}
		if (samr.obj$resp.type == samr.const.multiclass.response) {
			quant = c(0.75, 0.8, 0.85, 0.9, 0.95, 1)
		}
	}
	## estimate miss rate from sam object 'a'
	o = abs(delta.table[, 1] - del)
	oo = (1:nrow(delta.table))[o == min(o)]
	cut.lo = delta.table[oo, 7]
	cut.up = delta.table[oo, 8]
	ooo = samr.obj$tt > cut.lo & samr.obj$tt < cut.up
	cuts = quantile(samr.obj$tt[ooo], quant)
	ncuts <- length(cuts)
	ngenes <- rep(NA, ncuts)
	ngenes0 <- rep(NA, ncuts)
	ngenes2 <- rep(NA, ncuts)
	missrate <- rep(NA, ncuts)
	nperm = ncol(samr.obj$ttstar)
	for (j in 1:(ncuts - 1)) {
		ngenes2[j] <- sum(samr.obj$tt > cuts[j] & samr.obj$tt < 
			cuts[j + 1])
		ngenes0[j] <- sum(samr.obj$ttstar > cuts[j] & samr.obj$ttstar < 
			cuts[j + 1])/nperm
		missrate[j] <- (ngenes2[j] - samr.obj$pi0 * ngenes0[j])/ngenes2[j]
		missrate[j] <- max(missrate[j], 0)
	}
	cuts = round(cuts, 3)
	res = matrix(NA, ncol = 3, nrow = ncuts - 1)
	missrate = round(missrate, 4)
	for (i in 1:(ncuts - 1)) {
		res[i, 1] = paste(as.character(quant[i]), as.character(quant[i + 
			1]), sep = " -> ")
		res[i, 2] = paste(as.character(cuts[i]), as.character(cuts[i + 
			1]), sep = " -> ")
		res[i, 3] = 100 * missrate[i]
	}
	dimnames(res) = list(NULL, c("Quantiles", "Cutpoints", "Miss Rate(%)"))
	return(res)
}

varr <- function(x, meanx = NULL) {
	n <- ncol(x)
	p <- nrow(x)
	Y <- matrix(1, nrow = n, ncol = 1)
	if (is.null(meanx)) {
		meanx <- rowMeans(x)
	}
	ans <- rep(1, p)
	xdif <- x - meanx %*% t(Y)
	ans <- (xdif^2) %*% rep(1/(n - 1), n)
	ans <- drop(ans)
	return(ans)
}

samr.options <- list(debug=TRUE, #whether to turn on debugging or not
	err.file=ifelse(.Platform$OS.type=='windows', 'C:/samrtrace.txt', 'samrtrace.txt'),
	image.file=ifelse(.Platform$OS.type=='windows', 'C:/samrimage.Rdata', 'samrimage.Rdata'))

#
# Our error handler
#

.error.trace <- function() {
	err.message <- geterrmessage()
	if (!is.null(samr.options$image.file)) {
		save.image(samr.options$image.file)
	}
	if (!is.null(samr.options$err.file)) {
		sink(samr.options$err.file)
		print(err.message)
		traceback()
		sink()
	}
	winDialog(type = "ok", message = err.message)
}

##
## Upon loading, if we are in a windows environment, we use
#   the windows
## dialog mechanism to display errors. Useful for debugging
#   COM apps
##
.onLoad <- function(lib, pkg) {
	if (.Platform$OS.type == "windows") {
		# options(error=function() winDialog(type='ok',
		#   message=geterrmessage()))
		options(error = samr.xl.error.trace)
	}
}

##
## Upon unload, we set things back the way they were...
##
.onUnload <- function(libpath) {
	if (.Platform$OS.type == "windows") {
		options(error = NULL)
	}
}

samr.xl.build.data <- function(x, y, geneid, genenames, 
	logged2) {
	return(list(x = x, y = y, geneid = geneid, genenames = genenames, 
		logged2 = logged2))
}

insert.value <- function(vec, newval, pos) {
	if (pos == 1) 
		return(c(newval, vec))
	lvec <- length(vec)
	if (pos > lvec) 
		return(c(vec, newval))
	return(c(vec[1:pos - 1], newval, vec[pos:lvec]))
}

permute <- function(elem) {
	# generates all perms of the vector elem
	if (!missing(elem)) {
		if (length(elem) == 2) 
			return(matrix(c(elem, elem[2], elem[1]), nrow = 2))
		last.matrix <- permute(elem[-1])
		dim.last <- dim(last.matrix)
		new.matrix <- matrix(0, nrow = dim.last[1] * (dim.last[2] + 
			1), ncol = dim.last[2] + 1)
		for (row in 1:(dim.last[1])) {
			for (col in 1:(dim.last[2] + 1)) new.matrix[row + 
				(col - 1) * dim.last[1], ] <- insert.value(last.matrix[row, 
				], elem[1], col)
		}
		return(new.matrix)
	}
	else cat("Usage: permute(elem)\n\twhere elem is a vector\n")
}

sample.perms <- function(elem, nperms) {
	# randomly generates  nperms of the vector elem
	res = permute.rows(matrix(elem, nrow = nperms, ncol = length(elem), 
		byrow = T))
	return(res)
}

integer.base.b <- function(x, b = 2) {
	xi <- as.integer(x)
	if (xi == 0) {
		return(0)
	}
	if (any(is.na(xi) | ((x - xi) != 0))) 
		print(list(ERROR = "x not integer", x = x))
	N <- length(x)
	xMax <- max(x)
	ndigits <- (floor(logb(xMax, base = 2)) + 1)
	Base.b <- array(NA, dim = c(N, ndigits))
	for (i in 1:ndigits) {
		#i <- 1
		Base.b[, ndigits - i + 1] <- (x%%b)
		x <- (x%/%b)
	}
	if (N == 1) 
		Base.b[1, ]
	else Base.b
}

compute.block.perms = function(y, blocky, nperms) {
	# y are the data (eg class label 1 vs 2; or -1,1, -2,2 for
	#   paired data)
	# blocky are the block labels (abs(y) for paired daatr)
	ny = length(y)
	nblocks = length(unique(blocky))
	tab = table(blocky)
	total.nperms = prod(factorial(tab))
	# block.perms is a list of all possible permutations
	block.perms = vector("list", nblocks)
	# first enumerate all perms, when possible
	if (total.nperms <= nperms) {
		all.perms.flag = 1
		nperms.act = total.nperms
		for (i in 1:nblocks) {
			block.perms[[i]] = permute(y[blocky == i])
		}
		kk = 0:(factorial(max(tab))^nblocks - 1)
		#the rows of the matrix outerm runs through the 'outer
		#   product'
		# first we assume that all blocks have max(tab) members;
		#   then we remove rows of outerm that
		#  are illegal (ie when a block has fewer members)
		outerm = matrix(0, nrow = length(kk), ncol = nblocks)
		for (i in 1:length(kk)) {
			kkkk = integer.base.b(kk[i], b = factorial(max(tab)))
			if (length(kkkk) > nblocks) {
				kkkk = kkkk[(length(kkkk) - nblocks + 1):length(kkkk)]
			}
			outerm[i, (nblocks - length(kkkk) + 1):nblocks] = kkkk
		}
		outerm = outerm + 1
		# now remove rows that are illegal perms
		ind = rep(TRUE, nrow(outerm))
		for (j in 1:ncol(outerm)) {
			ind = ind & outerm[, j] <= factorial(tab[j])
		}
		outerm = outerm[ind, , drop = F]
		# finally, construct permutation matrix from outer product
		permsy = matrix(NA, nrow = total.nperms, ncol = ny)
		for (i in 1:total.nperms) {
			junk = NULL
			for (j in 1:nblocks) {
				junk = c(junk, block.perms[[j]][outerm[i, j], 
				  ])
			}
			permsy[i, ] = junk
		}
	}
	# next handle case when there are too many perms to
	#   enumerate
	if (total.nperms > nperms) {
		all.perms.flag = 0
		nperms.act = nperms
		permsy = NULL
		block.perms = vector("list", nblocks)
		for (j in 1:nblocks) {
			block.perms[[j]] = sample.perms(y[blocky == j], nperms = nperms)
		}
		for (j in 1:nblocks) {
			permsy = cbind(permsy, block.perms[[j]])
		}
	}
	return(list(permsy = permsy, all.perms.flag = all.perms.flag, 
		nperms.act = nperms.act))
}

getperms = function(y, nperms) {
	total.perms = factorial(length(y))
	if (total.perms <= nperms) {
		perms = permute(1:length(y))
		all.perms.flag = 1
		nperms.act = total.perms
	}
	if (total.perms > nperms) {
		perms = matrix(NA, nrow = nperms, ncol = length(y))
		for (i in 1:nperms) {
			perms[i, ] = sample(1:length(y), size = length(y))
		}
		all.perms.flag = 0
		nperms.act = nperms
	}
	return(list(perms = perms, all.perms.flag = all.perms.flag, 
		nperms.act = nperms.act))
}

parse.block.labels.for.2classes = function(y) {
	#this only works for 2 class case- having form jBlockn,
	#   where j=1 or 2
	n = length(y)
	y.act = rep(NA, n)
	blocky = rep(NA, n)
	for (i in 1:n) {
		blocky[i] = as.numeric(substring(y[i], 7, nchar(y[i])))
		y.act[i] = as.numeric(substring(y[i], 1, 1))
	}
	return(list(y.act = y.act, blocky = blocky))
}

parse.time.labels.and.summarize.data = function(x, 
	y, resp.type, time.summary.type) {
	# parse time labels, and summarize time data for each
	#   person, via a slope or area
	# does some error checking too
	n = length(y)
	last5char = rep(NA, n)
	last3char = rep(NA, n)
	for (i in 1:n) {
		last3char[i] = substring(y[i], nchar(y[i]) - 2, nchar(y[i]))
		last5char[i] = substring(y[i], nchar(y[i]) - 4, nchar(y[i]))
	}
	if (sum(last3char == "End") != sum(last5char == "Start")) {
		stop("Error in format of  time course data: a Start or End tag is missing")
	}
	y.act = rep(NA, n)
	timey = rep(NA, n)
	person.id = rep(NA, n)
	k = 1
	end.flag = FALSE
	person.id[1] = 1
	if (substring(y[1], nchar(y[1]) - 4, nchar(y[1])) != "Start") {
		stop("Error in format of  time course data: first cell should have a Start tag")
	}
	for (i in 1:n) {
		cat(i)
		j = 1
		while (substring(y[i], j, j) != "T") {
			j = j + 1
		}
		end.of.y = j - 1
		y.act[i] = as.numeric(substring(y[i], 1, end.of.y))
		timey[i] = substring(y[i], end.of.y + 5, nchar(y[i]))
		if (nchar(timey[i]) > 3 & substring(timey[i], nchar(timey[i]) - 
			2, nchar(timey[i])) == "End") {
			end.flag = TRUE
			timey[i] = substring(timey[i], 1, nchar(timey[i]) - 
				3)
		}
		if (nchar(timey[i]) > 3 & substring(timey[i], nchar(timey[i]) - 
			4, nchar(timey[i])) == "Start") {
			timey[i] = substring(timey[i], 1, nchar(timey[i]) - 
				5)
		}
		if (i < n & !end.flag) {
			person.id[i + 1] = k
		}
		if (i < n & end.flag) {
			k = k + 1
			person.id[i + 1] = k
		}
		end.flag = FALSE
	}
	timey = as.numeric(timey)
	# do a check that the format was correct
	tt = table(person.id, y.act)
	junk = function(x) {
		sum(x != 0)
	}
	if (sum(apply(tt, 1, junk) != 1) > 0) {
		num = (1:nrow(tt))[apply(tt, 1, junk) > 1]
		stop(paste("Error in format of  time course data, timecourse #", 
			as.character(num)))
	}
	npeople = length(unique(person.id))
	newx = matrix(NA, nrow = nrow(x), ncol = npeople)
	sd = matrix(NA, nrow = nrow(x), ncol = npeople)
	for (j in 1:npeople) {
		jj = person.id == j
		tim = timey[jj]
		xc = t(scale(t(x[, jj, drop = F]), center = TRUE, scale = FALSE))
		if (time.summary.type == "slope") {
			junk = quantitative.func(xc, tim - mean(tim))
			newx[, j] = junk$numer
			sd[, j] = junk$sd
		}
		if (time.summary.type == "signed.area") {
			junk = timearea.func(x[, jj, drop = F], tim)
			newx[, j] = junk$numer
			sd[, j] = junk$sd
		}
	}
	y.unique = y.act[!duplicated(person.id)]
	return(list(y = y.unique, x = newx, sd = sd))
}

check.format = function(y, resp.type, censoring.status = NULL) {
	# here i do some format checks for the input data$y
	# note that checks for time course data are done in the
	#   parse function for time course;
	#  we then check the output from the parser in this function
	if (resp.type == samr.const.twoclass.unpaired.response | 
		resp.type == samr.const.twoclass.unpaired.timecourse.response) {
		if (sum(y == 1) + sum(y == 2) != length(y)) {
			stop(paste("Error in input response data: response type ", 
				resp.type, " specified; values must be 1 or 2"))
		}
	}
	if (resp.type == samr.const.twoclass.paired.response | resp.type == 
		samr.const.twoclass.paired.timecourse.response) {
		if (sum(y) != 0) {
			stop(paste("Error in input response data: response type ", 
				resp.type, " specified; values must be -1, 1, -2, 2, etc"))
		}
		if (sum(table(y[y > 0]) != abs(table(y[y < 0])))) {
			stop(paste("Error in input response data:  response type ", 
				resp.type, " specified; values must be -1, 1, -2, 2, etc"))
		}
	}
	if (resp.type == samr.const.oneclass.response | resp.type == 
		samr.const.oneclass.timecourse.response) {
		if (sum(y == 1) != length(y)) {
			stop(paste("Error in input response data: response type ", 
				resp.type, " specified;  values must all be 1"))
		}
	}
	if (resp.type == samr.const.multiclass.response) {
		tt = table(y)
		nc = length(tt)
		if (sum(y <= nc & y > 0) < length(y)) {
			stop(paste("Error in input response data: response type ", 
				resp.type, " specified; values must be 1,2, ... number of classes"))
		}
		for (k in 1:nc) {
			if (sum(y == k) < 2) {
				stop(paste("Error in input response data: response type ", 
				  resp.type, " specified; there must be >1 sample per class"))
			}
		}
	}
	if (resp.type == samr.const.quantitative.response) {
		if (!is.numeric(y)) {
			stop(paste("Error in input response data: response type", 
				resp.type, " specified; values must be numeric"))
		}
	}
	if (resp.type == samr.const.survival.response) {
		if (is.null(censoring.status)) {
			stop(paste("Error in input response data: response type ", 
				resp.type, " specified; error in censoring indicator"))
		}
		if (!is.numeric(y) | sum(y < 0) > 0) {
			stop(paste("Error in input response data:  response type ", 
				resp.type, " specified; survival times  must be numeric and nonnegative"))
			if (sum(censoring.status == 0) + sum(censoring.status == 
				1) != length(censoring.status)) {
				stop(paste("Error in input response data: response type ", 
				  resp.type, " specified; censoring indicators must be 0 (censored) or 1 (failed)"))
			}
		}
		if (sum(censoring.status == 1) < 1) {
			stop(paste("Error in input response data:   response type ", 
				resp.type, " specified; there are no uncensored observations"))
		}
	}
	return()
}

mysvd <- function(x, n.components = NULL) {
	# finds PCs of matrix x
	p <- nrow(x)
	n <- ncol(x)
	# center the observations (rows)
	feature.means <- rowMeans(x)
	x <- t(scale(t(x), center = feature.means, scale = F))
	if (is.null(n.components)) {
		n.components = min(n, p)
	}
	if (p > n) {
		a <- eigen(t(x) %*% x)
		v <- a$vec[, 1:n.components, drop = FALSE]
		d <- sqrt(a$val[1:n.components, drop = FALSE])
		u <- scale(x %*% v, center = FALSE, scale = d)
		return(list(u = u, d = d, v = v))
	}
	else {
		junk <- svd(x, LINPACK = TRUE)
		nc = min(ncol(junk$u), n.components)
		return(list(u = junk$u[, 1:nc], d = junk$d[1:nc], v = junk$v[, 
			1:nc]))
	}
}

permute.rows <- function(x) {
	dd <- dim(x)
	n <- dd[1]
	p <- dd[2]
	mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
	matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

samr.assess.samplesize = function(samr.obj, data, 
	dif, samplesize.factors = c(1, 2, 3, 5), min.genes = 10, 
	max.genes = nrow(data$x)/2) {
	if (length(samplesize.factors) > 4) {
		stop("Length of samplesize.factors must be less than or equal to 4")
	}
	n.ssf = length(samplesize.factors)
	if (samr.obj$resp.type != samr.const.twoclass.unpaired.response & 
		samr.obj$resp.type != samr.const.twoclass.paired.response & 
		samr.obj$resp.type != samr.const.oneclass.response & 
		samr.obj$resp.type != samr.const.survival.response) {
		stop("Function only implemented for  twoclass.unpaired, twoclass.paired,\noneclass and survival data types")
	}
	m = nrow(data$x)
	n = ncol(data$x)
	if (samr.obj$resp.type == samr.const.twoclass.unpaired.response) {
		n1 = sum(data$y == 1)
		n2 = sum(data$y == 2)
	}
	if (samr.obj$resp.type == samr.const.twoclass.paired.response) {
		n1 = n/2
		n2 = n/2
	}
	nreps = 3
	klist = round(exp(seq(log(min.genes), log(max.genes), length = 10)))
	#power=rep(NA,length(klist))
	#type1=rep(NA,length(klist))
	fdr = matrix(NA, nrow = length(klist), ncol = n.ssf)
	fdr90 = matrix(NA, nrow = length(klist), ncol = n.ssf)
	fdr10 = matrix(NA, nrow = length(klist), ncol = n.ssf)
	fnr = matrix(NA, nrow = length(klist), ncol = n.ssf)
	fnr90 = matrix(NA, nrow = length(klist), ncol = n.ssf)
	fnr10 = matrix(NA, nrow = length(klist), ncol = n.ssf)
	cutp = matrix(NA, nrow = length(klist), ncol = n.ssf)
	#sd=samr.obj$sd-samr.obj$s0
	sd = samr.obj$sd
	if (samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
		samr.obj$resp.type == samr.const.twoclass.paired.response) {
		sigma = sd/sqrt(1/n1 + 1/n2)
		difm = dif/(sigma * sqrt(1/n1 + 1/n2))
	}
	if (samr.obj$resp.type == samr.const.oneclass.response | 
		samr.obj$resp.type == samr.const.survival.response) {
		sigma = (sqrt(n)) * sd
		difm = sqrt(n) * dif/sigma
	}
	# we only use the  20 of the perms, for speed
	nperms = min(20, samr.obj$nperms.act)
	perms.to.use = sample(1:samr.obj$nperms.act, size = nperms)
	#note: here I permute within each col of zstar, so that the
	# genes that are modified are different for each
	#   permutation
	zstar0 = t(permute.rows(t(samr.obj$ttstar0[, perms.to.use])))
	ii = 0
	for (k in klist) {
		ii = ii + 1
		oo = sample(1:m, size = k)
		temp = matrix(F, nrow = nrow(zstar0), ncol = ncol(zstar0))
		temp[oo, ] = T
		for (kk in 1:n.ssf) {
			zstar = zstar0
			zstar[oo, ] = zstar[oo, ] + difm[oo] * sqrt(samplesize.factors[kk])
			cutp[ii, kk] = quantile(abs(zstar), 1 - (k/m))
			temp[oo, ] = T
			#
			#   power[ii]=(sum(abs(zstar[oo,])>cutp[ii,kk])/samr.obj$nperms)/k
			#
			#   type1[ii]=(sum(abs(zstar[-oo,])>cutp[ii,kk])/samr.obj$nperms)/(m-k)
			u0 = colSums(abs(zstar) > cutp[ii, kk] & !temp)
			r0 = colSums(abs(zstar) > cutp[ii, kk])
			oo2 = !is.na(u0/r0)
			fdr[ii, kk] = median((u0/r0)[oo2])
			fdr90[ii, kk] = quantile((u0/r0)[oo2], 0.9)
			fdr10[ii, kk] = quantile((u0/r0)[oo2], 0.1)
			v0 = colSums(abs(zstar) < cutp[ii, kk] & temp)
			w0 = colSums(abs(zstar) < cutp[ii, kk])
			oo3 = !is.na(v0/w0)
			fnr[ii, kk] = median((v0/w0)[oo3])
			fnr90[ii, kk] = quantile((v0/w0)[oo3], 0.9)
			fnr10[ii, kk] = quantile((v0/w0)[oo3], 0.1)
		}
	}
	ng = round(klist)
	oo = !duplicated(ng)
	results = array(NA, c(sum(oo), 8, n.ssf))
	for (kk in 1:n.ssf) {
		results[, , kk] = cbind(ng, cutp[, kk], fdr[, kk], fdr90[, 
			kk], fdr10[, kk], fnr[, kk], fnr90[, kk], fnr10[, 
			kk])[oo, ]
	}
	dimnames(results) = list(NULL, c("number of genes", "cutpoint", 
		"FDR,1-power", "FDR90", "FDR10", "FNR,type1 error", "FNR90", 
		"FNR10"), as.character(samplesize.factors))
	return(list(results = results, dif.call = dif, difm = mean(difm), 
		samplesize.factors = samplesize.factors, n = ncol(data$x)))
}

samr.assess.samplesize.plot <- function(samr.assess.samplesize.obj, 
	logx = TRUE, call.win.metafile = FALSE) {
	n.ssf = length(samr.assess.samplesize.obj$samplesize.factors)
	if (call.win.metafile) {
		win.metafile()
	}
	if (n.ssf == 1) {
		par(mfrow = c(1, 1))
	}
	if (n.ssf == 2) {
		par(mfrow = c(1, 2))
	}
	if (n.ssf > 2) {
		par(mfrow = c(2, 2))
	}
	par(oma = c(0, 0, 2, 0))
	na.min = function(x) {
		min(x[!is.na(x)])
	}
	na.max = function(x) {
		max(x[!is.na(x)])
	}
	temp = samr.assess.samplesize.obj$results
	ymax = max(c(temp[, "FDR,1-power", ], temp[, "FDR90", ], 
		temp[, "FDR10", ], temp[, "FNR,type1 error", ], temp[, 
			"FDR90", ], temp[, "FDR10", ]))
	for (kk in 1:n.ssf) {
		results = samr.assess.samplesize.obj$results[, , kk]
		if (logx) {
			plot(results[, "number of genes"], results[, "FDR,1-power"], 
				log = "x", xlab = "Number of genes", ylab = "", 
				type = "n", ylim = c(0, ymax))
		}
		if (!logx) {
			plot(results[, "number of genes"], results[, "FDR,1-power"], 
				xlab = "Number of genes", ylab = "", type = "n", 
				ylim = c(0, ymax))
		}
		lines(results[, "number of genes"], results[, "FDR,1-power"], 
			col = 2, type = "b", pch = 19)
		lines(results[, "number of genes"], results[, "FDR90"], 
			col = 2, lty = 2, pch = 19)
		lines(results[, "number of genes"], results[, "FDR10"], 
			col = 2, lty = 2, pch = 19)
		lines(results[, "number of genes"], results[, "FNR,type1 error"], 
			col = 3, type = "b", pch = 19)
		lines(results[, "number of genes"], results[, "FNR90"], 
			col = 3, lty = 2, pch = 19)
		lines(results[, "number of genes"], results[, "FNR10"], 
			col = 3, lty = 2, pch = 19)
		mtext("FDR, 1-Power", side = 2, col = 2, cex = 0.8)
		mtext("FNR, Type 1 error", side = 4, col = 3, cex = 0.8)
		abline(h = 0.05, lty = 3)
		fac = samr.assess.samplesize.obj$samplesize.factors[kk]
		n = samr.assess.samplesize.obj$n
		title(paste("Sample size=", round(n * fac, 0)), cex = 0.7)
	}
	title(paste("Results for mean difference=", round(samr.assess.samplesize.obj$dif.call, 
		2)), outer = T)
	if (call.win.metafile) {
		dev.off()
	}
	return()
}

samr.pvalues.from.perms = function(tt, ttstar) {
	r = rank(c(abs(tt), abs(as.vector(ttstar))))[1:length(tt)]
	r2 = rank(c(abs(tt)))
	r3 = r - r2
	pv = (length(tt) - r3/ncol(ttstar) + 1)/length(tt)
	return(pv)
}

samr.tail.strength = function(samr.obj) {
	tt = samr.obj$tt
	ttstar = samr.obj$ttstar0
	pv = samr.pvalues.from.perms(tt, ttstar)
	m = length(pv)
	pvs = sort(pv)
	ts = (1/m) * sum((1 - pvs * (m + 1)/(1:m)))
	res = NULL
	nperms = min(ncol(ttstar), 20)
	ttstar.temp = ttstar[, 1:nperms]
	for (i in 1:nperms) {
		cat(i, fill = T)
		pvstar = samr.pvalues.from.perms(ttstar.temp[, i], ttstar.temp[, 
			-i])
		pvstar = sort(pvstar)
		tsstar = (1/m) * sum((1 - pvstar * (m + 1)/(1:m)))
		res = c(res, tsstar)
	}
	se.ts.perm = sqrt(var(res))
	return(list(ts = ts, se.ts = se.ts.perm))
}

quant <- function(x) {
	#dyn.load('/home/tibs/PAPERS/copa/quant.so')
	p = as.integer(nrow(x))
	n = as.integer(ncol(x))
	xx = t(x)
	storage.mode(xx) = "single"
	junk = .Fortran("quant", as.matrix(xx), n, p, integer(n), 
		out = single(2 * p), PACKAGE = "samr")
	return(matrix(junk$out, ncol = 2))
}

################
#new functions for SAMseq
######################################################################
#\t\tEstimate sequencing depths
#\tArguments:
#\t\tx: data matrix. nrow=#gene, ncol=#sample
#\tValue:
# depth: estimated sequencing depth. a vector with len
#   #sample.
######################################################################
samr.estimate.depth <- function(x) {
	iter <- 5
	cmeans <- colSums(x)/sum(x)
	for (i in 1:iter) {
		n0 <- rowSums(x) %*% t(cmeans)
		prop <- rowSums((x - n0)^2/(n0 + 1e-08))
		qs <- quantile(prop, c(0.25, 0.75))
		keep <- (prop >= qs[1]) & (prop <= qs[2])
		cmeans <- colMeans(x[keep, ])
		cmeans <- cmeans/sum(cmeans)
	}
	depth <- cmeans/mean(cmeans)
	return(depth)
}

######################################################################
#\t\tResampling
#\tArguments:
#\t\tx: data matrix. nrow=#gene, ncol=#sample
#\t\td: estimated sequencing depth
#\t\tnresamp: number of resamplings
#\tValue:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
######################################################################
resample <- function(x, d, nresamp = 20) {
	ng <- nrow(x)
	ns <- ncol(x)
	dbar <- exp(mean(log(d)))
	xresamp <- array(0, dim = c(ng, ns, nresamp))
	for (k in 1:nresamp) {
		for (j in 1:ns) {
			xresamp[, j, k] <- rpois(n = ng, lambda = (dbar/d[j]) * 
				x[, j]) + runif(ng) * 0.1
		}
	}
	for (k in 1:nresamp) {
		xresamp[, , k] <- t(rankcol(t(xresamp[, , k])))
	}
	return(xresamp)
}

######################################################################
#\t\tTwoclass Wilcoxon statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of values 1 and 2
#\tValue:
#\t\ttt: the statistic.
######################################################################
#
wilcoxon.unpaired.seq.func <- function(xresamp, y) {
	tt <- rep(0, dim(xresamp)[1])
	for (i in 1:dim(xresamp)[3]) {
		tt <- tt + rowSums(xresamp[, y == 2, i]) - sum(y == 2) * 
			(length(y) + 1)/2
	}
	tt <- tt/dim(xresamp)[3]
	return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

######################################################################
#\t\tTwoclass paired Wilcoxon statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of values +1, -1, +2, -2, ...
#\tValue:
#\t\ttt: the statistic.
######################################################################
wilcoxon.paired.seq.func <- function(xresamp, y) {
	tt <- rep(0, dim(xresamp)[1])
	for (i in 1:dim(xresamp)[3]) {
		tt <- tt + rowSums(xresamp[, y > 0, i]) - sum(y > 0) * 
			(length(y) + 1)/2
	}
	tt <- tt/dim(xresamp)[3]
	return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

######################################################################
#\t\tMulticlass statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of values 1, 2, ..., K
#\tValue:
#\t\ttt: the statistic.
######################################################################
multiclass.seq.func <- function(xresamp, y)
{
	# number of classes and number of samples in each class
	K <- max(y)
	n.each <- rep(0, K)
	for (k in 1 : K)
	{
		n.each[k] <- sum(y == k)
	}
	# the statistic
	tt <- temp <- rep(0, dim(xresamp)[1])
	stand.contrasts <- matrix(0, dim(xresamp)[1], K)
	
	for (i in 1 : dim(xresamp)[3])
	{
		for (k in 1 : K)
		{
			temp <- rowSums(xresamp[, y == k, i])
			tt <- tt + temp ^2 / n.each[k]
			stand.contrasts[, k] <- stand.contrasts[, k] + temp
		}
	}
	# finalize
	nresamp <- dim(xresamp)[3]
	ns <- dim(xresamp)[2]
	tt <- tt / nresamp * 12 / ns / (ns + 1) - 3 * (ns + 1)
	stand.contrasts <- stand.contrasts / nresamp
	stand.contrasts <- scale(stand.contrasts, center=n.each * (ns + 1) / 2, 
		scale=sqrt(n.each * (ns - n.each) * (ns + 1) / 12))
	return(list(tt = tt, numer = tt, sd = rep(1, length(tt)), 
		stand.contrasts = stand.contrasts))
}

## Jun commented this function
#######################################################################
##\t\tQuantitative statistics
##\tArguments:
##\t\txresamp: an rank array with dim #gene*#sample*nresamp
##\t\ty: outcome vector of real values
##\tValue:
##\t\ttt: the statistic.
#######################################################################
#quantitative.seq.func <- function(xresamp, y)
#{
#\ty.ranked <- rank(y) - (dim(xresamp)[2] + 1) / 2
#
#\ttt <- rep(0, dim(xresamp)[1])
#
#\tfor (i in 1 : dim(xresamp)[3])
#\t{
# tt <- tt + (xresamp[, , i] - (dim(xresamp)[2] + 1) / 2)
#   %*% y.ranked
#\t}
#
#\tns <- dim(xresamp)[2]
#\ttt <- tt / (dim(xresamp)[3] * (ns ^ 3 - ns) / 12)
#
#
#   return(list(tt=as.vector(tt),numer=as.vector(tt),sd=rep(1,length(tt))))
#}

## Jun added starts
######################################################################
#\t\tQuantitative statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of real values
#\tValue:
#\t\ttt: the statistic.
######################################################################
quantitative.seq.func <- function(xresamp, y) {
	tt <- rep(0, dim(xresamp)[1])
	for (i in 1:dim(xresamp)[3]) {
		y.ranked <- rank(y, ties.method = "random") - (dim(xresamp)[2] + 
			1)/2
		tt <- tt + (xresamp[, , i] - (dim(xresamp)[2] + 1)/2) %*% 
			y.ranked
	}
	ns <- dim(xresamp)[2]
	tt <- tt/(dim(xresamp)[3] * (ns^3 - ns)/12)
	return(list(tt = as.vector(tt), numer = as.vector(tt), sd = rep(1, 
		length(tt))))
}
## Jun added ends

######################################################################
#\t\tSurvival statistics
#\tArguments:
#\t\txresamp: an rank array with dim #gene*#sample*nresamp
#\t\ty: outcome vector of real values
#\t\tcensoring.status: 1=died, 0=censored
#\tValue:
#\t\ttt: the statistic.
######################################################################
cox.seq.func <- function(xresamp, y, censoring.status) {
	# get the dimensions
	ng <- dim(xresamp)[1]
	ns <- dim(xresamp)[2]
	# prepare for the calculation
	# find the index matrix
	Dn <- sum(censoring.status == 1)
	Dset <- c(1:ns)[censoring.status == 1]  # the set of died
	ind <- matrix(0, ns, Dn)
	# get the matrix
	for (i in 1:Dn) {
		ind[y >= y[Dset[i]] - 1e-08, i] <- 1/sum(y >= y[Dset[i]] - 
			1e-08)
	}
	ind.sums <- rowSums(ind)
	# calculate the score statistic
	tt <- apply(xresamp, 3, function(x, cen.ind, ind.para, ind.sums.para) {
		dev1 <- x %*% cen.ind
		x.ind <- x %*% ind.para
		dev2 <- (x * x) %*% ind.sums.para - rowSums(x.ind * x.ind)
		dev1/(sqrt(dev2) + 1e-08)
	}, (censoring.status - ind.sums), ind, ind.sums)
	tt <- rowMeans(tt)
	return(list(tt = tt, numer = tt, sd = rep(1, length(tt))))
}

rankcol = function(x) {
	# ranks the elements within each col of the matrix x
	# and returns these ranks in a matrix
	n = nrow(x)
	p = ncol(x)
	mode(n) = "integer"
	mode(p) = "integer"
	mode(x) = "single"
	if (!is.loaded("rankcol")) {
		#dyn.load('/home/tibs/PAPERS/jun2/test/rankcol.so')
	}
	junk = .Fortran("rankcol", x, n, p, xr = integer(n * p), 
		integer(n), PACKAGE = "samr")
	xr = matrix(junk$xr, nrow = n, ncol = p)
	return(xr)
}

## Jun added starts
######################################################################
#\t\tfoldchange of twoclass unpaired sequencing data
######################################################################
foldchange.seq.twoclass.unpaired <- function(x, y, depth)
{
	require("matrixStats")
	x.norm <- scale(x, center = F, scale = depth) + 1e-08
	fc <- rowMedians(x.norm[, y == 2])/rowMedians(x.norm[, y == 
		1])
	return(fc)
}

######################################################################
#\t\tfoldchange of twoclass paired sequencing data
######################################################################
foldchange.seq.twoclass.paired <- function(x, y, depth) {
	require("matrixStats")
	nc <- ncol(x)/2
	o1 <- o2 <- rep(0, nc)
	for (j in 1:nc) {
		o1[j] <- which(y == -j)
		o2[j] <- which(y == j)
	}
	x.norm <- scale(x, center = F, scale = depth) + 1e-08
	d <- x.norm[, o2, drop = F]/x.norm[, o1, drop = F]
	fc <- rowMedians(d, na.rm = T)
	return(fc)
}
## Jun added ends

## Jun added starts
#######################################################################
#\tcompute the delta table for sequencing data
#######################################################################
samr.compute.delta.table.seq <- function(samr.obj, 
	min.foldchange = 0, dels = NULL) {
	res1 <- NULL
	flag <- T
	## check whether any gene satisfies the foldchange
	#   restrictions
	if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
		samr.obj$resp.type == samr.const.twoclass.paired.response) & 
		(min.foldchange > 0)) {
		sat.up <- (samr.obj$foldchange >= min.foldchange) & (samr.obj$evo > 
			0)
		sat.dn <- (samr.obj$foldchange <= 1/min.foldchange) & 
			(samr.obj$evo < 0)
		if (sum(sat.up) + sum(sat.dn) == 0) {
			flag <- F
		}
	}
	if (flag) {
		if (is.null(dels)) {
			dels <- generate.dels(samr.obj, min.foldchange = min.foldchange)
		}
		cat("Number of thresholds chosen (all possible thresholds) =", 
			length(dels), fill = T)
		if (length(dels) > 0) {
			## sort delta to make the fast calculation right
			dels <- sort(dels)
			## get the upper and lower cutoffs
			cat("Getting all the cutoffs for the thresholds...\n")
			slabs <- samr.seq.detec.slabs(samr.obj, dels, min.foldchange)
			cutup <- slabs$cutup
			cutlow <- slabs$cutlow
			g2 <- slabs$g2
			## get the number of errors under the null hypothesis
			cat("Getting number of false positives in the permutation...\n")
			errnum <- samr.seq.null.err(samr.obj, min.foldchange, 
				cutup, cutlow)
			res1 <- NULL
			gmed <- apply(errnum, 2, median)
			g90 = apply(errnum, 2, quantile, 0.9)
			res1 <- cbind(samr.obj$pi0 * gmed, samr.obj$pi0 * 
				g90, g2, samr.obj$pi0 * gmed/g2, samr.obj$pi0 * 
				g90/g2, cutlow, cutup)
			res1 <- cbind(dels, res1)
			dimnames(res1) <- list(NULL, c("delta", "# med false pos", 
				"90th perc false pos", "# called", "median FDR", 
				"90th perc FDR", "cutlo", "cuthi"))
		}
	}
	return(res1)
}

######################################################################
#\tget the number of significance in the null distribution
######################################################################
samr.seq.null.err <- function(samr.obj, min.foldchange, 
	cutup, cutlow) {
	errup = matrix(NA, ncol = length(cutup), nrow = ncol(samr.obj$ttstar0))
	errlow = matrix(NA, ncol = length(cutlow), nrow = ncol(samr.obj$ttstar0))
	cutup.rank <- rank(cutup, ties.method = "min")
	cutlow.rank <- rank(-cutlow, ties.method = "min")
	for (jj in 1:ncol(samr.obj$ttstar0)) {
		#cat(jj, fill=TRUE)
		keep.up <- keep.dn <- samr.obj$ttstar0[, jj]
		if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
			samr.obj$resp.type == samr.const.twoclass.paired.response) & 
			(min.foldchange > 0)) {
			keep.up <- keep.up[samr.obj$foldchange.star[, jj] >= 
				min.foldchange]
			keep.dn <- keep.dn[samr.obj$foldchange.star[, jj] <= 
				1/min.foldchange]
		}
		errup[jj, ] <- length(keep.up) - (rank(c(cutup, keep.up), 
			ties.method = "min")[1:length(cutup)] - cutup.rank)
		errlow[jj, ] <- length(keep.dn) - (rank(c(-cutlow, -keep.dn), 
			ties.method = "min")[1:length(cutlow)] - cutlow.rank)
	}
	errnum <- errup + errlow
	return(errnum)
}

######################################################################
#\tdetect multiple slabs
######################################################################
samr.seq.detec.slabs <- function(samr.obj, dels, min.foldchange) {
	## initialize calculation
	tag <- order(samr.obj$tt)
	if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
		samr.obj$resp.type == samr.const.twoclass.paired.response) & 
		(min.foldchange > 0)) {
		res.mat <- data.frame(tt = samr.obj$tt[tag], fc = samr.obj$foldchange[tag], 
			evo = samr.obj$evo, dif = samr.obj$tt[tag] - samr.obj$evo)
		res.up <- res.mat[res.mat$evo > 0, ]
		res.lo <- res.mat[res.mat$evo < 0, ]
		res.up <- res.up[res.up$fc >= min.foldchange, ]
		res.lo <- res.lo[res.lo$fc <= 1/min.foldchange, ]
	}
	else {
		res.mat <- data.frame(tt = samr.obj$tt[tag], evo = samr.obj$evo, 
			dif = samr.obj$tt[tag] - samr.obj$evo)
		res.up <- res.mat[res.mat$evo > 0, ]
		res.lo <- res.mat[res.mat$evo < 0, ]
	}
	## begin calculating
	cutup <- rep(1e+10, length(dels))
	cutlow <- rep(-1e+10, length(dels))
	g2.up <- g2.lo <- rep(0, length(dels))
	if (nrow(res.up) > 0) {
		res.up <- data.frame(dif = res.up$dif, tt = res.up$tt, 
			num = nrow(res.up):1)
		## get the upper part
		j <- 1
		ii <- 1
		while (j <= nrow(res.up) & ii <= length(dels)) {
			if (res.up$dif[j] > dels[ii]) {
				cutup[ii] <- res.up$tt[j]
				g2.up[ii] <- res.up$num[j]
				ii <- ii + 1
			}
			else {
				j <- j + 1
			}
		}
	}
	if (nrow(res.lo) > 0) {
		res.lo <- data.frame(dif = res.lo$dif, tt = res.lo$tt, 
			num = 1:nrow(res.lo))
		## get the lower part
		j <- nrow(res.lo)
		ii <- 1
		while (j >= 1 & ii <= length(dels)) {
			if (res.lo$dif[j] < -dels[ii]) {
				cutlow[ii] <- res.lo$tt[j]
				g2.lo[ii] <- res.lo$num[j]
				ii <- ii + 1
			}
			else {
				j <- j - 1
			}
		}
	}
	g2 <- g2.up + g2.lo
	return(list(cutup = cutup, cutlow = cutlow, g2 = g2))
}

#######################################################################
#\tcompute the delta table for sequencing data
#######################################################################
generate.dels <- function(samr.obj, min.foldchange = 0) {
	dels <- NULL
	## initialize calculation
	tag <- order(samr.obj$tt)
	if ((samr.obj$resp.type == samr.const.twoclass.unpaired.response | 
		samr.obj$resp.type == samr.const.twoclass.paired.response) & 
		(min.foldchange > 0)) {
		res.mat <- data.frame(tt = samr.obj$tt[tag], fc = samr.obj$foldchange[tag], 
			evo = samr.obj$evo, dif = samr.obj$tt[tag] - samr.obj$evo)
		res.up <- res.mat[res.mat$evo > 0, ]
		res.lo <- res.mat[res.mat$evo < 0, ]
		res.up <- res.up[res.up$fc >= min.foldchange, ]
		res.lo <- res.lo[res.lo$fc <= 1/min.foldchange, ]
	}
	else {
		res.mat <- data.frame(tt = samr.obj$tt[tag], evo = samr.obj$evo, 
			dif = samr.obj$tt[tag] - samr.obj$evo)
		res.up <- res.mat[res.mat$evo > 0, ]
		res.lo <- res.mat[res.mat$evo < 0, ]
	}
	## for the upper part
	up.vec <- rep(NA, nrow(res.up))
	if (nrow(res.up) > 0) {
		st <- 1e-08
		i.cur <- 1
		for (i in 1:nrow(res.up)) {
			if (res.up$dif[i] > st) {
				st <- res.up$dif[i]
				up.vec[i.cur] <- st
				i.cur <- i.cur + 1
			}
		}
	}
	## for the lower part
	lo.vec <- rep(NA, nrow(res.lo))
	if (nrow(res.lo) > 0) {
		st <- -1e-08
		i.cur <- 1
		for (i in nrow(res.lo):1) {
			if (res.lo$dif[i] < st) {
				st <- res.lo$dif[i]
				lo.vec[i.cur] <- st
				i.cur <- i.cur + 1
			}
		}
	}
	## combine them
	vec <- c(up.vec, -lo.vec)
	vec <- vec[!is.na(vec)]
	vec <- vec - 1e-08
	dels <- sort(unique(vec))
	return(dels)
}

#######################################################################
#\tgenerate normalized data
#######################################################################
samr.norm.data <- function(x, depth = NULL) {
	if (is.null(depth)) {
		depth <- samr.estimate.depth(x)
	}
	scaling.factors <- (prod(depth)^(1/length(depth)))/depth
	x <- 2 * sqrt(x + 3/8)
	x <- scale(x, center = F, scale = sqrt(1/scaling.factors))
	return(x)
}
## Jun added ends
