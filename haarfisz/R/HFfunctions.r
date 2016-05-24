"SDF1" <- 
structure(.Data = list()
, class = "data.frame"
, row.names = character(0)
)

"connect.tree" <- 
function(weights)
{
	TT <- length(weights) + 1.
	J <- logb(TT, 2.)
	k <- TT/2.
	for(j in (J - 2.):0.) {
		for(i in (k + 1.):(k + 2.^(j)))
			weights[i] <- max(weights[i], weights[2. * i - k - 2.^(j + 1.) - 1.], weights[
				2. * i - k - 2.^(j + 1.)])
		k <- k + TT * 2.^(j - J)
	}
	return(weights)
}

"denoise.poisson" <- 
function(y, meth.1 = hf.bt, cs.1 = 50, meth.2 = hf.cv, cs.2 = 50, hybrid = TRUE)
{
	n <- length(y)
	est <- rep(0., n)
	for(i in 1.:cs.1) {
		y.temp <- shift.sequence(y, i)
		x.temp <- hft(y.temp)
		x.temp.clean <- meth.1(x.temp)
		y.temp.clean <- hft.inv(x.temp.clean)
		y.clean <- shift.sequence(y.temp.clean,  - i)
		est <- (1. - 1./i) * est + 1./i * y.clean
	}
	if(hybrid) {
		est.2 <- rep(0., n)
		for(i in 1.:cs.2) {
			y.temp <- shift.sequence(y, i)
			x.temp <- hft(y.temp)
			x.temp.clean <- meth.2(x.temp)
			y.temp.clean <- hft.inv(x.temp.clean)
			y.clean <- shift.sequence(y.temp.clean,  - i)
			est.2 <- (1. - 1./i) * est.2 + 1./i * y.clean
		}
		est <- (est + est.2)/2.
	}
	return(est)
}

"hf.bt" <- 
function(x, filter.number = 1, family = "DaubExPhase", min.level = 3, noise.level = NULL)
{
	x.w <- wd(x, filter.number, family)
	coeffs <- x.w$D
	TT <- length(coeffs) + 1.
	J <- floor(logb(TT, 2.))
	y <- coeffs[1.:(TT/2.)]
	if(is.null(noise.level))
		noise.level <- median(abs(diff(y) - median(diff(y))))/(sqrt(2.) * 0.67449999999999999
			)
	thresh <- noise.level * sqrt(2. * log(TT))
	weights <- rep(1., TT - 1.)
	weights[1.:(TT - 2.^min.level)] <- abs(coeffs[1.:(TT - 2.^min.level)]) >= thresh
	weights <- connect.tree(weights)
	x.w$D <- x.w$D * weights
	x.w.r <- wr(x.w)
	return(x.w.r)
}

"hf.cv" <- 
function(x, filter.number = 1, family = "DaubExPhase", min.level = 3, type = "hard")
{
	x.w <- wd(x, filter.number, family)
	x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1.), policy = "cv", type = type)
	return(wr(x.w.t))
}

"hf.tiu" <- 
function(x, filter.number = 1, family = "DaubExPhase", min.level = 3, noise.level = 1)
{
	TT <- length(x)
	thresh <- noise.level * sqrt(2. * log(TT))
	x.w <- wd(x, filter.number, family, type = "station")
	x.w.t <- threshold(x.w, policy = "manual", value = thresh, type = "hard")
	x.w.t.r <- AvBasis(convert(x.w.t))
	return(x.w.t.r)
}

"hf.u" <- 
function(x, filter.number = 10, family = "DaubLeAsymm", min.level = 3, type = "hard")
{
	x.w <- wd(x, filter.number, family)
	x.w.t <- threshold(x.w, levels = (min.level):(x.w$nlevels - 1.), policy = "universal", type
		 = type)
	return(wr(x.w.t))
}

"hft" <- 
function(data)
{
	a <- 2.
	n <- length(data)
	nhalf <- n/2.
	J <- logb(n, 2.)
	res <- data
	sm <- rep(0., nhalf)
	det <- sm
	for(i in 1.:J) {
		sm[1.:nhalf] <- (res[2. * (1.:nhalf) - 1.] + res[2. * (1.:nhalf)])/a
		det[1.:nhalf] <- (res[2. * (1.:nhalf) - 1.] - res[2. * (1.:nhalf)])/a
		det[sm > 0.] <- det[sm > 0.]/sqrt(sm[sm > 0.])
		res[1.:nhalf] <- sm[1.:nhalf]
		res[(nhalf + 1.):n] <- det[1.:nhalf]
		n <- n/2.
		nhalf <- nhalf/2.
		sm <- 0.
		det <- 0.
	}
	nhalf <- 1.
	n <- 2.
	for(i in 1.:J) {
		sm[1.:nhalf] <- res[1.:nhalf]
		det[1.:nhalf] <- res[(nhalf + 1.):n]
		res[2. * (1.:nhalf) - 1.] <- a/2. * (sm[1.:nhalf] + det[1.:nhalf])
		res[2. * (1.:nhalf)] <- a/2. * (sm[1.:nhalf] - det[1.:nhalf])
		n <- 2. * n
		nhalf <- 2. * nhalf
	}
	return(res)
}

"hft.inv" <- 
function(data)
{
	a <- 2.
	n <- length(data)
	nhalf <- n/2.
	J <- logb(n, 2.)
	res <- data
	sm <- rep(0., nhalf)
	det <- sm
	for(i in 1.:J) {
		sm[1.:nhalf] <- (res[2. * (1.:nhalf) - 1.] + res[2. * (1.:nhalf)])/a
		det[1.:nhalf] <- (res[2. * (1.:nhalf) - 1.] - res[2. * (1.:nhalf)])/a
		res[1.:nhalf] <- sm[1.:nhalf]
		res[(nhalf + 1.):n] <- det[1.:nhalf]
		n <- n/2.
		nhalf <- nhalf/2.
	}
	nhalf <- 1.
	n <- 2.
	for(i in 1.:J) {
		sm[1.:nhalf] <- res[1.:nhalf]
		det[1.:nhalf] <- res[(nhalf + 1.):n]
		res[2. * (1.:nhalf) - 1.] <- a/2. * (sm[1.:nhalf] + det[1.:nhalf] * sqrt(sm[1.:nhalf]
			))
		res[2. * (1.:nhalf)] <- a/2. * (sm[1.:nhalf] - det[1.:nhalf] * sqrt(sm[1.:nhalf]))
		res[1.:n][res[1.:n] < 0.] <- 0.
		n <- 2. * n
		nhalf <- 2. * nhalf
	}
	return(res)
}

