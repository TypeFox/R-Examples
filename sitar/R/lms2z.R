	LMS2z <- function(x, y, sex, measure, ref, toz=TRUE) {
#	converts measurement y to/from z-score adjusted for x & sex
#		using LMS reference 'ref' for 'measure'
#	x		age
#	y		measurement (or z-score if toz FALSE)
#	sex		sex variable (male=1, female=2)
#	measure label for measurement, one of:
#		'ht' 'wt' 'bmi' 'head' 'sitht' 'leglen' 'waist' 'bfat'
#	ref		name of reference, one of: 'uk90' 'who06'
#	toz		if TRUE returns measurement converted to z-score using ref
#		  	if FALSE returns z-score converted to measurement using ref
	if (!length(sex) %in% c(1, length(x))) stop('sex wrong length for x')
	lms <- paste(c('L', 'M', 'S'), measure, sep='.')
	v <- matrix(nrow=length(x), ncol=4)
	v[, 4] <- sex
	ref <- get(ref)
	for (i in 1:3) {
		for (ix in 1:2) {
			sexvar <- as.numeric(v[, 4]) == ix
			sexref <- as.numeric(ref$sex) == ix
			if (any(sexvar)) v[sexvar, i] <- spline(ref$years[sexref],
				ref[sexref, lms[i]], method='natural', xout=x[sexvar])$y
		}
	}
	cz <- if (toz) zLMS(y, v[, 1], v[, 2], v[, 3])
		else cLMS(y, v[, 1], v[, 2], v[, 3])
	if (!is.null(dim(cz))) {
	  cz <- as.data.frame(cz)
	  if (!toz) names(cz) <- z2cent(y)
	  rownames(cz) <- x
	}
	cz
}

	zLMS <- function(x, L = 1, M, S) {
  	L0 <- L + 1e-7 * (L == 0)
  	LMS <- data.frame(cbind(L0, M, S))
  	if (length(x) == nrow(LMS) || min(length(x), nrow(LMS)) == 1)
  	  drop(with(LMS, ((x / M) ^ L0 - 1) / L0 / S))
  	else
  	  drop(with(LMS, (t(x %*% t(1 / M)) ^ L0 - 1) / L0 / S))
	}

	cLMS <- function(z, L = 1, M, S) {
	  L0 <- L + 1e-7 * (L == 0)
	  LMS <- data.frame(cbind(L0, M, S))
	  if (length(z) == nrow(LMS) || min(length(z), nrow(LMS)) == 1)
	    drop(with(LMS, M * (1 + L0 * S * z) ^ (1/L0)))
	  else
	    drop(with(LMS, M * (1 + L0 * S %*% t(z)) ^ (1/L0)))
	}

	z2cent <- function(z) {
#	z is z-score
#	returns corresponding centile as label
	np <- ifelse(abs(z) < 2.33, 0, 1)
	ct <- round(pnorm(z) * 100, np)
	mod10 <- ifelse(np == 1, 0, floor(ct %% 10))
	th <- ifelse(mod10 == 0 | mod10 > 4, 4, mod10)
	th <- paste(ct, c('st','nd','rd','th')[th], sep='')
	th[th == '0th'] <- paste('SDS', round(z[th == '0th'], 1), sep='')
	th[th == '100th'] <- paste('SDS', round(z[th == '100th'], 1), sep='+')
	th
}
