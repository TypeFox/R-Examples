###################################################
#		Permutation and get the object
###################################################
Score.Stat <- function(dat, para)
{
	cat("Calculating the score statistics...\n")
	
	type <- match(dat$type, c('twoclass', 'multiclass', 'quant'), nomatch=0)
	
	if (type == 0)
	{
		stop('data type is not \'twoclass\', \'multiclass\' or \'quant\'!')
	}
	
	ps.obj <- switch(type,
		Class.Score(n=dat$n, y=dat$y, npermu=para$npermu, pair=dat$pair),
		Class.Score(n=dat$n, y=dat$y, npermu=para$npermu, pair=F),
		Quant.Score(n=dat$n, y=dat$y, npermu=para$npermu))
	
	return(ps.obj)
}

###################################################
#		Score statistics for twoclass or multiclass response
###################################################
Class.Score <- function(n, y, npermu=100, pair=F)
{
	SMALL.VAL <- 1e-8
	
	cat("twoclass or multiclass data...\n")
	
	###################################################
	### warm up
	n0 <- rowSums(n) %*% t(Est.Depth(n)$cmeans)
	n_n0 <- n - n0
	
	###################################################
	### calculate for the original data
	K <- max(y)
	ind.ck <- matrix(0, length(y), K)
	for (k in 1 : K)
	{
		ind.ck[y == k, k] <- 1
	}
	
	tt <- rowSums((n_n0 %*% ind.ck) ^ 2 / (n0 %*% ind.ck + SMALL.VAL))
	tt.signed <- sqrt(tt) * sign((n_n0 %*% ind.ck)[, 1])

	###################################################
	### calculate for the permutation data
	### initialize permutation
	ind <- Permu.Ind(y=y, npermu=npermu, pair=pair)
	npermu.act <- nrow(ind)
	ttstar0 <- matrix(NA, nrow=nrow(n), ncol=npermu.act)
	
	### permutation
	for (i in 1 : npermu.act)
	{
		ind.ck.star <- ind.ck[ind[i, ], ]
		ttstar0[, i] <- rowSums((n_n0 %*% ind.ck.star) ^ 2 / (n0 %*% ind.ck.star + SMALL.VAL))
	}
	
	return(list(tt=tt, ttstar0=ttstar0, tt.signed=tt.signed))
}

###################################################
#		Score statistics for Quantitative response
###################################################
Quant.Score <- function(n, y, npermu=100)
{
	SMALL.VAL <- 1e-8
	
	cat("quantitative data...\n")
	
	###################################################
	### make the mean of y to be zero
	y <- y - mean(y)
	
	###################################################
	### warm up
	n0 <- rowSums(n) %*% t(Est.Depth(n)$cmeans)
	n_n0 <- n - n0
	y2 <- y ^ 2
	
	###################################################
	### calculate for the original data
	tt.signed <- (n_n0 %*% y) / (sqrt(n0 %*% y2) + SMALL.VAL)
	
	tt <- tt.signed ^ 2
	
	###################################################
	### calculate for the permutation data
	### initialize permutation
	ind <- Permu.Ind(y=y, npermu=npermu, pair=F)
	npermu.act <- nrow(ind)
	ttstar0 <- matrix(NA, nrow=nrow(n), ncol=npermu.act)
	
	### permutation
	for (i in 1 : npermu.act)
	{
		ttstar0[, i] <- ((n_n0 %*% y[ind[i, ]]) / (sqrt(n0 %*% y2[ind[i, ]]) + SMALL.VAL)) ^ 2
	}
	
	return(list(tt=tt, ttstar0=ttstar0, tt.signed=tt.signed))
}
