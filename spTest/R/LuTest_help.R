#' @keywords internal
make_htestIso_LZ <- function(HTout, df)
{
	mmeth <- "Test of reflection and complete symmetry from Lu and Zimmerman (2005) for sampling locations on the integer grid using the perdiodogram."
	names(mmeth) <- "method"

	mpv1 <- HTout$pvalue.refl
	names(mpv1) <- "p.value.refl"
	
	if( !is.null(HTout$pvalue.comp))
	{
		 mpv2 <- HTout$pvalue.comp
		 names(mpv2) <- "p.value.comp"
	}

	obj <- list(method = mmeth, p.value.refl = mpv1, p.value.comp = mpv2)

	class(obj) <- c("htestIso")		
	return(obj)
}

#' @keywords internal
scale_coords = function(spdata)
{
	min.x <- min(spdata[,1])
	min.y <- min(spdata[,2])
	
	coords.new <- spdata[,1:2]
	coords.new[,1] <- (coords.new[,1] - min.x)+1
	coords.new[,2] <- (coords.new[,2] - min.y)+1
	
	spdata.shift <-cbind(coords.new, spdata[,3])
	
	return(spdata.shift)
	
}
#' @keywords internal
cov_lags = function(nrows, ncols)
{
	jlags <- c(0:(ncols-1), -1*(1:(ncols-1)))
	jlags <- sort(jlags)
	klags <- 0:(nrows-1)
	lags <- expand.grid(jlags, klags)
	return(cbind(lags[,1],lags[,2]))
}
#' @keywords internal
lag_dist = function(locs1, locs2)
{
	n1 <- dim(locs1)[1]
	n2 <- dim(locs2)[1]
	index <- expand.grid(1:n2,1:n1)
	index <- cbind(index[,1], index[,2])
	index <- index[,2:1]
	splags <- c()
	for(i in 1:n1)
	{ 
		lags.d <- cbind(locs1[i,1] - locs2[,1], locs1[i,2] - locs2[,2])
		splags <- rbind(splags, lags.d)
	}
	
	return(cbind(index, splags))
}
#' @keywords internal
chat_jk = function(lagvec, spdata, nrows, ncols)
{	
	N <- dim(spdata)[1]
	cj <- lagvec[1]
	ck <- lagvec[2]
	mysum <- 0
	if(cj >= 0)
	{
		for(u in 1:(ncols-cj))
		{
			for(v in 1:(nrows-ck))
			{	
				loc1 <- which(spdata[,1] == u & spdata[,2] == v)
				z1 <- spdata[loc1,3]
				loc2 <- which(spdata[,1] == u+cj & spdata[,2] == v+ck)
				z2 <- spdata[loc2,3]	
				mysum <- mysum + (z1*z2)
			}	
		}
		chat <- mysum/N	
	}
	
	if(cj < 0)
	{
		for(u in (-cj+1):ncols)
		{
			for(v in 1:(nrows-ck))
			{	
				loc1 <- which(spdata[,1] == u & spdata[,2] == v)
				z1 <- spdata[loc1,3]
				loc2 <- which(spdata[,1] == u+cj & spdata[,2] == v+ck)
				z2 <- spdata[loc2,3]	
				mysum <- mysum + (z1*z2)
			}	
		}
		chat <- mysum/N	
	}
	
	return(chat)
}
#' @keywords internal
est_cov = function(lagsmat, spdata, nrows, ncols)
{
	chat <- apply(lagsmat, 1, FUN = chat_jk, spdata, nrows, ncols)
	return(chat)	
}
#' @keywords internal
cov_complete = function(chat, lags, ret.mat = FALSE)
{
	njlags <- length(unique(lags[,1]))
	klags.computed <- length(unique(lags[,2]))
	nklags <- 2*length(unique(lags[,2])) - 1
	Cmat <- matrix(data = NA, ncol = njlags, nrow = nklags)
	jlags <- unique(lags[,1])
	klags <- unique(lags[,2])
	zero <- which(klags == 0)
	tmp <- klags[-zero]
	klags <- sort( c(klags, -1*tmp) )
	klags <- rev(klags)
	mat.cols <- 1:njlags
	mat.rows <- 1:nklags
	
	for(i in 1:length(chat))
	{
		cj <- lags[i,1]
		ck <- lags[i,2]
		cloc <- which(jlags == cj)
		rloc <- which(klags == ck)
		row.coord <- mat.rows[rloc]
		col.coord <- mat.cols[cloc]
		Cmat[row.coord,col.coord] <-  chat[i]
	}
	
	lags2 <-  c()
	chat2 <-  c()
	for(a in 1:njlags)
	{
		for(b in (klags.computed+1):nklags) 
		{
			a.map <- which(mat.cols == a)
			j.lag <- jlags[a.map]
			b.map <- which(mat.rows == b)
			k.lag <- klags[b.map]
			d.map <- which(klags == -k.lag)
			d <- mat.rows[d.map]
			e.map <- which(jlags == -j.lag)
			e <- mat.cols[e.map]
			Cmat[b,a] <- Cmat[d,e]
			chat2 <- c(chat2, Cmat[d,e])
			lags2 <- rbind(lags2, c(j.lag, k.lag))
		}
	}
	
	chat.full <- cbind(rbind(lags, lags2), c(chat, chat2))
	
	if(ret.mat == T)
	{
		list("covmat" = Cmat, "lags" = rbind(lags, lags2))
	}
	else
	{
		return(chat.full)
	}
}
#' @keywords internal
get_Fourier_freqs = function(nrowss, ncolss)
{
	s <- 1:ncolss
	p <- c( -1*floor((s-1)/2 ), floor(s/2) )
	p <- sort( unique(p) )
	om1 <- 2*pi*p/ncolss
	
	s <- 1:nrowss
	q <- c( -1*floor((s-1)/2 ), floor(s/2) )
	q <- sort( unique(q) )
	om2 <- 2*pi*q/nrowss
	
	ffs <- expand.grid(om1, om2)
	return(cbind(ffs[,1], ffs[,2]))
}
#' @keywords internal
periodogram = function(cdata, freqs)
{
	nfreqs <- dim(freqs)[1]
	ivec <- c()
	for(i in 1:nfreqs)
	{
		mysum <- 0
		mysum <- sum( cdata[,3]*cos(freqs[i,1]*cdata[,1] + freqs[i,2]*cdata[,2]) )
		I <- mysum/((2*pi)^2)
		ivec <- c(ivec, I)
	}
	
	return(cbind(freqs, ivec))
}
#' @keywords internal
test_reflection_sym = function(pdata, nrows, ncols, nsim = 5000)
{
	V <- rbind(c(0,0), c(0, pi), c(pi,0), c(pi, pi))
	S <- pdata[which(pdata[,2] >= 0),c(1,2)]
	bad1 <- which(S[,1] < 0 & S[,2] == 0)
	bad2 <- which(S[,1] < 0 & S[,2] == pi)
	bad <- c(bad1, bad2)
	S <- S[-bad,]
	bad <-  c()
	for(i in 1:dim(V)[1])
	{
		bad <- c(bad, which(S[,1] == V[i,1] & S[,2] == V[i,2]))
	}
	Splus <- S[-bad,]
	
	bad1 <- which(Splus[,1] == 0 | Splus[,2] == 0)
	bad2 <- which(Splus[,1] == pi | Splus[,2] == pi)
	Splus <- Splus[-c(bad1, bad2),]
	
	if( dim(Splus)[1]/2 != floor((nrows-1)/2)*floor((ncols-1)/2))
	{stop("error: wrong Splus values")}
	
	good <- c()
	for(i in 1:dim(Splus)[1] )
	{
		good <- c(good, which(Splus[i,1] == pdata[,1] & Splus[i,2] == pdata[,2]) )
	}

	Pj <- pdata[good,]
	
	pos <- which(Pj[,1] > 0 & Pj[,2] > 0)
	Pjpos <- Pj[pos,]
	Rj <- c()
	for(i in 1:dim(Pjpos)[1])
	{
		om1 <- Pjpos[i,1]
		om2 <- Pjpos[i,2]
		matchom <- which(Pj[,1] == -om1 & Pj[,2] == om2)
		rjnew <- Pjpos[i,3] / Pj[matchom,3]
		Rj <- c(Rj, rjnew)
	}
	
	Nr <- length(Rj)
	Rj <- sort(Rj)
	jseq <- 1:Nr
	Frj <- Rj/(1+Rj)
	ins <- ((2*jseq - 1)/(2*Nr) - Frj )^2
	CvM <- 1/(12*Nr) + sum(ins)

	CvMvec <- c()
	for(i in 1:nsim)
	{
		mysamp <- rf(Nr, 2, 2)
		mysamp <- sort(mysamp)
		Fx <- mysamp/(1+mysamp)
		ins <- ((2*jseq - 1)/(2*Nr) - Fx )^2
		cvm.sim <- 1/(12*Nr) + sum(ins)
		CvMvec <- c(CvMvec, cvm.sim) 
	}

	sum(CvMvec > CvM)/nsim
}
#' @keywords internal
test_complete_sym = function(pdata, nsim = 5000)
{
	V <- rbind(c(0,0), c(0, pi), c(pi,0), c(pi, pi))
	bad <- which( pdata[,2] < 0)
	Sstar <- pdata[-bad,1:2]
	bad <- c()
	for(i in 1:dim(V)[1])
	{
		bad <- c(bad, which(Sstar[,1] == V[i,1] & Sstar[,2] == V[i,2]))
	}
	Sstar <- Sstar[-bad,]
	bad <- which(Sstar[,1] == Sstar[,2])
	if(length(bad) > 0)
	{
		Sstar <- Sstar[-bad,]
	}
	bad <- which(Sstar[,1] < 0 & (Sstar[,1] == -1*Sstar[,2]) )
	if(length(bad > 0))
	{
			Sstar <- Sstar[-bad,]
	}
	bad <- which(Sstar[,1] < 0 & Sstar[,2] == 0)
	if(length(bad)>0)
	{
			Sstar <- Sstar[-bad,]
	}
	
	if(dim(Sstar)[1] == 0)
	{stop("No frequencies in S*, cannot do test of complete symmetry, test stopped")}
	
	good <-  c()
	for(i in 1:dim(Sstar)[1] )
	{
		good <-  c(good, which(pdata[,1] == Sstar[i,1] & pdata[,2] == Sstar[i,2] ) )
	}

	Pj <- pdata[good,]
	
	good1 <-  which(Pj[,1] >= 0 & Pj[,1] > Pj[,2])
	good2 <- which(Pj[,1] < 0 & Pj[,2] < -Pj[,1])
	Pj.loop <-  Pj[c(good1, good2),]
	Rj <- c()
	for(i in 1:dim(Pj.loop)[1])
	{
		om1 <- Pj.loop[i,1]
		om2 <- Pj.loop[i,2]
		if(om1 > 0)
		{
			matchom <- which(Pj[,1] == om2 & Pj[,2] == om1)
			rjnew <- Pj.loop[i,3] / Pj[matchom,3]
			Rj <- c(Rj, rjnew)
		}
		if(om1 < 0)
		{
			matchom <- which(Pj[,1] == -1*om2 & Pj[,2] == -1*om1)
			rjnew <- Pj.loop[i,3] / Pj[matchom,3]
			Rj <- c(Rj, rjnew)
		}	
	}
	
	Nr <- length(Rj)
	if(Nr == 0)
	{
		print("Cannot perform test of complete symmetry -- no frequencies to compare")
		return(1)
	}
	Rj <- sort(Rj)
	jseq <- 1:Nr
	Frj <- Rj/(1+Rj)
	ins <- ((2*jseq - 1)/(2*Nr) - Frj )^2
	CvM <- 1/(12*Nr) + sum(ins)

	CvMvec <- c()
	for(i in 1:nsim)
	{
		mysamp <- rf(Nr, 2, 2)
		mysamp <- sort(mysamp)
		Fx <- mysamp/(1+mysamp)
		ins <- ((2*jseq - 1)/(2*Nr) - Fx )^2
		cvm.sim <- 1/(12*Nr) + sum(ins)
		CvMvec <- c(CvMvec, cvm.sim) 
	}

	pvalue <- sum(CvMvec > CvM)/nsim
	return(pvalue)
}