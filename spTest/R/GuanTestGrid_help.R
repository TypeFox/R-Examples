#' @keywords internal
make_htestIso_GG <- function(HTout, df)
{
	nv <- 0
	names(nv) <- "difference in directional semivariograms"

	mya <- "two-sided"
	names(mya) <- "difference in directional semivariogram"

	mmeth <- "Test of isotropy from Guan et. al. (2004) for gridded sampling locations using the sample semivariogram."
	names(mmeth) <- "method"

	myest <- c(HTout$gamma.hat[,3])
	nlags <- length(c(myest))
	tmp <- character()
	for(i in 1:nlags)
	{
		tmp <- c(tmp, paste("(",as.character(HTout$gamma.hat[i,1]),",",as.character(HTout$gamma.hat[i,2]) ,")", sep = ""))
	}
	names(myest) <- tmp

	myts <- HTout$test.stat
	names(myts) <- "Chi-sq"

	mparms <- df
	names(mparms) <- "df"

	mpv1 <- HTout$pvalue.chisq
	names(mpv1) <- "p.value.chisq"

	mpv2 <- HTout$pvalue.finite
	names(mpv2) <- "p.value.finite"

	msh <- HTout$sigma.hat

	mblk <- HTout$n.subblocks
	names(mblk) <- "No. of Subblocks"

	obj <- list(null.value = nv, alternative = "two.sided",
method = mmeth, estimate = myest,
statistic = myts, parameter = mparms, p.value = mpv1, p.value.finite = mpv2, sigma.hat = msh, n.subblocks = mblk, p.value.refl = NULL, p.value.comp = NULL) 

		class(obj) <- c("htestIso")		
		return(obj)
}

#' @keywords internal
scale_coords_guan = function(spdata, delta)
{
	spdata.orig <- spdata
	min.x <-  min(spdata[,1])
	min.y <-  min(spdata[,2])
	new.xcoords <- (spdata[,1] - min.x)/delta
	new.ycoords <- (spdata[,2] - min.y)/delta
	t.coords <- cbind(new.xcoords, new.ycoords)
	spdata[,1:2] <- t.coords
	return(spdata)
}
#' @keywords internal
lag_dist_diff_reg = function(spdata)
{
	locs <-  spdata[,1:2]
	z <-  spdata[,3]
	z <-  z-mean(z)
	n <-  dim(spdata)[1]
	index <-  sort(rep(1:n,n))
	splags <-  c()
	prod <-  c()
	origin.pts <-  matrix(data = NA, nrow = 0, ncol = 2)
	for(i in 1:n)
	{
		cpt.mat <-  cbind(rep(locs[i,1], n), rep(locs[i,2], n))
		origin.pts <-  rbind(origin.pts, cpt.mat)
		lags.d <-  -1*cbind(locs[i,1] - locs[,1], locs[i,2] - locs[,2])
		splags <-  rbind(splags, lags.d)
		prod <-  c(prod, (z[i] - z)^2)
	}
	x.coord <- origin.pts[,1]
	y.coord <- origin.pts[,2]
	x.lag <- splags[,1]
	y.lag <- splags[,2]
	xts <- prod
	rv <- cbind(x.coord, y.coord, index, x.lag, y.lag, xts)
	row.names(rv) <- NULL
	return(rv)
}
#' @keywords internal
raw_data_sub = function(rawdata, lagmat)
{
	nlags <-  dim(lagmat)[1]
	good <-  c()
	for(i in 1:nlags)
	{
		gd <- which(rawdata[,4] == lagmat[i,1] & rawdata[,5] == lagmat[i,2])
		good <- c(good, gd)
	}
	
	return(rawdata[good,])
}
#' @keywords internal
est_gamma2 = function(rawdata, lagmat, edge.eff = T)
{
	nlags <-  dim(lagmat)[1]
	n.bin <-  c()
	gamma.hat <-  c()
	
	if(edge.eff == F)
	{
		for(i in 1:nlags)
		{
			clag <-  lagmat[i,]
			good <-  which(rawdata[,4] == clag[1] & rawdata[,5] == clag[2])
			n.clag <-  length(good)
			n.bin <-  c(n.bin, n.clag)
			gh <-  sum(rawdata[good,6])/(2*n.clag)
			gamma.hat <-  c(gamma.hat, gh)
		}
		
		gamma.est <-  cbind(lagmat, gamma.hat, n.bin)
		colnames(gamma.est) <-  c("lag.x","lag.y","gamma.hat", "n.bin")
		rv <- list("gamma.hat" = gamma.est)
		return(rv)
	}

	if(edge.eff == T)
	{
		index.set <-  c()
		clag <-  lagmat[1,]
		good <-  which(rawdata[,4] == clag[1] & rawdata[,5] == clag[2])
		index.set <-  rawdata[good,3]
		
		for(i in 2:nlags)
		{
			clag <-  lagmat[i,]
			good <-  which(rawdata[,4] == clag[1] & rawdata[,5] == clag[2])
			index.tmp <-  rawdata[good,3]
			index.set <-  intersect(index.set, index.tmp)
		}	
		
		good <-  which(rawdata[,3] %in% index.set)
		rawdata <-  rawdata[good,]
		for(i in 1:nlags)
		{
			clag <-  lagmat[i,]
			good <-  which(rawdata[,4] == clag[1] & rawdata[,5] == clag[2])
			n.clag <-  length(good)
			n.bin <-  c(n.bin, n.clag)
			gh <-  sum(rawdata[good,6])/(2*n.clag)
			gamma.hat <-  c(gamma.hat, gh)
		}

		good <-  which(rawdata[,3] %in% index.set)
		gamma.est <-  cbind(lagmat, gamma.hat, n.bin)
		colnames(gamma.est) <-  c("lag.x","lag.y","gamma.hat", "n.bin")
		rv <- list("gamma.hat" = gamma.est, "edge.data" = rawdata[good,])
		return(rv)
	}
}
#' @keywords internal
est_sigma_reg = function(spdata, lagmat, blk.width, blk.height, edge.eff = T, finite.sample = T)
{	
	npts <- dim(spdata)[1]
	sb.card <- blk.width*blk.height
	fn <- 1 - (sb.card/npts)
	
	n.lags <- dim(lagmat)[1]
	
	if(edge.eff == T)
	{
		rawdata <- lag_dist_diff_reg(spdata)
		rawdata <- raw_data_sub(rawdata, lagmat)
		gh.data <- est_gamma2(rawdata, lagmat, edge.eff = T)
		ea.data <- gh.data[[2]]
		blk.data.list <- block_data_xts(ea.data, blk.width, blk.height)
		block.ghat.list <- lapply(blk.data.list, ghat_block_ea, lagmat)
		block.ghats <- do.call(rbind, block.ghat.list)
	}
	
	if(edge.eff == F)
	{
		blk.data.list <- block_data_zt(spdata, blk.width, blk.height)
		blk.data.diff.list <- lapply(blk.data.list, lag_dist_diff_reg)
		blk.data.diff.list <- lapply(blk.data.diff.list, raw_data_sub, lagmat)
		block.ghats.list <- lapply(blk.data.diff.list, est_gamma2, lagmat, edge.eff = F)
		block.ghats <- do.call(rbind, block.ghats.list)
	}
	
	ghat.mean <-  apply(block.ghats, 2, mean)
	block.ghats.mc <-  block.ghats
	for(i in 1:dim(block.ghats)[1])
	{
		block.ghats.mc[i,] <-  block.ghats[i,] - ghat.mean
	}
	
	if(finite.sample == T)
	{
		kn <- dim(block.ghats)[1]
		sigma.hat <-  sb.card/(kn*fn)*t(block.ghats.mc)%*%block.ghats.mc
	}
	
	if(finite.sample == F)
	{
		sigma.hat <- cov(block.ghats)
	}
	
	return(list("sigma.hat" = sigma.hat, "block.ghats" = block.ghats))
}
#' @keywords internal
block_data_xts = function(eadata, bw, bh)
{
	min.x <- min(eadata[,1])
	max.x <- max(eadata[,1])
	min.y <- min(eadata[,2])
	max.y <- max(eadata[,2])
	
	blk.x <- min.x:(max.x-bw+1)
	blk.y <- min.y:(max.y-bh+1)
	
	blk.coords <- expand.grid(blk.x, blk.y)
	blk.coords <- cbind(blk.coords[,1], blk.coords[,2])

	block.data <-  list()
	n.blks <- dim(blk.coords)[1]
	for(i in 1:n.blks)
	{
		cblock <- blk.coords[i,]
		inb <- which(eadata[,1] >= cblock[1] & eadata[,1] < cblock[1]+bw & eadata[,2] >= cblock[2] & eadata[,2] < cblock[2]+bh)
		block.data[[i]] <- eadata[inb,]
	}
	
	return(block.data)
}
#' @keywords internal
block_data_zt = function(spdata, bw, bh)
{
	min.x <- min(spdata[,1])
	max.x <- max(spdata[,1])
	min.y <- min(spdata[,2])
	max.y <- max(spdata[,2])
	
	blk.x <- min.x:(max.x-bw+1)
	blk.y <- min.y:(max.y-bh+1)
	
	blk.coords <- expand.grid(blk.x, blk.y)
	blk.coords <- cbind(blk.coords[,1], blk.coords[,2])
	
	block.data <-  list()
	n.blks <- dim(blk.coords)[1]
	for(i in 1:n.blks)
	{
		cblock <- blk.coords[i,]
		inb <- which(spdata[,1] >= cblock[1] & spdata[,1] < cblock[1]+bw & spdata[,2] >= cblock[2] & spdata[,2] < cblock[2]+bh)
		block.data[[i]] <- spdata[inb,]
	}
	
	return(block.data)
}
#' @keywords internal
ghat_block_ea = function(blk.data, lagmat)
{
	nlags <- dim(lagmat)[1]
	n.bin <- c()
	gamma.hat <- c()
	
	for(i in 1:nlags)
	{
		clag <-  lagmat[i,]
		good <-  which(blk.data[,4] == clag[1] & blk.data[,5] == clag[2])
		n.clag <-  length(good)
		n.bin <-  c(n.bin, n.clag)
		gh <-  sum(blk.data[good,6])/(2*n.clag)
		gamma.hat <-  c(gamma.hat, gh)
	}
	
	return(gamma.hat)
}
