#' @keywords internal
make_htestIso_MS <- function(HTout, df)
{
	nv <- 0
	names(nv) <- "difference in directional covariograms"

	mya <- "two-sided"
	names(mya) <- "difference in directional covariogram"

	mmeth <- "Test of isotropy from Maity and Sherman (2012) for sampling locations with general sampling design using the sample covariogram."
	names(mmeth) <- "method"

	myest <- c(HTout$C.hat[,3])
	nlags <- length(c(myest))
	tmp <- character()
	for(i in 1:nlags)
	{
		tmp <- c(tmp, paste("(",as.character(HTout$C.hat[i,1]),",",as.character(HTout$C.hat[i,2]) ,")", sep = ""))
	}
	names(myest) <- tmp

	myts <- HTout$test.stat
	names(myts) <- "Chi-sq"

	mparms <- df
	names(mparms) <- "df"

	mpv1 <- HTout$pvalue.chisq
	names(mpv1) <- "p.value.chisq"

	msh <- HTout$V.hat

	obj <- list(null.value = nv, alternative = "two.sided",
method = mmeth, estimate = myest,
statistic = myts, parameter = mparms, p.value = mpv1, sigma.hat = msh, p.value.refl = NULL, p.value.comp = NULL) 

		class(obj) <- c("htestIso")		
		return(obj)
}

#' @keywords internal
lag_dist_prod = function(spdata)
{
	locs <- spdata[,1:2]
	z <- spdata[,3]
	z <- z-mean(z)
	n <- dim(spdata)[1]
	splags <- c()
	prod <- c()
	for(i in 1:n)
	{
		lags.d <- cbind(locs[,1] - locs[i,1], locs[,2] - locs[i,2])
		splags <- rbind(splags, lags.d)
		prod <- c(prod, z[i]*z)
	}
	
	rv <- cbind(splags, prod)
	colnames(rv) <- c("x.lag", "y.lag", "zij")
	return(rv)
}
#' @keywords internal
epkern = function(u)
{
	rv <- rep(0, length(u))
	good <- which(u < 1 & u > -1)
	rv[good] <- 0.75*(1 - u[good]^2)
	return(rv)	
}
#' @keywords internal
est_chat_MS = function(raw_cov_data,lagmat, kappa = 1, user.bandwidth = T, bandwidth = c(1,1))
{
	chat <- apply(lagmat, 1, est_chat_t, raw_cov_data, kappa, user.bandwidth, bandwidth)
	chat.mat <- cbind(lagmat, chat)
	colnames(chat.mat) <- c("lag.x", "lag.y", "C.hat")
	return(chat.mat)
}
#' @keywords internal
est_chat_t = function(lag, raw_cov_data, kappa, usr.band, bw)
{
	lag.x <- lag[1]
	lag.y <- lag[2]
	
	Dx <- raw_cov_data[,1]
	Dy <- raw_cov_data[,2]
	
	if(usr.band == T)
	{
		hx <- bw[1]
		hy <- bw[2]
	}
	
	if(usr.band == F)
	{	
		sDx <- sd(Dx)
		sDy <- sd(Dy)
		N <-  dim(raw_cov_data)[1]
		hx <- kappa*sDx*N^(-1/5)
		hy <- kappa*sDy*N^(-1/5)
	}
	
	xarg <- (lag.x - Dx)/hx
	yarg <- (lag.y - Dy)/hy
	
	top <- epkern(xarg)*epkern(yarg)*raw_cov_data[,3]
	bot <- epkern(xarg)*epkern(yarg)
	
	chat <- sum(top)/sum(bot)
	return(chat)
}
#' @keywords internal
get_block_coords = function(blk.dims, xlims, ylims, grid = c(1,1))
{	
	x.grid <- seq(xlims[1], xlims[2], by = grid[1])
	y.grid <- seq(ylims[1], ylims[2], by = grid[2])

	bad.x <- which(x.grid + blk.dims[1] > xlims[2])
	bad.y <- which(y.grid + blk.dims[2] > ylims[2])
	if(length(bad.x) > 0)
	{
		x.grid <- x.grid[-bad.x]
	}
	if(length(bad.y) > 0)
	{
		y.grid <- y.grid[-bad.y]
	}

	Ln <- expand.grid(y.grid, x.grid)
	Ln <- cbind(Ln[,2], Ln[,1])

	kx.grid <- seq(xlims[1], xlims[2], by = blk.dims[1])
	ky.grid <- seq(ylims[1], ylims[2], by = blk.dims[2])
	bad.kx <- which( (kx.grid+blk.dims[1]) > xlims[2])
	bad.ky <- which( (ky.grid+blk.dims[2]) > ylims[2])
	if(length(bad.kx)>0)
	{
		kx.grid <- kx.grid[-bad.kx]
	}
	if(length(bad.ky)>0)
	{
		ky.grid <- ky.grid[-bad.ky]
	}

	Kn <- expand.grid(ky.grid, kx.grid)
	Kn <- cbind(Kn[,2], Kn[,1])
	
	rv <- list("Ln" = Ln, "Kn" = Kn)
}
#' @keywords internal
block_data_list = function(spdata, blk.dims, Ln.blk.coords)
{
	blk.data.list <- list()
	nblks <- dim(Ln.blk.coords)[1]
	for(i in 1:nblks)
	{
		xl <- Ln.blk.coords[i,1]
		yl <- Ln.blk.coords[i,2]
		xu <- xl + blk.dims[1]
		yu <- yl + blk.dims[2]
		inblk <- which(spdata[,1] >= xl & spdata[,1] < xu & spdata[,2] >= yl & spdata[,2] < yu)
		blk.data <- spdata[inblk,]	
		blk.data.list[[i]] <- blk.data
	}
	
	return(blk.data.list)
}
#' @keywords internal
spatial_boot = function(Kn, Ln, blk.data.list)
{
	nK.blks <- dim(Kn)[1]
	nL.blks <- dim(Ln)[1]
	nL.blk.indeces <- 1:nL.blks
	blk.dl.lengths <- lapply(blk.data.list, dim)
	blk.dl.length <- do.call(rbind, blk.dl.lengths)
	bad <- which(blk.dl.length[,1] == 0)
	if(length(bad) > 0)
	{
		nL.blk.indeces <- nL.blk.indeces[-bad]
		if(length(bad) > 1)
		{warning(paste("There are ", length(bad), " empty subblocks; consider increasing block size if there are many empty subblocks."))}
	}
	bsample <- sample(nL.blk.indeces, nK.blks, replace = T)
	
	sp.boot.data = matrix(data = NA, nrow = 0, ncol = 3)
	for(i in 1:nK.blks)
	{
		sp.boot.data <- rbind(sp.boot.data, change_coords(Kn[i,], Ln[ bsample[i] , ], blk.data.list[[ bsample[i] ]]) )
	}

	return(sp.boot.data)
}
#' @keywords internal
change_coords = function(blk1.coords, blk2.coords, blk2.data)
{
	xdiff <- blk1.coords[1] - blk2.coords[1]
	ydiff <- blk1.coords[2] - blk2.coords[2]
	blk2.data <- matrix(blk2.data, ncol = 3, byrow = F)
	shifted.blk.data <- cbind(blk2.data[,1]+xdiff,blk2.data[,2]+ydiff, blk2.data[,3])
	return(shifted.blk.data)
}
#' @keywords internal
est_block_chats = function(lags, spdata, nBoot, blk.dims, xlims, ylims, grid, kappa = 1, usr.band = F, bw = c(1,1))
{
	blk.dims <- c(grid[1]*blk.dims[1], grid[2]*blk.dims[2])
	blk.coords <- get_block_coords(blk.dims, xlims, ylims, grid)
	Kn <- blk.coords$Kn
	Ln <- blk.coords$Ln
	blk.dl <- block_data_list(spdata, blk.dims, blk.coords$Ln)
	
	boot.data.list <- list()
	for(i in 1:nBoot)
	{
		boot.data.list[[i]] <- spatial_boot(Kn, Ln, blk.dl)
	}
	
	rd.list <- lapply(boot.data.list, lag_dist_prod)
	boot.chat.list <- lapply(rd.list, est_chat_MS, lags, kappa, usr.band, bw)
	
	boot.Chats <- matrix(data = NA, nrow = nBoot, ncol = dim(lags)[1])
	for(i in 1:nBoot)
	{
		boot.Chats[i,] <- boot.chat.list[[i]][,3]
	}
	
	return(boot.Chats)	
}