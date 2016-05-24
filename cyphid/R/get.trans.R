get.trans <-
function(dat){
		lambda <- .9 
		norder <- 6
		samples <- seq(1,length(dat), length=length(dat))
		nbasis <- length(samples) + norder-2
		basis <- create.bspline.basis(c(1,length(dat)), nbasis, norder, samples)
		fdPar <- fdPar(basis, 4, lambda)
		fd <- smooth.basis(samples, dat, fdPar)$fd
		accfd <- deriv.fd(fd,2)
		acc <- eval.fd(samples, accfd)
		frames <- seq(1,length(dat))
    minindex <- get.peaks(-acc, 3)
		minframes <- frames[minindex]
		mindat <- acc[minindex]
		transmin <- which.min(mindat)
		return <- list(transmin=transmin, minframes=minframes, acc=acc)
	}

