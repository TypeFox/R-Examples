AFactMod <- function(dat, demean = TRUE
			, add.effects = c("none", "individual", "time", "twoways")
			, dim.criterion = c("PC1", "PC2", "PC3", "IC1", "IC2" 
			, "IC3", "IPC1", "IPC2", "IPC3", "ED", "ER", "GR")
			, factor.dim, d.max, sig2.hat
			, restrict.mode= c("restrict.factors","restrict.loadings")
			, allow.dual = TRUE){

 	is.regular.panel(dat, stopper = TRUE) 
      nr   <- nrow(dat)
	nc   <- ncol(dat)
  	with.trans <- match.arg(add.effects)
	dat.trans <- FUN.with.trans(dat, N = nc, T = nr, is.intercept = demean, effect = with.trans)
	dat <- dat.trans$TDM
  ## missing parameters

	if(missing(factor.dim)) factor.dim  <- NULL
	if(missing(d.max))      d.max 	<- NULL
	if(missing(sig2.hat))   sig2.hat    <- NULL

  ## pca.fit

	obj.pca.fit <- pca.fit(dat, given.d       = factor.dim
					  , restrict.mode = restrict.mode
					  , allow.dual    = allow.dual)

  ## dimension selection
	dim.criterion <- match.arg(dim.criterion)
	est.dim       <- switch(dim.criterion,
				PC1 = B.OptDim(obj.pca.fit, criteria = c("PC1")
					, d.max = d.max, sig2.hat = sig2.hat),
				PC2 = B.OptDim(obj.pca.fit, criteria = c("PC2")
					, d.max = d.max, sig2.hat = sig2.hat),
				PC3 = B.OptDim(obj.pca.fit, criteria = c("PC3")
					, d.max = d.max, sig2.hat = sig2.hat),
				IC1 = B.OptDim(obj.pca.fit, criteria = c("IC1")
					, d.max = d.max, sig2.hat = sig2.hat),
				IC2 = B.OptDim(obj.pca.fit, criteria = c("IC2")
					, d.max = d.max, sig2.hat = sig2.hat),
				IC3 = B.OptDim(obj.pca.fit, criteria = c("IC3")
					, d.max = d.max, sig2.hat = sig2.hat),
				IPC1 = B.OptDim(obj.pca.fit, criteria = c("IPC1")
					, d.max = d.max, sig2.hat = sig2.hat),
				IPC2 = B.OptDim(obj.pca.fit, criteria = c("IPC2")
					, d.max = d.max, sig2.hat = sig2.hat),
				IPC3 = B.OptDim(obj.pca.fit, criteria = c("IPC3")
					, d.max = d.max, sig2.hat = sig2.hat),
				ED   = O.OptDim(obj.pca.fit, d.max = d.max),
				
				ER   = RH.OptDim(obj.pca.fit, criteria = c("ER")
					, d.max = d.max),
				GR   = RH.OptDim(obj.pca.fit, criteria = c("GR")
					, d.max = d.max),
                                )
	opt.d  <- est.dim[1,2]
	used.d <- ifelse(is.null(factor.dim), opt.d, factor.dim)

  ## factors and loadings parameters

	if(used.d!=0){
		factors <- obj.pca.fit$factors[, 1:used.d, drop= FALSE]
		loadings  <- obj.pca.fit$loadings[, 1:used.d, drop= FALSE]
		dat.fit <- tcrossprod(factors, loadings)
		sd2 <- obj.pca.fit$Sd2[used.d+1]
		R <- list(fitted.values =dat.fit, factors = factors
			, loadings = loadings, sd2= sd2, given.fdim = factor.dim
			, optimal.fdim = opt.d, used.fdim = used.d)
		}
	else{
		factors <- matrix(0, nr, 1)
		scores  <- matrix(0, nc, 1)
		dat.fit <- tcrossprod(factors, scores)
		sd2 <- obj.pca.fit$Sd2[used.d+1]
		R <- list(fitted.values =dat.fit, factors = factors
			, loadings = loadings, resid.sd2= sd2, given.fdim = factor.dim
			, optimal.fdim = opt.d, used.fdim = used.d)
		}
	R
      }
