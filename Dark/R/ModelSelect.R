ModelSelect <- function(obj, P) {

	.GlobalEnv$x <- obj$time
	.GlobalEnv$y <- obj$thrs

	param <- matrix(0, 3, 7)
	idx <- c(3, 5, 7)
	AIC <- numeric(7)
	mFn <- c(1, 1, P3, 1, P5c, 1, P7c)
	jj = 1

	for (ii in idx) {
		Fn = mFn[[ii]]

		Out <- apply(P, 1, FUN = Fn)
		oP <- P[which(Out == min(Out)), ]
		oPval <- Out[which(Out == min(Out))]
		Opt <- optim(oP[1:ii], Fn)
		while (Opt$con) Opt <- optim(Opt$par, Fn)

		init <- oP[1:Fn(1)$Pn]
		opt <- Opt$par[1:Fn(1)$Pn]
		param[jj, 1:Fn(1)$Pn] <- opt
		val <- Opt$val
		Pn <- Fn(1)$Pn
		Mod <- Fn(1)$Mod
		# builds the dark object to pass to AICc
		obj$init = init
		obj$opt = opt
		obj$val = val
		obj$Pn = Pn
		obj$Mod = Mod
		AIC[ii] <- AICc(obj)

		jj = jj + 1
	}

	on.exit(rm(list = c("x", "y"), envir = .GlobalEnv))

	list(AIC = AIC, param = param)
}
