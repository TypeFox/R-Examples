mf.flux <- 
function(x, var.par, method = "r2", time.unit = "S", all.through = TRUE, iv = 1, wndw = 0.1, pdk = 0.5, min.dp = 20, nrmse.lim = 0.1, r2.qual = 0.9, range.lim = 5, out.unit = "auto", elementar = FALSE, hardflag = list(range = TRUE), consecutive = FALSE){
	## which method to use
	METHODS <- c("r2", "rmse", "AIC")
	method <- pmatch(method, METHODS)
	
	## prepare input
	n <- nrow(x)
	# function vp that realises the compilation of dat for further use
	vp <- function(x, sel){
		if(is.character(sel)){out <- x[,sel]}
		else{out <- rep(sel, n)}
		return(out)
	}
	# using vp to extract variables from x and combine with fixed parameters
	dat <- data.frame(sapply(var.par, function(y) vp(x, y)))
	
	## prepare and organize handthrough via all.through
	if(all.through){
		x.out <- x[ , -match(unlist(var.par)[c("ghg", "time")], names(x))]
		do.means <- vector(length = ncol(x.out))
		for(i in c(1:ncol(x.out))){
			do.means[i] <- is.numeric(x.out[,i])
		}	
		dat.out <- lapply(which(!do.means), function(x) as.character(x.out[1,x]))
		dat.out <- c(dat.out, lapply(which(do.means), function(x) mean(x.out[,x], na.rm=TRUE)))
		names(dat.out) <- c(names(x.out)[!do.means], names(x.out)[do.means])
	}
	
	## prepare and organize handthrough via var.par (ATTENTION there seems to be a problem here)
	else{
		stv <- match(c("ghg", "time", "area", "volume"), names(dat))
		handthrough <- names(which(sapply(var.par[-stv], is.character)))
		dat.out <- lapply(c("area", "volume", handthrough), function(y) mean(dat[,y], na.rm=TRUE))
		names(dat.out) <- paste("htd", names(dat.out), sep=".")
	}
	
	## extract variables from x(dat) for flux calculation
	x <- dat[,c("ghg", "time", "area", "volume", "t.air", "p.air")]
	# prepare time entries
	x$time <- as.vector(x$time - min(x$time))
	x <- x[order(x$time),]
	# if time follows a fixed interval, create times
	if(iv!=1){
		x$time <- seq(length(x$time))*iv
	}
	
	## other preparations
	# put original data aside for later reporting
	inn <- x
	# get original ghg.range
	ghg.range <- range(x$ghg, na.rm=TRUE)
	
	## consecutive approach
	if(consecutive){
		ghg.n <- x$ghg-min(x$ghg)+1
		ghg.n <- ghg.n/max(ghg.n)
		time.n <- seq(n)
		wndw.n <- n*wndw
		indizes <- embed(time.n, wndw.n)
		sel <- which(runsd(apply(indizes, 1, function(y) coef(lm(ghg.n[y] ~ time.n[y]))[2]), wndw.n) < 0.005)
		sel <- c(sel, max(sel)+seq(floor(wndw.n)))
		lm4flux <- lm(ghg ~ time, data=x[sel,])
	}
	
	## check for high frequency strong fluctuations and remove data points
	# first check for unbalanced concentration measurements
	# in this case hff is defined all which is not the prevailing value
	xghg.fac <- as.factor(x$ghg)
	unik <- summary(xghg.fac, maxsum=length(levels(xghg.fac)))
	if(max(unik)/length(x$ghg) > pdk){
		conz <- as.numeric(names(sort(unik, decreasing=TRUE)))[1]
		x$hff <- x$ghg!=conz
	}
	# if there is no prevailing value check hff according to running sd compared to range
	# but only when range.lim <= ghg.range
	else{
		if(diff(ghg.range) > range.lim){
			rsd <- runsd(x$ghg, ceiling(n*wndw))/diff(ghg.range)
			x$hff <- rsd >= wndw
		}
		else{
			x$hff <- FALSE
		}
	}
	x <- x[!x$hff,]
		
	## check range limit and do shortcut if range of concentrations is below range.lim
	if(hardflag$range & (diff(ghg.range) <= range.lim)){
		flux <- 0
		x$ghg <- rep(mean(x$ghg), length(x$ghg))
		lm4flux <- lm(ghg ~ time, data = x)
		if(out.unit == "auto"){ out.unit <- "g" }
		r2 <- 1
		range.m <- diff(ghg.range)
		nrmse <- 0
		n.inn <- length(inn[,1])
		n.out <- n.inn
		podpu <- 1
	}
	else{
		## getting the best linear part
		# derive minimum number of measurements
		n <- nrow(x)
		ndk <- ifelse((n*pdk) >= min.dp, ceiling(n*pdk), min.dp)
		# with moving windows of various lengths
		if(ndk >= n){  
			ndk <- n-1
			warning("Number of entries <= min.dp and min.dp has been automatically adjusted to n-1")
		}
		winds <- c(ndk:n)
		indizes <- lapply(winds, function(y) lapply(c(1:n), function(x) seq(x,(x+y)))[1:(n-y)])
		indizes <- unlist(indizes, recursive=FALSE)
		if(method==1){
			rs <- sapply(indizes, function(y) summary(lm(ghg ~ time, data=x[y,]))$r.squared)
			sel <- indizes[[order(rs, decreasing = TRUE)[1]]]
		}
		if(method==2){
			lms <- lapply(indizes, function(y) lm(ghg ~ time, data=x[y,]))
			rmse <- unlist(lapply(lms, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])))
			sel <- indizes[[order(rmse)[1]]]
		}
		if(method==3){
			AIC <- sapply(indizes, function(y) AIC(lm(ghg ~ time, data=x[y,])))
			sel <- indizes[[order(AIC)[1]]]
		}				
		# actually get the best model
		lm4flux <- lm(ghg ~ time, data = x[sel,])

		## flux calculation
		# prepare
		time.range <- range(lm4flux$model[,2])
		dc <- diff(predict(lm4flux, newdata = data.frame(time = time.range)))
		m <- 44.01
		T <- mean(x$t.air, na.rm=TRUE)
		V <- x$volume[1]
		A <- x$area[1]
		t <- diff(time.range)
		t <- switch(time.unit, S = t/60/60, M = t/60, H = t)
		p <- mean(x$p.air, na.rm=TRUE)
		# do calculation
		flux <- gflux(ct=dc, T=T, V=V, A=A, M=m, t=t, p=p)
		flux <- flux/1e+6

		## according to the output.unit the unit is changed
		## per default the function tries to guess a unit that best reflects the actual value
		fluxes <- unlist(list(fg = flux*1e+15, pg = flux*1e+12, ng = flux*1e+9, mug = flux*1e+6, mg = flux*1e+3, g = flux, kg = flux/1e+3))
		if(out.unit == "auto"){
			flux <- fluxes[(abs(fluxes) < 10) & (abs(fluxes) >= 0.01)]
			out.unit <- names(flux)
		}
		else{
			flux <- fluxes[out.unit]
		}

		## extract model parameters for flagging and output
		# extract model r2.adj, range, and nrmse
		r2 <- summary(lm4flux)$r.squared
		range.m <- diff(range(lm4flux$model[,1], na.rm=TRUE))
		if(range.m == 0){ range.m <- 0.001 }
		nrmse <- sqrt(sum(residuals(lm4flux)^2)/summary(lm4flux)$df[2])/range.m
		# extract model n, original n, and percentage of data points used
		n.out <- length(lm4flux$model[,1])
		n.inn <- length(inn[,1])
		podpu <- n.out/n.inn		
	}
	
	## flagging
	# set the r2-quality flag and hardflag flux (set to NA) if set so
	r2.f <- ifelse(r2 >= r2.qual, TRUE, FALSE)
	if(exists("r2", hardflag)){
		if(hardflag$r2){
			flux <- ifelse(r2.f, flux, NA)
		}
	}
	# set the range-quality flag and hardflag (set to 0) flux if set so
	if(is.na(range.lim)){range.lim <- 0}
	range.f <- ifelse(diff(ghg.range) >= range.lim, TRUE, FALSE)
	if(hardflag$range){
		flux <- ifelse(range.f, flux, 0)
	}
	# set the nrmse-quality flag
	nrmse.f <- ifelse(nrmse <= nrmse.lim, TRUE, FALSE)
	if(exists("nrmse", hardflag)){
		if(hardflag$nrmse){
			flux <- ifelse(nrmse.f, flux, NA)
		}
	}
	# check nomba
	# ambient values from Mace Head Ireland and global average (CO2)
	# via http://cdiac.ornl.gov/pns/current_ghg.html as of August 1st, 2011
	ghg <- var.par$ghg
	ambient <- switch(ghg, CH4 = 1870, N2O = 323, CO2 = 388.5)
	nomba <- sum(dat[,1] <= ambient)
	if(exists("nomba", hardflag)){
		flux <- ifelse(nomba > hardflag$nomba, NA, flux)
	}
	# give back as element value
	if(elementar){
		flux <- switch(ghg, CO2 = flux*12/44, CH4 = flux*12/16, N2O = flux*28/44)
		ghg <- paste(ghg, switch(ghg, CO2 = "C", CH4 = "C", N2O = "N"), sep="-")
	}
	
	## prepare output
	fluss <- list(ghg = ghg, flux = flux, r2 = r2, nrmse = nrmse, r2.f = r2.f, range.f = range.f, nrmse.f = nrmse.f, nomba.f = nomba, unit = out.unit, podpu = podpu)
	# changing names to reflect on hardflags (somehwat ugly)
	hf <- pmatch(names(hardflag)[sapply(hardflag, is.logical)], names(fluss)[5:8])
	nms <- sapply(strsplit(names(fluss), split="\\."), function(x) x[1])
	nms[5:8][hf] <- paste(nms[5:8][hf], ".hf", sep="")
	nms[5:8][-hf] <- paste(nms[5:8][-hf], ".f", sep="")
	names(fluss) <- nms
	# put output together
	res <- list(fluss = fluss, mod = lm4flux, out = dat.out, inn = inn)
	class(res) <- "fluxx"
	cat(".")
	return(res)
}