flux.conv <-
function(fl.dat, ghg = "CH4", r2.qual = 0.8, nrmse.lim = 0.2, out.unit = "auto", elementar = FALSE, hardflag = list(range = TRUE)){
	lm4flux <- fl.dat$lm4flux
	dat <- fl.dat$orig.dat
	ct <- as.numeric(coef(lm4flux)[1] + coef(lm4flux)[2]*max(lm4flux$model[2]))
	c0 <- as.numeric(coef(lm4flux)[1])
	## set molar weight according to ghg
	m <- switch(ghg, CH4 = 16.04, N2O = 44.01, CO2 = 44.01)
	## average temperature in chamber during measurement
	T <- mean(dat$t.air, na.rm=TRUE)
	V <- dat$volume[1]
	A <- dat$area[1]
	t <- max(lm4flux$model[2])/60
	p <- mean(dat$p.air, na.rm=TRUE)
	## extract range limit. if no range limit is handed over, rl = NULL
	rl <- mean(dat$rl)
	## calculating the flux using the ideal gas law plus relations to chamber area and volume
	## details in gflux()!
	flux <- gflux(ct=ct, c0=c0, T=T, V=V, A=A, M=m, t=t, p=p)
	## according to the ghg the flux is transferred to the
	## common base g/m2*h
	flux <- switch(ghg, CO2 = flux/1e+6, CH4 = flux/1e+9, N2O = flux/1e+9)
	## according to the output.unit the unit is changed
	## per default the function tries to guess a unit that best reflects the actual value
	fluxes <- unlist(list(non = flux*0, ng = flux*1e+9, mug = flux*1e+6, mg = flux*1e+3, g = flux, t = flux/1e+3))
	if(out.unit == "auto"){
		if(flux==0){
			out.unit <- "non"
		}
		else{
			flux <- fluxes[(abs(fluxes) < 10) & (abs(fluxes) >= 0.01)]
			out.unit <- names(flux)
		}
	}
	else{
		flux <- fluxes[out.unit]
	}
	## provide easy access to model r2.adj
	r2 <- summary(lm4flux)$r.squared
	## provide easy access to model range
	range.m <- diff(range(lm4flux$model[,1], na.rm=TRUE))
	## provide easy access to model nrmse
	nrmse <- sqrt(sum(residuals(lm4flux)^2)/summary(lm4flux)$df[2])/diff(range(lm4flux$model[1], na.rm=TRUE))
	## set the r2-quality flag and hardflag flux (set to NA) if set so
	if(is.null(hardflag)) { hardflag <- list(dummy = 1) }
	r2.f <- ifelse(r2 >= r2.qual, TRUE, FALSE)
	if(exists("r2", hardflag)){
		if(hardflag$r2){
			flux <- ifelse(r2.f, flux, NA)
		}
	}
	## set the range-quality flag and hardflag (set to 0) flux if set so
	range.lim <- rl
	if(is.na(range.lim)){range.lim <- 0}
	range.f <- ifelse(range.m >= range.lim, TRUE, FALSE)
	if(exists("range", hardflag)){
		if(hardflag$range){
			flux <- ifelse(range.f, flux, 0)
		}
	}
	## set the nrmse-quality flag
	nrmse.f <- ifelse(nrmse <= nrmse.lim, TRUE, FALSE)
	if(exists("nrmse", hardflag)){
		if(hardflag$nrmse){
			flux <- ifelse(nrmse.f, flux, NA)
		}
	}
	## check nomba
	## ambient values from Mace Head Ireland and global average (CO2)
	## via http://cdiac.ornl.gov/pns/current_ghg.html as of August 16th, 2013
	ambient <- switch(ghg, CH4 = 1874, N2O = 324, CO2 = 392.6)
	nomba <- sum(dat[,1] <= ambient)
	if(exists("nomba", hardflag)){
		flux <- ifelse(nomba > hardflag$nomba, NA, flux)
	}
	## give back as element value
	if(elementar){
		flux <- switch(ghg, CO2 = flux*12/44, CH4 = flux*12/16, N2O = flux*28/44)
		ghg <- paste(ghg, switch(ghg, CO2 = "C", CH4 = "C", N2O = "N"), sep="-")
	}
	## prepare output
	fluss <- list(ghg = ghg, flux = flux, r2 = r2, nrmse = nrmse, r2.f = r2.f, range.f = range.f, nrmse.f = nrmse.f, nomba.f = nomba, leak.f = NA)
	## changing names to reflect on hardflags (somewhat ugly)
	hf <- pmatch(names(hardflag)[sapply(hardflag, is.logical)], names(fluss)[5:9])
	nms <- sapply(strsplit(names(fluss), split="\\."), function(x) x[1])
	nms[5:9][hf] <- paste(nms[5:9][hf], ".hf", sep="")
	nms[5:9][-hf] <- paste(nms[5:9][-hf], ".f", sep="")
	names(fluss) <- nms
	## put output together
	res <- list(fluss = fluss, fl.dat = fl.dat, unit = out.unit)
	class(res) <- "flux"
	return(res)
	}

