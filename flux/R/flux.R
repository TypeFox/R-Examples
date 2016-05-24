flux <-
function(x, var.par, co2ntrol = list(leak = TRUE, relay = FALSE), min.allowed = 3, max.nrmse = 0.1, nrmse.lim = 0.2, r2.qual = 0.8, range.lim = 30, out.unit = "auto", elementar = FALSE, hardflag = list(range = TRUE), asterisks = TRUE){
	## extract the name vector from x
	nmes <- x$nmes
	## extract the data tables from x
	x <- x$tables
	## provide a leak.flag in case no co2ntrol is wanted
	leak.flag <- rep(2, length(x))
	if(is.null(co2ntrol)){co2ntrol <- list(leak = FALSE, relay = FALSE)}
	## provide correct quality parameters for each gas
	if(length(max.nrmse)==1){max.nrmse <- list(CO2 = max.nrmse, CH4 = max.nrmse, N2O = max.nrmse)}
	if(length(nrmse.lim)==1){nrmse.lim <- list(CO2 = nrmse.lim, CH4 = nrmse.lim, N2O = nrmse.lim)}
	if(length(r2.qual)==1){r2.qual <- list(CO2 = r2.qual, CH4 = r2.qual, N2O = r2.qual)}
	if(!is.null(range.lim)){
		if(class(range.lim)=="numeric"){range.lim <- list(CO2 = range.lim, CH4 = range.lim, N2O = range.lim)}
		else{
			sel <- match(names(range.lim), c("CO2", "CH4", "N2O"))
			rl <- list(CO2 = 0, CH4 = 0, N2O = 0)
			rl[sel] <- range.lim
			range.lim <- rl
		}
	}
	## prepare out.unit handling
	if(length(out.unit)==1){out.unit <- list(CO2 = out.unit, CH4 = out.unit, N2O = out.unit)}
	## prepare output. necessary to have empty slots later
	## put the original results into a list
	flux.res <- list(CO2 = NA, CH4 = NA, N2O = NA)
	## prepare all in one big results table
	flux.table <- data.frame(nmes)
	## run the outlier detection and elimination routines for the ghg and estimate the fluxes
	## prepare selectors
	co <- grep("CO2", names(var.par))
	ch <- grep("CH4", names(var.par))
	no <- grep("N2O", names(var.par))
	vp <- var.par[-c(co, ch, no)]
	cat(".",sep="")
	## CO2 stuff
	if(length(co)!=0){
		var.par.CO2 <- var.par
		## check for gcq and add 0 if necessary
		if(sum(names(var.par.CO2)=="CO2.gcq")==0) {
			var.par.CO2$CO2.gcq = 0
			warning("CO2 GC quality flags have been set to zero")
		}
		# update co
		co <- grep("CO2", names(var.par.CO2))
		##
		names(var.par.CO2)[names(var.par.CO2)=="CO2"] <- "ghg"
		names(var.par.CO2)[names(var.par.CO2)=="CO2.gcq"] <- "gc.qual"
		## check about range.lim and attach calibration data to the data
		if(!is.null(range.lim)){
			CO2.rl <- vector("numeric", length(x))
			for(i in c(1:length(x))){
				CO2.rl[i] <- ifelse(length(range.lim$CO2)==1, range.lim$CO2, range.lim$CO2[i])
				x[[i]]$CO2.rl <- CO2.rl[i]
			}
		}
		else{
			CO2.rl <- sapply(x, function(y) y$CO2.rl[1])
		}	
		CO2.pre <- lapply(x, function(x) flux.odae(x, var.par = c(var.par.CO2[co], vp), min.allowed = min.allowed, max.nrmse = max.nrmse$CO2, rl="CO2.rl"))
		cat(".",sep="")
		## do prep for CO2 and co2ntrol if wanted
		if(!is.null(co2ntrol)){
			CO2.range.check <- sapply(x, function(x) ifelse(diff(range(x[,var.par.CO2$ghg])) >= mean(x$CO2.rl), TRUE, FALSE))
			if(co2ntrol$leak){
				leak.flag <- sapply(CO2.pre, function(x) coef(x$lm4flux)[2]) <= 0
				leak.flag <- !as.logical(leak.flag*CO2.range.check)
			}
			if(co2ntrol$relay){
				x.sub <- lapply(which(CO2.range.check), function(y) x[[y]][CO2.pre[[y]]$row.select,])
				x[which(CO2.range.check)] <- x.sub
			}
		}
		## run the CO2 flux estimation via flux.conv and gflux
		CO2.res <- lapply(CO2.pre, function(x) flux.conv(x, ghg = "CO2", r2.qual = r2.qual$CO2, nrmse.lim = nrmse.lim$CO2, out.unit = out.unit$CO2, elementar = elementar, hardflag = hardflag))
		if(co2ntrol$leak){
			for(i in c(1:length(CO2.res))){
				CO2.res[[i]]$fluss$leak.f <- leak.flag[i]
			}
		}		
		flux.res$CO2 <- CO2.res
		## when pv values are to be reported… extract them
		if(asterisks){
			CO2.pv <- sapply(CO2.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
			CO2.pv <- as.vector(symnum(CO2.pv, corr=FALSE, cutpoints = c(0,.001,.01,.05,.1,1), symbols = c("***","**","*","."," ")))
		}
		## make table for CO2
		CO2.table <- t(sapply(CO2.res, function(x) unlist(x$fluss[2:8])))
		CO2.units <- sapply(CO2.res, function(x) x$unit)
		CO2.table <- data.frame(CO2.units, CO2.pv, CO2.table)
		CO2.table$leak.f <- leak.flag*1
		CO2.table$rl <- CO2.rl
		names(CO2.table)[-c(1:2)] <- paste("CO2.", names(CO2.table)[-c(1:2)], sep="")
		flux.table <- data.frame(flux.table, CO2.table)
	}
	cat(".",sep="")
	## CH4 stuff
	if(length(ch)!=0){
		## do prep for CH4
		var.par.CH4 <- var.par
		## check for gcq and add 0 if necessary
		if(sum(names(var.par.CH4)=="CH4.gcq")==0) {
			var.par.CH4$CH4.gcq = 0
			warning("CH4 GC quality flags have been set to zero")
		}
		# update ch
		ch <- grep("CH4", names(var.par.CH4))
		##		
		names(var.par.CH4)[names(var.par.CH4)=="CH4"] <- "ghg"
		names(var.par.CH4)[names(var.par.CH4)=="CH4.gcq"] <- "gc.qual"
		## ATTENTION, check this: #names(var.par.CH4)[names(var.par.CH4)=="CH4.rl"] <- "rl"
		## check about range.lim and attach calibration data to the data
		if(!is.null(range.lim)){
			CH4.rl <- vector("numeric", length(x))
			for(i in c(1:length(x))){
				CH4.rl[i] <- ifelse(length(range.lim$CH4)==1, range.lim$CH4, range.lim$CH4[i])
				x[[i]]$CH4.rl <- CH4.rl[i]
			}
		}
		else{
			CH4.rl <- sapply(x, function(y) y$CH4.rl[1])
		}		
		CH4.pre <- lapply(x, function(x) flux.odae(x, var.par = c(var.par.CH4[ch], vp), min.allowed = min.allowed, max.nrmse = max.nrmse$CH4, rl="CH4.rl"))
		## run the CH4 flux estimation via flux.conv and gflux
		CH4.res <- lapply(CH4.pre, function(x) flux.conv(x, ghg = "CH4", r2.qual = r2.qual$CH4, nrmse.lim = nrmse.lim$CH4, out.unit = out.unit$CH4, elementar = elementar, hardflag = hardflag))
		if(co2ntrol$leak){
			for(i in c(1:length(CH4.res))){
				CH4.res[[i]]$fluss$leak.f <- leak.flag[i]
			}
		}
		if(exists("leak", hardflag)){
			if(hardflag$leak){
				if(length(co)==0){warning("leakage can only be evaluated and hardflagged if CO2 data are given")}
				else{
					for(i in c(1:length(CH4.res))){
						CH4.res[[i]]$fluss$flux <- ifelse(leak.flag[i], CH4.res[[i]]$fluss$flux, NA)
					}
				}	
			}	
		}
		flux.res$CH4 <- CH4.res
		## when pv values are to be reported… extract them
		if(asterisks){	
			CH4.pv <- sapply(CH4.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
			CH4.pv <- as.vector(symnum(CH4.pv, corr=FALSE, cutpoints = c(0,.001,.01,.05,.1,1), symbols = c("***","**","*","."," ")))
		}
		## make table for CH4
		CH4.table <- t(sapply(CH4.res, function(x) unlist(x$fluss[2:8])))
		CH4.units <- sapply(CH4.res, function(x) x$unit)
		CH4.table <- data.frame(CH4.units = unlist(CH4.units), CH4.pv, data.frame(CH4.table))
		CH4.table$rl <- CH4.rl
		names(CH4.table)[-c(1:2)] <- paste("CH4.", names(CH4.table)[-c(1:2)], sep="")
		flux.table <- data.frame(flux.table, CH4.table)	}
	cat(".",sep="")
	## N2O stuff
	if(length(no)!=0){
		## do prep for N2O
		var.par.N2O <- var.par
		## check for gcq and add 0 if necessary
		if(sum(names(var.par.N2O)=="N2O.gcq")==0) {
			var.par.N2O$N2O.gcq = 0
			warning("N2O GC quality flags have been set to zero")
		}
		# update no
		no <- grep("N2O", names(var.par.N2O))
		##
		names(var.par.N2O)[names(var.par.N2O)=="N2O"] <- "ghg"
		names(var.par.N2O)[names(var.par.N2O)=="N2O.gcq"] <- "gc.qual"
		## check about range.lim and attach calibration data to the data
		if(!is.null(range.lim)){
			N2O.rl <- vector("numeric", length(x))
			for(i in c(1:length(x))){
				N2O.rl[i] <- ifelse(length(range.lim$N2O)==1, range.lim$N2O, range.lim$N2O[i])
				x[[i]]$N2O.rl <- N2O.rl[i]
			}
		}
		else{
			N2O.rl <- sapply(x, function(y) y$N2O.rl[1])
		}
		N2O.pre <- lapply(x, function(x) flux.odae(x, var.par = c(var.par.N2O[no], vp), min.allowed = min.allowed, max.nrmse = max.nrmse$N2O, rl="N2O.rl"))
		## run the N2O flux estimation via flux.conv and gflux
		N2O.res <- lapply(N2O.pre, function(x) flux.conv(x, ghg = "N2O", r2.qual = r2.qual$N2O, nrmse.lim = nrmse.lim$N2O, out.unit = out.unit$N2O, elementar = elementar, hardflag = hardflag))
		if(co2ntrol$leak){
			for(i in c(1:length(N2O.res))){
				N2O.res[[i]]$fluss$leak.f <- leak.flag[i]
			}
		}
		if(exists("leak", hardflag)){
			if(hardflag$leak){
				if(length(co)==0){warning("leakage can only be evaluated and hardflagged if CO2 data are given")}
				else{
					for(i in c(1:length(N2O.res))){
						N2O.res[[i]]$fluss$flux <- ifelse(leak.flag[i], N2O.res[[i]]$fluss$flux, NA)
					}
				}	
			}	
		}		
		flux.res$N2O <- N2O.res
		## when pv values are to be reported… extract them
		if(asterisks){
		N2O.pv <- sapply(N2O.res, function(x) coef(summary(x$fl.dat$lm4flux))[2,4])
		N2O.pv <- as.vector(symnum(N2O.pv, corr=FALSE, cutpoints = c(0,.001,.01,.05,.1,1), symbols = c("***","**","*","."," ")))
		}
		## make table for N2O
		N2O.table <- t(sapply(N2O.res, function(x) unlist(x$fluss[2:8])))
		N2O.units <- sapply(N2O.res, function(x) x$unit)
		N2O.table <- data.frame(N2O.units, N2O.pv, N2O.table)
		N2O.table$rl <- N2O.rl
		names(N2O.table)[-c(1:2)] <- paste("N2O.", names(N2O.table)[-c(1:2)], sep="")
		flux.table <- data.frame(flux.table, N2O.table)
	}
	cat(".",sep="")
	## handthrough data
	sel <- c(!is.na(flux.res))
	htd <- data.frame(t(sapply(flux.res[[ which(sel==TRUE)[1] ]], function(x) x$fl.dat$dat.out)))
	htf <- data.frame(t(sapply(flux.res[[ which(sel==TRUE)[1] ]], function(x) x$fl.dat$fac.out)))
	nms <- c(names(flux.table), names(htd)[-c(1:2)], names(htf[-c(1:2)]))
	flux.table <- data.frame(flux.table, htd[,-c(1:2)], htf[,-c(1:2)])
	names(flux.table) <- nms
	## compile results for output
	res <- list(flux.res = flux.res, flux.table = flux.table, range.lim = range.lim)
	class(res) <- "fluxes"
	return(res)
}
