clean.int <-
function(dat, v.avg.min, v.avg.max, dir.clean, turb.clean, icing, rep, n.rep) {
###	internal function for cleaning data
	
	if(!is.null(v.avg.min) && !is.null(dat$v.avg)) {
		replaced <- length(dat$v.avg[dat$v.avg<v.avg.min & !is.na(dat$v.avg)])
		if(dir.clean && !is.null(dat$dir.avg)) dat$dir.avg[dat$v.avg<v.avg.min] <- NA
		if(!is.null(dat$turb.int)) dat$turb.int[dat$v.avg<v.avg.min] <- NA
		dat$v.avg[dat$v.avg<v.avg.min] <- NA
		if(replaced>0) message(replaced, " samples lower than ", v.avg.min, " replaced by 'NA' in average wind speed")
	}
	if(!is.null(v.avg.max) && !is.null(dat$v.avg)) {
		replaced <- length(dat$v.avg[dat$v.avg>v.avg.max & !is.na(dat$v.avg)])
		#if(dir.clean && !is.null(dat$dir.avg)) dat$dir.avg[dat$v.avg>v.avg.max] <- NA
		dat$v.avg[dat$v.avg>v.avg.max] <- NA
		if(replaced>0) message(replaced, " samples higher than ", v.avg.max, " replaced by 'NA' in average wind speed")
	}
	if(dir.clean && !is.null(dat$dir.avg)) {
		replaced <- length(dat$dir.avg[(dat$dir.avg<0 | dat$dir.avg>360) & !is.na(dat$dir.avg)])
		dat$dir.avg[dat$dir.avg<0 | dat$dir.avg>360] <- NA
		if(replaced>0) message(replaced, " samples outside the range of 0-360 replaced by 'NA' in average wind direction")
	}
	if(!is.null(turb.clean) && !is.null(dat$turb.int)) {
		replaced <- length(dat$turb.int[dat$v.avg<turb.clean & !is.na(dat$v.avg) & !is.na(dat$turb.int)])
		dat$turb.int[dat$v.avg<turb.clean] <- NA
		if(replaced>0) message(replaced, " samples with average wind speed lower than ", turb.clean, " m/s replaced by 'NA' in turbulence intensity")
	}
	if(icing && !is.null(dat$dir.avg) && !is.null(dat$dir.std)) {
		replaced <- length(dat$dir.avg[dat$dir.std==0])
		dat$dir.avg[dat$dir.std==0] <- NA
		if(replaced>0) message(replaced, " samples with wind direction standard deviation = 0 replaced by 'NA' in average wind direction assuming icing")
	}
	if(!is.null(rep)) {
		replaced <- 0
		for(i in 1:length(rep)) {
			if(any(names(dat)==rep[i])) {
				clmn <- which(names(dat)==rep[i])
				l <- length(dat[,clmn])
				if(n.rep>l) warning("Cannot clean repetitions - 'n.rep' too large", call.=FALSE)
				else {
					for(n in 1:(l-n.rep+1)) {
						if(!is.na(dat[n,clmn]) && !is.na(dat[n+1,clmn])) if(dat[n,clmn]==dat[n+1,clmn]) {
							x <- 2
							while(n+x<=l) {
								if(is.na(dat[n+x,clmn])) break
								if(dat[n,clmn]!=dat[n+x,clmn]) break
								x <- x+1
							}
							if(x>=n.rep) {
								dat[(n+1):(n+x-1),clmn] <- NA
								replaced <- replaced + x-1
							}
						}
					}
				}		
			}
			if(replaced>0) message(replaced, " repetitions replaced by 'NA' in ", rep[i])
		}
	}
	
	if(is.null(attr(dat, "clean"))) {
		if(is.null(rep)) attr(dat, "clean") <- list(v.avg.min=v.avg.min, v.avg.max=v.avg.max, dir.clean=dir.clean, turb.clean=turb.clean, icing=icing)
		else attr(dat, "clean") <- list(v.avg.min=v.avg.min, v.avg.max=v.avg.max, dir.clean=dir.clean, turb.clean=turb.clean, icing=icing, rep=rep, n.rep=n.rep-1)
	} else {
		if(is.null(attr(dat, "clean")$v.avg.min)) attr(dat, "clean")$v.avg.min <- v.avg.min
		else if(!is.null(v.avg.min)) if(v.avg.min>attr(dat, "clean")$v.avg.min) attr(dat, "clean")$v.avg.min <- v.avg.min
		if(is.null(attr(dat, "clean")$v.avg.max)) attr(dat, "clean")$v.avg.max <- v.avg.max
		else if(!is.null(v.avg.max)) if(v.avg.max<attr(dat, "clean")$v.avg.max) attr(dat, "clean")$v.avg.max <- v.avg.max
		attr(dat, "clean")$dir.clean <- attr(dat, "clean")$dir.clean | dir.clean
		if(is.null(attr(dat, "clean")$turb.clean)) attr(dat, "clean")$turb.clean <- turb.clean
		else if(!is.null(turb.clean)) if(turb.clean>attr(dat, "clean")$turb.clean) attr(dat, "clean")$turb.clean <- turb.clean
		attr(dat, "clean")$icing <- attr(dat, "clean")$icing | icing
		if(is.null(attr(dat, "clean")$rep)) attr(dat, "clean")$rep <- rep
		else if(!is.null(rep)) attr(dat, "clean")$rep <- unique(attr(dat, "clean")$rep, rep)
		if(!is.null(rep) && is.null(attr(dat, "clean")$n.rep)) attr(dat, "clean")$n.rep <- n.rep-1
		else if(!is.null(rep) && !is.null(n.rep)) if(n.rep-1<attr(dat, "clean")$n.rep) attr(dat, "clean")$n.rep <- n.rep-1
	}
	
	return(dat)
}
