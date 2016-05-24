flux.odae <-
function(dat, var.par, min.allowed = 3, max.nrmse = 0.1, rl = NULL){
	## extract range limits (if they are there)
	if(is.null(rl) & !exists("rl", dat)){
		rl <- 0
		warning("No range limit defined. Has been set to 0.") 
	}
	if(!is.numeric(rl)){ rl <- dat[,rl]}
	if(exists("rl", dat)){rl <- dat$rl}
	## function vp that realises the compilation of dat for further use
	vp <- function(dat, sel){
		if(is.character(sel)){out <- dat[,sel]}
		else{out <- rep(sel, nrow(dat))}
		return(out)
	}
	## using vp to extract variables from dat and combine with fixed parameters
	dat <- lapply(var.par, function(x) vp(dat, x))
	## organize handthrough
	stv <- match(c("ghg", "time", "gc.qual", "area", "volume"), names(dat))
	handthrough <- names(which(sapply(var.par[-stv], is.character)))
	hdt.fac.sel <- sapply(dat[handthrough], is.numeric)
	fac.out <- c(CH4 = "methane",CO22 = "carbon dioxide")
	fac.out <- c(fac.out, sapply(dat[handthrough][!hdt.fac.sel], function(x) as.character(x[1])))
	hdt <- c(dat[c("area", "volume")], dat[handthrough][hdt.fac.sel])
	dat.out <- round(sapply(hdt, mean), 3)
	## extract variables from dat for later flux calculation
	dat <- sapply(dat[c("ghg", "time", "gc.qual", "area", "volume", "t.air", "p.air")], function(x) x[])
	dat <- data.frame(dat)
	## sort dat according to time column (just in case the dara came weird)
	dat <- dat[order(dat$time),]
	## getting all possible combinations of at least three x
	vers <- unlist(lapply(c(min.allowed:nrow(dat)), function(x) combn(c(1:nrow(dat)), x, simplify=FALSE)), recursive=FALSE)
	vers <- lapply(vers, "sort")
	## when 2 or more measurements were taken at a time the number
	## of min.allowed has to be extended to the number of unique
	## time steps. This is tested for and corrected in the following
	sel <- unlist(lapply(vers, function(x) length(unique(dat[x,2]))))
	vers <- vers[sel >= min.allowed]
	## determine the number of concentration measurements per version
	vers.n <- sapply(vers, "length")
	## now we can calculate the regression models
	vers.lm <- lapply(vers, function(x) lm(ghg ~ time, data=dat[x,]))
	## and extract the statistic nrmse
	nrmse <- unlist(lapply(vers.lm, function(x) sqrt(sum(residuals(x)^2)/summary(x)$df[2])/diff(range(x$model[1], na.rm=TRUE))))
	## rank according to nrmse
	ranks <- order(vers.n, 1-nrmse, decreasing=TRUE)
	m2t <- ranks[nrmse[ranks] <= max.nrmse][1]
	## avoid problems with NA's
	if(is.na(m2t)){	
		m2t <- order(nrmse)[1]
	}
	## select and store the best model
	lm4flux <- vers.lm[[m2t]]
	## store measurements that are used in the best model
	row.select <- vers[[m2t]]
	names(dat)[1] <- var.par$ghg
	dat$rl <- rl
	res <- list(lm4flux = lm4flux, row.select = row.select, orig.dat = dat, dat.out = dat.out, fac.out = fac.out)
	return(res)
	}

