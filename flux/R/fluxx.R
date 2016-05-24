fluxx <-
function(x, var.par, subset, asterisks = FALSE, loop = "auto", ...){
	## extract the name vector from x
	nmes <- data.frame(x$nmes)
	## extract the data tables from x
	x <- x$tables
	## make subsets
	mf <- match.call(expand.dots = FALSE)
	m <- match("subset", names(mf), 0L)
	if(m!=0){
		x <- x[subset]
		nmes <- nmes[subset,]
	}
	
	## which gas?
	ghgs <- c("CO2", "CH4", "N2O")
	ghg.sel <- match(ghgs, names(var.par))
	ghgs <- ghgs[!is.na(ghg.sel)]
	names(var.par)[names(var.par)==ghgs] <- "ghg"
	
	## check loop
	if(loop == "auto"){
		loop <- ifelse(length(x) >= 100, TRUE, FALSE)
	}

	## run the CO2 flux estimation via mf.flux
	if(loop){
		n <- length(x)
		ghg.res <- vector("list", n)
		for(i in c(1:n)){
			ghg.res[[i]] <- mf.flux(x[[i]], var.par = var.par, ...)
		}
	}
	else{
		ghg.res <- lapply(x, function(y) mf.flux(y, var.par = var.par, ...))
	}
	
	## when pv values are to be reported, extract them
	pv <- round(sapply(ghg.res, function(x) coef(summary(x$mod))[2,4]), 3) 
	if(asterisks){
		pv <- as.vector(symnum(pv, corr=FALSE, cutpoints = c(0,.001,.01,.05,.1,1), symbols = c("***","**","*","."," ")))
	}
	
	## make table for CO2
	isnum <- sapply(ghg.res[[1]]$fluss, is.numeric)
	islog <- sapply(ghg.res[[1]]$fluss, is.logical)
	ischar <- sapply(ghg.res[[1]]$fluss, is.character)
	ghg.tab.num <- t(sapply(ghg.res, function(x) unlist(x$fluss[isnum])))
	ghg.tab.chr <- t(sapply(ghg.res, function(x) unlist(x$fluss[ischar])))
	ghg.tab.log <- data.frame(t(sapply(ghg.res, function(x) unlist(x$fluss[islog]))))*1
	ghg.table <- data.frame(pv, ghg.tab.log, ghg.tab.chr, ghg.tab.num)
	names(ghg.table) <- paste(ghgs, names(ghg.table), sep=".")
	
	## prepare all in one big results table
	isnum <- sapply(ghg.res[[1]]$out, is.numeric)
	htd.chr <- t(sapply(ghg.res, function(x) unlist(x$out[!isnum])))
	htd.num <- t(sapply(ghg.res, function(x) unlist(x$out[isnum])))
	flux.table <- data.frame(nmes, ghg.table, htd.chr, htd.num)
	flux.table <- flux.table[,-grep("\\.1", names(flux.table))]
	
	## compile results for output
	flux.res <- vector("list", 1)
	flux.res[[1]] <- ghg.res
	names(flux.res) <- ghgs
	res <- list(flux.res = flux.res, flux.table = flux.table)
	class(res) <- "fluxxes"
	return(res)
}