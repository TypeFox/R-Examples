

turf <- function(data, n, k, combos, ...) {
	
	##TURF Analysis for R
	##version 0.8-7 (built on R 3.1.2), depends: dplyr (>= 0.3.0), magrittr
	##updated 2014-12-01, Jack Horne (jack@jackhorne.net)
	
	#get data
	stime <- unclass(Sys.time())
	if (is.data.frame(data)) { datX <- data }
    else { datX <- read.table(data, header = TRUE, sep = "\t") }

	#error check initial arguments
	if(min(k) <= 1) stop("min(k) must be greater than or equal to 1")
	if(n <= max(k)) stop("max(k) must be less than the total number of items n")
	if((ncol(datX) - 2) < n) stop("Data must have respondent id, weight and at least n additional columns") 
	
	#collect and error check any additional arguments
	ctrlX <- turf.args()
	argsX <- list(...)
	
	if (length(argsX)) {
        ctrlargs <- names(formals(turf.args))
        indx <- match(names(argsX), ctrlargs, nomatch = 0L)
        if (any(indx == 0L)) {
            stop(gettextf("Argument %s not matched", names(argsX)[indx == 0L]), domain = NA) }
			
		ctrlX[names(argsX)] <- argsX
		
		##error check ctrlX
		if(length(ctrlX$depth) != 1 | ctrlX$depth < 1 | !is.numeric(ctrlX$depth)) {
			stop("depth must be a scalar greater than or equal to 1") }	
		if(length(ctrlX$keep) != 1 | ctrlX$keep < 0 | !is.numeric(ctrlX$keep)) {
				stop("keep must be a scalar greater than or equal to 0") }
		if(!is.logical(ctrlX$mc)) stop("mc must be TRUE/FALSE")
		if(!(ctrlX$sort %in% c("a", "d", "n"))) stop("sort must be %in% c('a', 'd', 'n')")
		
		ctrlX$depth <- floor(ctrlX$depth)
		ctrlX$keep <- floor(ctrlX$keep)
    }
	
	depth <- ctrlX$depth
	keep <- ctrlX$keep
	mc <- ctrlX$mc
	nsims <- ctrlX$nsims
	psims <- ctrlX$psims
	sort <- ctrlX$sort
	
	#generate TURF combos as required
	if(missing(combos)) {
		if(is.null(psims)) {
			psims <- colSums(datX[,-c(1:2)]) / sum(datX[,2])
			psims <- psims / sum(psims)
		}
	
		combos <- turf.combos(n, k, mc=mc, nsims=nsims, psims=psims)
	
	}
	
	#error check structure of combos
	if(!is.list(combos)) stop("combos must be a list structure")
	if(length(combos) != length(k)) stop("combos must have k components")
	
	
	#begin TURF analysis
	turf.agg <- list()
	
	for(i in 1:length(combos)) {
		gtime <- unclass(Sys.time())
		combos[[i]] <- as.matrix(combos[[i]])
		nR <- ifelse(sort == "n" | keep == 0, nrow(combos[[i]]), min(keep, nrow(combos[[i]])))
		
		if(ncol(combos[[i]]) != n) stop("combos[[i]] must have n columns")
		if(any(combos[[i]] != 0 & combos[[i]] != 1)) stop("combos[[i]] may contain only 0s and 1s")
		
		comb.expand <- kronecker(combos[[i]], rep(1, nrow(datX)))
		datX.expand <- kronecker(rep(1, nrow(combos[[i]])), as.matrix(datX))
		
		frqX <- rowSums(comb.expand * datX.expand[,-c(1:2)] * datX.expand[,2])
		
		rch.ind <- which(frqX >= depth * datX.expand[,2])
		rchX <- vector("numeric", length(frqX))
		rchX[rch.ind] <- datX.expand[rch.ind, 2]
		
		dagg <- as.data.frame(cbind(rchX, frqX))
		combo.ind <- kronecker(1:nrow(combos[[i]]), rep(1, nrow(datX)))
		
		dagg <- cbind(as.factor(combo.ind), dagg)
		names(dagg)[1] <- "combo.ind"
		
		turf.grp <- group_by(dagg, combo = combo.ind)
		turf.agg[[i]] <- summarise(turf.grp, rchX = sum(rchX), frqX = sum(frqX))
		turf.agg[[i]][,2:3] <- turf.agg[[i]][,2:3] / sum(datX[,2])
		turf.agg[[i]] <- cbind(turf.agg[[i]], combos[[i]])
		if(sort == "a" | sort == "d") {
			turf.agg[[i]] <- turf.agg[[i]][order(turf.agg[[i]][,2], turf.agg[[i]][,3], decreasing=(ctrlX$sort=="d")),] }
		turf.agg[[i]] <- turf.agg[[i]][1:nR,]
		rownames(turf.agg[[i]]) <- NULL
		
		cat(k[i], " of ", n, ": ", unclass(Sys.time()) - gtime, " sec", "\n", sep="")
		flush.console()
		
	}
	
	etime <- unclass(Sys.time())
	cat("total time elapsed:", etime - stime, "sec", "\n")
	
	return(list(turf=turf.agg, call=match.call()))
	
}


turf.args <- function(depth=1L, keep=0, mc=FALSE, nsims=10000, psims=NULL, sort="d") {

	#This function is present so arguments can be documented
	#Error checking occurs in turf() and turf.combos() functions
	
	return(list(depth=depth, keep=keep, mc=mc, nsims=nsims, psims=psims, sort=sort))

}


turf.combos <- function(n, k, ...) {

	if(min(k) <= 1) stop("min(k) must be greater than or equal to 1")
	if(n <= max(k)) stop("max(k) must be less than the total number of items n")
	
	#collect and error check any additional arguments
	ctrlX <- turf.args()
	argsX <- list(...)
	
	if (length(argsX)) {
        ctrlargs <- names(formals(turf.args))
        indx <- match(names(argsX), ctrlargs, nomatch = 0L)
        if (any(indx == 0L)) {
            stop(gettextf("Argument %s not matched", names(argsX)[indx == 0L]), domain = NA) }
			
		ctrlX[names(argsX)] <- argsX
		
		##error check ctrlX
		if(!is.logical(ctrlX$mc)) stop("mc must be TRUE/FALSE")
    }
	
	mc <- ctrlX$mc
	nsims <- ctrlX$nsims
	psims <- ctrlX$psims
	
	if(is.null(psims)) psims <- rep(1/n, n)
	if(length(psims) != n) stop("psims must be a vector of length n")
	if(any(psims < 0)) stop("psims must contain values greater than or equal to zero")

	combX <- list()
	
	for(j in 1:length(k)) {
	
		k[j] <- floor(k[j])
		comb.mat <- t(combn(n, k[j]))
	
		if(mc & nrow(comb.mat) > nsims) {
			#Monte Carlo simulated candidate set based on individual item probabilities
			if(nsims < 1) stop("nsims must be greater than zero for Monte Carlo sampling")
			if(sum(psims > 0) < k[j]) stop("psims must contain at least k[j] non-zero values")
			nsims <- floor(nsims)
			comb.mat <- matrix(0, nsims, k)
			for(m in 1:nsims) {
				xrow <- sample(1:n, k, prob=psims)
				comb.mat[m,] <- xrow
			}
		}
		
		dummy.mat <- matrix(0, nrow(comb.mat), n)

		for(i in 1:n) {
			dd <- which(comb.mat == i, arr.ind=TRUE)
			dummy.mat[dd[,1], i] <- 1
		}
		
		combX[[j]] <- dummy.mat
		
	}
	
	return(combX)

}







