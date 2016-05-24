linkage.power <- function(x, N=100, available=x$available, afreq=c(0.5, 0.5), 
                          loop_breakers=NULL, threshold=NULL, seed=NULL, verbose=FALSE) {
	if(inherits(x,"singleton")) stop("This function is not applicable to singleton objects.")
	if (length(available)==0) available = x$orig.ids
   sims = linkageSim(x, N=N, available=available, afreq=afreq, loop_breakers=loop_breakers, unique=FALSE, seed=seed, verbose=verbose)
	if(verbose) cat("\n")
	lds = lod(sims, theta=0, loop_breakers=loop_breakers, verbose=FALSE)
	res = structure( list(sim = sims, lod = lds, maxlod = max(lds), elod = mean(lds)), class="powres")
   summary(res, threshold=threshold)
   invisible(res)
}

summary.powres <- function(object, threshold=NULL, ...) {
	x <- object
   afrq = attr(x$sim$markerdata[[1]], 'afreq')
   cat(x$sim$nMark, "markers simulated with allele frequencies", .prettycat(afrq, 'and'), "\n")
   cat("Highest singlepoint LOD score:", x$maxlod, "\n")
   cat("Markers obtaining maximum score:", mean((x$maxlod - x$lod) < 1e-4)*100, "%\n")
   cat("ELOD:", x$elod, "\n")
   if (!is.null(threshold)) {
      cat("\nThreshold:", threshold, "\n")
      cat("Markers with score above threshold:", mean(x$lod >= threshold)*100, "%\n")
   }
}


power.varyPar <- function(x, N=100, varyPar, values, all=FALSE, loop_breakers=NULL, seed=NULL) {
	if(inherits(x,"singleton")) stop("This function is not applicable to singleton objects.")
	if (is.null(x$model)) stop("No model set.")
	if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
	
	if (all) x = setAvailable(x, x$orig.ids)
	
	models = switch(x$model$chrom,
		AUTOSOMAL = {switch(varyPar,
			f0 = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances[1]=k; mod}),
			f1 = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances[2]=k; mod}),
			f2 = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances[3]=k; mod}),
			dfreq = lapply(values, function(k) {mod=unclass(x$model); mod$dfreq=k; mod}),
			stop("Argument 'varyPar' must be one of 'f0', 'f1', 'f2', 'dfreq'.") ) },
		X = {switch(varyPar,
			f0_m = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$male[1]=k; mod}),
			f1_m = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$male[2]=k; mod}),
			f0_f = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$female[1]=k; mod}),
			f1_f = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$female[2]=k; mod}),
			f2_f = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$female[3]=k; mod}),
			dfreq = lapply(values, function(k) {mod=unclass(x$model); mod$dfreq=k; mod}),
			stop("Argument 'varyPar' must be one of 'f0_m', 'f1_m', 'f0_f', 'f1_f', 'f2_f', 'dfreq'.") ) } )

	sims = .SNPsim(x, N=N, loop_breakers=loop_breakers, seed=seed, unique=TRUE) 
	cat(N, "markers simulated;", ncol(sims$markerdata)/2, "unique.\n")

	lods = sapply(models, function(mod) max(lod(setModel(sims, model=mod), theta=0, loop_breakers=loop_breakers, verbose=FALSE)))
	plot(values, lods, type='l', xlab=varyPar, ylab='Maximum LOD')
	data.frame(value=values, LOD=lods)
}


.lod.varyParam <- function(x, markers=1, f0, f1, f2, dfreq, afreq1, theta, loop_breakers=NULL, plot=T, ...) {
	if(inherits(x, 'singleton')) stop("This function is not applicable to singleton objects.")
	if (is.null(x$model)) stop("No model set.")
	
	if(missing(f0)) f0 = x$model$penetrances[1]
	if(missing(f1)) f1 = x$model$penetrances[2]
	if(missing(f2)) f2 = x$model$penetrances[3]
	if(missing(dfreq)) dfreq = x$model$dfreq
	if(missing(afreq1)) {
		vary_a = FALSE
		afreq1 = -99  # dummy...anything of length 1
	} else {
		if(length(afreq1)>1 && length(markers)>1) stop("When varying marker allele frequencies, a single marker must be chosen")
		vary_a = TRUE
	}
	if(missing(theta)) theta=0
	
	arglist = list(f0,f1,f2,dfreq,afreq1,theta)
	arglength = lengths(arglist)
	varyP = which(arglength > 1)
	if(any(arglength < 1) || length(varyP) == 0 || length(varyP) > 2) stop("Something is wrong with the input parameters: See ?lod.varyParam")
	
	lods = numeric(0)
	for(f_0 in f0) for(f_1 in f1) for(f_2 in f2) for(d_frq in dfreq) for(a_frq in afreq1) {
		y = setModel(x, penetrances=c(f_0, f_1, f_2), dfreq=d_frq)
		if(vary_a) {
			nn = attr(y$markerdata[[markers]], 'nalleles')
			y = modifyMarker(y, marker=markers, afreq=c(a_frq, rep(1-a_frq, nn-1)/(nn-1)))
		}
		ld = lod(y, markers=markers, theta=theta, loop_breakers=loop_breakers, max.only=TRUE)  # theta can be of length > 1!
		lods = c(lods, ld)
	}
	shortnames = c('f0','f1','f2','d','a','th')
	longnames = c('Phenocopy rate', 'Penetrance parameter f1', 'Penetrance parameter f2', 'Frequency of disease allele', 'Frequency of marker allele 1', 'Recombination fraction (theta)')
	values1 = arglist[[varyP[1]]]
	names1 = paste(shortnames[varyP[1]], round(values1, 4), sep="_")
		
	switch(length(varyP), {
		result = structure(lods, names=names1)
		if(plot) plot(values1, result, xlab = longnames[varyP[1]], ylab="LOD", ...)
	},{
		values2 = arglist[[varyP[2]]]
		names2 = paste(shortnames[varyP[2]], round(values2, 4), sep="_")
		result = matrix(lods, byrow=T, ncol=arglength[varyP[2]], dimnames=list(names1, names2))
		if(plot)	contour(values1, values2, result, xlab = longnames[varyP[1]], ylab=longnames[varyP[2]], ...)
	})
	result	
}

