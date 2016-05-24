lod <-
function(x, markers=seq_len(x$nMark), theta=0, loop_breakers=NULL, max.only=FALSE, verbose=FALSE, tol=0.01) {
	stopifnot(inherits(x, 'linkdat'))
   if(inherits(x, 'singleton')) stop("This function is not applicable to singleton objects.")
	if (is.null(x$model)) 	stop("No model set.")
	if (x$nMark==0) 		stop("No marker data indicated.")
	if (x$hasLoops)		{	
		if(is.null(loop_breakers)) stop("The pedigree has loops. Please indicate loop breakers.")
		x = breakLoops(x, loop_breakers)
	}
	if (min(markers)<1 || max(markers)>x$nMark)	stop("Nonexistent marker indicated.")
   chrom = x$model$chrom
	
 	maxall = max(unlist(lapply(x$markerdata[markers], attr, 'nalleles')))
   if(maxall > 4) {
		if(is.character(theta) && theta=="max") stop("The option theta='max' is not implemented for markers with more than 4 alleles.")
		return(.lodm(x, markers=markers, theta=theta, loop_breakers=loop_breakers, max.only=max.only, verbose=verbose))
	}

	ilink <- function(x, marker) {
      m = x$markerdata[[marker]]
      nall = attr(m, 'nalleles')
      afreq = attr(m, 'afreq')
      initialCalc = .initialCalc(x, afreq=afreq, chrom)
      
      log_denom = stopl = likelihood_LINKAGE(x, m, logbase=10, initialCalc=initialCalc, TR.MATR=TR05[[nall]])
		if (log_denom==-Inf) return(c(NaN, NaN))
   
      startl = likelihood_LINKAGE(x, m, logbase=10, TR.MATR=TRZ[[nall]])
		optimal = optimize(likelihood_LINKAGE, c(0, 0.5), x=x, marker=m, afreq=NULL, initialCalc=initialCalc, logbase=10, tol=tol, maximum=TRUE)
		if (optimal$objective > max(startl,stopl)) {
			log_numer <- optimal$objective; theta_max <- optimal$maximum
		} else {
		log_numer = max(startl, stopl); theta_max <- c(0,.5)[which.max(c(startl, stopl))]
		}
		c(log_numer - log_denom, theta_max)
	}
	
	TR05 = lapply(1:maxall, function(n) .TRmatrNEW(0.5, n, chrom))
	
	map = .getMap(x, na.action=1, verbose=F)[markers, , drop=F]

	if (is.numeric(theta)) {
		stopifnot(max(theta)<=0.5, min(theta)>=0)
		if (verbose) 
			cat("Computing singlepoint LOD scores at each marker\nfor the following recombination ", ifelse(length(theta)==1, "value:\n", "values:\n"),
				paste("  theta =",theta,collapse="\n"), "\n", sep="")
		
		markerdata_list = x$markerdata[markers]
		trm_list = lapply(theta, function(th) lapply(1:maxall, function(n) .TRmatrNEW(th, n, chrom)))
		denoms = unlist(lapply(markerdata_list, function(m) 
                      likelihood_LINKAGE(x, m, theta=NULL, logbase=10, TR.MATR=TR05[[attr(m, 'nalleles')]]) ))
		numers = vapply(markerdata_list, function(m)  
                  unlist(lapply(trm_list, function(TR) 
                         likelihood_LINKAGE(x, marker=m, theta=NULL, logbase=10, TR.MATR=TR[[attr(m, 'nalleles')]]))), 
                  FUN.VALUE=numeric(length(theta)))
		res = numers - rep(denoms, each=length(theta))
		res = structure(res, dim=c(length(theta), length(markers)), dimnames = list(theta, map$MARKER), analysis="mlink", map=map, class="linkres")
	} 
	else if (identical(theta,"max")) {
		if (verbose) cat("Computing singlepoint LOD scores for each marker,\nmaximizing over all recombination values.\n\n") 
		TRZ = lapply(1:maxall, function(n) .TRmatrNEW(0, n, chrom))
      res = unlist(lapply(markers, ilink, x=x))
		res = structure(res, dim=c(2, length(markers)), dimnames = list(c("LOD", "t_max"), map$MARKER), analysis="ilink", map=map, class="linkres")
	}
	
	if (verbose) summary.linkres(res)
	if (max.only) ifelse(all(is.na(res)), NA, max(res, na.rm=TRUE))
	else res
}




.lodm <- function(x, markers=seq_len(x$nMark), theta=0, loop_breakers=NULL, max.only=FALSE, verbose=FALSE) {

	stopifnot(inherits(x,"linkdat"))
	if (is.null(x$model)) 	stop("No model set.")
	if (x$nMark==0) 		stop("No marker data indicated.")
	if (x$hasLoops)		{	
		if(is.null(loop_breakers)) stop("The pedigree has loops. Please indicate loop breakers.")
		x = breakLoops(x, loop_breakers)
	}

	theta_list = as.list(theta)
	lods = vapply(x$markerdata[markers], function(marker) {
		attr(marker, 'chrom') = x$model$chrom
		start.dat = .startdata_MD(x, marker)
		denom = likelihood.linkdat(x, locus1=marker, locus2="disease", theta=0.5, startdata = start.dat, logbase=10)
		unlist(lapply(theta_list, function(tht) likelihood.linkdat(x, locus1=marker, locus2="disease", theta=tht, startdata = start.dat, logbase=10) - denom))
	}, FUN.VALUE=numeric(length(theta)))
	
	map = .getMap(x, na.action=1, verbose=F)[markers, ,drop=FALSE]

	res = structure(lods, dim=c(length(theta), length(markers)), dimnames = list(theta, map$MARKER), analysis="mlink", map=map, class="linkres")
	if (verbose) summary(res)
	if (max.only) ifelse(all(is.na(res)), NA, max(res, na.rm=TRUE))
	else res
}

