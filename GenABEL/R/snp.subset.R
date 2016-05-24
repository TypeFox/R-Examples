snp.subset <- function(data,snpsubset) {
	if (missing(snpsubset)) {
		warning("snpsubset missing")
		return(data)
	}
	if (is(data,"scan.gwaa")) {
		out <- snp.subset.scan.gwaa(data,snpsubset)
	} else if (is(data,"check.marker")) {
		out <- snp.subset.check.marker(data,snpsubset)
	} else {stop("data should be of type check.marker-class or scan.gwaa-class")}
	out
}

"snp.subset.scan.gwaa" <- 
function(data,snpsubset) {
	norig <- length(snpnames(data))
	tokeep <- rep(FALSE,norig)
	if (is.numeric(snpsubset) || is.logical(snpsubset)) {
		tokeep[snpsubset] <- TRUE
	} else if (is.character(snpsubset)) {
		tokeep <- snpnames(data) %in% snpsubset
	} else {stop("snpsubset should have type numeric, logical or character")}
	out <- new("scan.gwaa",
			results=data@results[tokeep,],
			annotation = data@annotation[tokeep,], 
			lambda = data@lambda,
			idnames = data@idnames, 
			call = match.call(), 
			family = data@family
	) 
	out 
}

"snp.subset.check.marker" <- 
function(mcobj,subs) {
	o <- list()
# processing mcobj$call
	cidx <- match(subs,mcobj$call$name)
	o$call$name <- mcobj$call$name[cidx]
	o$call$map <- mcobj$call$map[cidx]
	o$call$chromosome <- mcobj$call$chromosome[cidx]
	o$call$call <- match.call()
# processing mcobj$details.redundancy
	for (i in names(mcobj$details.redundancy)) {
		if (any(subs==i))
		o$details.redundancy[i] <- mcobj$details.redundancy[i]
	}

# processing simple lists
	o$redundant <- subs[!is.na(match(subs,mcobj$redundant))]
	o$nofreq <- subs[!is.na(match(subs,mcobj$nofreq))]
	o$nocall <- subs[!is.na(match(subs,mcobj$nocall))]
	o$ok <- subs[!is.na(match(subs,mcobj$ok))]

# out hwe & chi2
	cidx <- !is.na(match(mcobj$nohwe,subs))
	o$nohwe <- mcobj$nohwe[cidx]
	o$chi2.hwe <- mcobj$chi2.hwe[cidx]

# output
	class(o) <- "check.marker"
	o
}
