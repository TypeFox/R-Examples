# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date: December 2009
# Version 0.1
# Licence GPL v3

setClass('MaxEnt',
	contains = 'DistModel',
	representation (
		lambdas  = 'vector',
		results = 'matrix',
		path = 'character',
		html = 'character'
	),	
	prototype (	
		lambdas = as.vector(NA),
		results = as.matrix(NA),
		path = '',
		html = ''
	),
)



setClass('MaxEntReplicates',
	representation (
		models  = 'list',
		results = 'matrix',
		html = 'character'
	),	
	prototype (	
		models = list(),
		results = as.matrix(NA),
		html = ''
	),
)


setMethod ('show' , 'MaxEntReplicates', 
	function(object) {
		cat('class     :' , class(object), '\n')
		cat('replicates:', length(object@models), '\n')
		if (file.exists(object@html)) {
			browseURL( paste("file:///", object@html, sep='') )
		} else {
			cat('output html file no longer exists\n')
		}
	}
)	

		
		
		

setMethod ('show' , 'MaxEnt', 
	function(object) {
		cat('class    :' , class(object), '\n')
		cat('variables:', colnames(object@presence), '\n')
		# cat('lambdas\n')
		# print(object@lambdas)
#		pp <- nrow(object@presence)
#		cat('\npresence points:', pp, '\n')
#		if (pp < 5) { 
#			print(object@presence)
#		} else {
#			print(object@presence[1:5,])
#			cat('  (... ...  ...)\n')
#			cat('\n')
#		}
#		pp <- nrow(object@absence)
#		cat('\nabsence points:', pp, '\n')
#		if (pp < 5) {
#			print(object@absence)
#		} else {
#			print(object@absence[1:5,])
#			cat('  (... ...  ...)\n')
#			cat('\n')
#		}
#		cat('\nmodel fit\n')
#		print(object@results)
#		cat('\n')

		if (file.exists(object@html)) {
			browseURL( paste("file:///", object@html, sep='') )
		} else {
			cat('output html file no longer exists\n')
		}
	}
)	


if (!isGeneric("maxent")) {
	setGeneric("maxent", function(x, p, ...)
		standardGeneric("maxent"))
}	


.rJava <- function() {
	if (is.null(getOption('dismo_rJavaLoaded'))) {
		# to avoid trouble on macs
		Sys.setenv(NOAWT=TRUE)
		if ( requireNamespace('rJava') ) {
			rJava::.jpackage('dismo')
			options(dismo_rJavaLoaded=TRUE)
		} else {
			stop('rJava cannot be loaded')
		}
	}
}

.getMeVersion <- function() {
	jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
	if (!file.exists(jar)) {
		stop('file missing:\n', jar, '.\nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/')
	}
	.rJava()
	mxe <- rJava::.jnew("meversion") 
	v <- try(rJava::.jcall(mxe, "S", "meversion") )
	if (class(v) == 'try-error') {
		stop('"dismo" needs a more recent version of Maxent (3.3.3b or later) \nPlease download it here: http://www.cs.princeton.edu/~schapire/maxent/
		\n and put it in this folder:\n',
		system.file("java", package="dismo"))
	} else if (v == '3.3.3a') { 
		stop("please update your maxent program to version 3.3.3b or later. This version is no longer supported. \nYou can download it here: http://www.cs.princeton.edu/~schapire/maxent/'")
	}
	return(v)
}


setMethod('maxent', signature(x='missing', p='missing'), 
	function(x, p, silent=FALSE, ...) {
		v <- .getMeVersion()
		if (!silent) {
			cat('This is MaxEnt version', v, '\n' )
		}
		invisible(TRUE)
	}
)

setMethod('maxent', signature(x='SpatialGridDataFrame', p='ANY'), 
	function(x, p, a=NULL,...) {
		factors = NULL
		for (i in 1:ncol(x@data)) {
			if (is.factor(x@data[,i]) | is.character(x@data[,i])) { 
				factors = c(factors, colnames(x@data)[i]) 
			}
		}
		x <- brick(x)
		p <- .getMatrix(p)
		if (! is.null(a) ) { 
			a <- .getMatrix(a) 
		}
		# Signature = raster, ANY
		maxent(x, p, a, factors=factors, ...)
	}
)


.getMatrix <- function(x) {
	if (inherits(x, 'SpatialPoints')) {
		x <- data.frame(coordinates(x))
	} else if (inherits(x, 'matrix')) {
		x <- data.frame(x)
	}
	if (! class(x) == 'data.frame' ) {
		stop('data should be  a matrix, data.frame, or SpatialPoints* object')
	}
	if (dim(x)[2] != 2) {
		stop('presence or absence coordinates data should be a matrix or data.frame with 2 columns' ) 	
	}
	colnames(x) <- c('x', 'y')
	return(x)
} 


setMethod('maxent', signature(x='Raster', p='ANY'), 
	function(x, p, a=NULL, factors=NULL, removeDuplicates=TRUE, nbg=10000, ...) {

		p <- .getMatrix(p)
		if (removeDuplicates) {
			cells <- unique(cellFromXY(x, p))
			pv <- data.frame(extract(x, cells))
		} else {
			pv <- data.frame(extract(x, p))
		}

		lpv <- nrow(pv)
		pv <- stats::na.omit(pv)
		nas <- lpv - nrow(pv)
		if (nas > 0) {
			if (nas >= 0.5 * lpv) {
				stop('more than half of the presence points have NA predictor values')
			} else {
				warning(nas, ' (', round(100*nas/lpv,2), '%) of the presence points have NA predictor values')
			}
		} 
		
		if (! is.null(a) ) {
			a <- .getMatrix(a)
			av <- data.frame(extract(x, a))
			avr <- nrow(av)
			av <- stats::na.omit(av)
			nas <- length(as.vector(attr(av, "na.action")))
			if (nas > 0) {
				if (nas >= 0.5 * avr) {
					stop('more than half of the absence points have NA predictor values')
				} else {
					warning(nas, ' (', round(100*nas/avr, 2), '%) of the presence points have NA predictor values')
				}
			}
		} else { 
		# random absence
			if (is.null(nbg)) {
				nbg <- 10000 
			} else {
				if (nbg < 100) {
					stop('number of background points is very low')
				} else if (nbg < 1000) {
					warning('number of background points is very low')
				}
			}

			if (nlayers(x) > 1) {
				xy <- randomPoints( raster(x,1), nbg, p, warn=0 )
			} else {
				xy <- randomPoints(x, nbg, p, warn=0 )			
			}
			av <- data.frame(extract(x, xy))
			av <- stats::na.omit(av)
			if (nrow(av) == 0) {
				stop('could not get valid background point values; is there a layer with only NA values?')
			}
			if (nrow(av) < 100) {
				stop('only got:', nrow(av), 'random background point values; is there a layer with many NA values?')
			}
			if (nrow(av) < 1000) {
				warning('only got:', nrow(av), 'random background point values; Small exent? Or is there a layer with many NA values?')
			}
		}
		
		# Signature = data.frame, missing

		x <- rbind(pv, av)
		
		if (!is.null(factors)) {
			for (f in factors) {
				x[,f] <- factor(x[,f])
			}
		}
		
		p <- c(rep(1, nrow(pv)), rep(0, nrow(av)))
		maxent(x, p, ...)	
	}
)


.getreps <- function(args) {
	if (is.null(args)) { return(1) }
	args <- trim(args)
	i <- which(substr(args,1,10) == 'replicates')
	if (! isTRUE(i > 0)) {
		return(1)
	} else {
		i <- args[i]
		i <- strsplit(i, '=')[[1]][[2]]
		return(as.integer(i))
	}
}



setMethod('maxent', signature(x='data.frame', p='vector'), 
	function(x, p, args=NULL, path, silent=FALSE, ...) {
	
		MEversion <- .getMeVersion()

		x <- cbind(p, x)
		x <- stats::na.omit(x)
		x[is.na(x)] <- -9999  # maxent flag for NA, unless changed with args(nodata= ), so we should check for that rather than use this fixed value.

		p <- x[,1]
		x <- x[, -1 ,drop=FALSE]

		factors <- NULL
		for (i in 1:ncol(x)) {
			if (class(x[,i]) == 'factor') {
				factors <- c(factors, colnames(x)[i])
			}
		}
		
		if (!missing(path)) {
			path <- trim(path)
			dir.create(path, recursive=TRUE, showWarnings=FALSE)
			if (!file.exists(path)) {
				stop('cannot create output directory: ', path)
			}
			dirout <- path			
		} else {
			dirout <- .meTmpDir()
			f <- paste(round(runif(10)*10), collapse="")
			dirout <- paste(dirout, '/', f, sep='')
			dir.create(dirout, recursive=TRUE, showWarnings=FALSE)
			if (! file.exists(dirout)) {
				stop('cannot create output directory: ', f)
			}
		}
		
		pv <- x[p==1, ,drop=FALSE]
		av <- x[p==0, ,drop=FALSE]
		me <- new('MaxEnt')
		me@presence <- pv
		me@absence <- av
		me@hasabsence <- TRUE
		me@path <- dirout

		pv <- cbind(data.frame(species='species'), x=1:nrow(pv), y=1:nrow(pv), pv)
		av <- cbind(data.frame(species='background'), x=1:nrow(av), y=1:nrow(av), av)
		
		pfn <- paste(dirout, '/presence', sep="")
		afn <- paste(dirout, '/absence', sep="")
		write.table(pv, file=pfn, sep=',', row.names=FALSE)
		write.table(av, file=afn, sep=',', row.names=FALSE)

		mxe <- rJava::.jnew("mebridge")
		
		
		replicates <- .getreps(args) 
		args <- c("-z", args)

		if (is.null(factors)) {
			str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args)) 
		} else {
			str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args), rJava::.jarray(factors))
		}
		if (!is.null(str)) {
			stop("args not understood:\n", str)
		}

	
		if (replicates > 1) {
		
			mer <- new('MaxEntReplicates')
			d <- t(read.csv(paste(dirout, '/maxentResults.csv', sep='') ))
			d1 <- d[1,]
			d <- d[-1, ,drop=FALSE]
			dd <- matrix(as.numeric(d), ncol=ncol(d))
			rownames(dd) <- rownames(d)
			colnames(dd) <- d1
			mer@results <- dd
			f <- paste(dirout, "/species.html", sep='')
			html <- readLines(f)
			html[1] <- "<title>Maxent model</title>"
			html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
			html[3] <- sub("model for species", "model result", html[3])
			newtext <- paste("using 'dismo' version ", packageDescription('dismo')$Version, "& Maxent version")
			html[3] <- sub("using Maxent version", newtext, html[3])
			f <- paste(dirout, "/maxent.html", sep='')
			writeLines(html, f)	
			mer@html <- f
			
			for (i in 0:(replicates-1)) {	
				mex <- me
				mex@lambdas <- unlist( readLines( paste(dirout, '/species_', i, '.lambdas', sep='') ) )
					
				f <- paste(mex@path, "/species_", i, ".html", sep='')
				html <- readLines(f)
				html[1] <- "<title>Maxent model</title>"
				html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
				html[3] <- sub("model for species", "model result", html[3])
				newtext <- paste("using 'dismo' version ", packageDescription('dismo')$Version, "& Maxent version")
				html[3] <- sub("using Maxent version", newtext, html[3])
				f <- paste(mex@path, "/maxent_", i, ".html", sep='')
				writeLines(html, f)
				mex@html <- f
				mer@models[[i+1]] <- mex
				mer@models[[i+1]]@results <- dd[, 1+1, drop=FALSE]				
			}
			
			return(mer)
			
		} else {
			
			me@lambdas <- unlist( readLines( paste(dirout, '/species.lambdas', sep='') ) )
			d <- t(read.csv(paste(dirout, '/maxentResults.csv', sep='') ))
			d <- d[-1, ,drop=FALSE]
			dd <- matrix(as.numeric(d))
			rownames(dd) <- rownames(d)
			me@results <- dd
			
			f <- paste(me@path, "/species.html", sep='')
			html <- readLines(f)
			html[1] <- "<title>Maxent model</title>"
			html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
			html[3] <- sub("model for species", "model result", html[3])
			newtext <- paste("using 'dismo' version ", packageDescription('dismo')$Version, "& Maxent version")
			html[3] <- sub("using Maxent version", newtext, html[3])
			f <- paste(me@path, "/maxent.html", sep='')
			writeLines(html, f)	
			me@html <- f
		}
		
		me
	}
)


.meTmpDir <- function() {
	return( paste(raster::tmpDir(), 'maxent', sep="") )
}


.maxentRemoveTmpFiles <- function() {
	d <- .meTmpDir()
	if (file.exists(d)) {
		unlink(paste(d, "/*", sep=""), recursive = TRUE)
	}
}

setMethod("plot", signature(x='MaxEnt', y='missing'), 
	function(x, sort=TRUE, main='Variable contribution', xlab='Percentage', ...) {
		r <- x@results
		rnames <- rownames(r)
		i <- grep('.contribution', rnames)
		r <- r[i, ]
		names(r) <- gsub('.contribution', '', names(r))
		if (sort) {
			r <- sort(r)
		}
		dotchart(r, main=main, xlab=xlab, ...)
		invisible(r)
	}
)

