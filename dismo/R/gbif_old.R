# Author: Robert J. Hijmans
# Date : December 2009-2011
# Version 1.0
# Licence GPL v3

# 2011-12-04
# implemented trycatch to deal with poor response from GBIF server
# suggestion and changed code provided by John Baumgartner

# 2013-01-15
# translate ISO2 codes to full country names
# add "cloc"
# 2013-06-19
# added 'species concept' option
# suggested by Aaron Dodd

.GBIFKey <- function(species) {
#	if (! require(XML)) { stop('You need to install the XML package to use this function') }
	url <- "http://data.gbif.org/ws/rest/taxon/list?dataproviderkey=1&rank=species&scientificname="
	url <- paste0(url, species)
   	doc <- XML::xmlInternalTreeParse(url)
	
	node <- XML::getNodeSet(doc, "//gbif:summary")
	m <- XML::xmlGetAttr(node[[1]], 'totalReturned')
	if (!(as.integer(m) > 0)) {
		return(NA)
	} else {
		node <- XML::getNodeSet(doc, "//tc:TaxonConcept")
		XML::xmlGetAttr(node[[1]], 'gbifKey')
	}
}


.gbif_old <- function(genus, species='', concept=FALSE, ext=NULL, args=NULL, geo=TRUE, sp=FALSE, removeZeros=FALSE, download=TRUE, getAlt=TRUE, returnConcept=FALSE,ntries=5, nrecs=1000, start=1, end=NULL, feedback=3) {
	
	
#	if (! require(XML)) { stop('You need to install the XML package to use this function') }

	gbifxmlToDataFrame <- function(s) {
		# this sub-funciton was hacked from xmlToDataFrame in the XML package by Duncan Temple Lang
		
		doc <- try(XML::xmlInternalTreeParse(s))

		nodes <- XML::getNodeSet(doc, "//to:TaxonOccurrence")
		if(length(nodes) == 0)   return(data.frame())
		varNames <- c("continent", "country", "stateProvince", "county", "locality",  "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters", "maximumElevationInMeters", "minimumElevationInMeters", "maximumDepthInMeters", "minimumDepthInMeters", "institutionCode", "collectionCode", "catalogNumber",  "basisOfRecordString", "collector", "earliestDateCollected", "latestDateCollected",  "gbifNotes")
		dims <- c(length(nodes), length(varNames)) 
   # create an empty data frame with as many rows and columns as needed.
		ans <- as.data.frame(replicate(dims[2], rep(as.character(NA), dims[1]), simplify = FALSE), stringsAsFactors = FALSE)
		names(ans) <- varNames
    # Fill in the rows based on the names.
		for(i in seq(length = dims[1])) {
			ans[i,] <- XML::xmlSApply(nodes[[i]], XML::xmlValue)[varNames]
		}

		nodes <- XML::getNodeSet(doc, "//to:Identification")
		varNames <- c("taxonName")
		dims <- c(length(nodes), length(varNames)) 
		tax <- as.data.frame(replicate(dims[2], rep(as.character(NA), dims[1]), simplify = FALSE), stringsAsFactors = FALSE)
		names(tax) <- varNames
    # Fill in the rows based on the names.
		for(i in seq(length = dims[1])) {
			tax[i,] <- XML::xmlSApply(nodes[[i]], XML::xmlValue)[varNames]
		}
		cbind(tax, ans)
	}

	tmpfile <- paste(tempfile(), '.XML')

	
	if (!is.null(ext)) { 
		ext <- round(extent(ext), 5)
		global <- extent(-180,180,-90,90)
		ex <- intersect(ext, global)
		if (!is.null(ex)) {
			ex <- paste('&minlatitude=', ex@ymin,'&maxlatitude=', 
			ex@ymax, '&minlongitude=', ex@xmin, '&maxlongitude=', ex@xmax, sep='')
		} else {
			warning('invalid extent')
		}
	} else {
		ex <- NULL
	}
	
	getkey <- TRUE
	if (is.logical(concept)) {
		genus <- trim(genus)
		species <- trim(species)
		gensp <- paste(genus, species)
		spec <- gsub("   ", " ", species) 
		spec <- gsub("  ", " ", spec) 	
		spec <- gsub(" ", "%20", spec)  # for genus species var. xxx
		spec <- paste(genus, '+', spec, sep='')
	} else {
		key <- round(as.numeric(concept))
		if (key < 1) stop('concept should be a positive integer')
		gensp <- concept
		concept <- TRUE
		getkey <- FALSE
	}
	
	if (sp) geo <- TRUE
	if (geo) { cds <- '&coordinatestatus=true' 
	} else { cds <- '' }
    base <- 'http://data.gbif.org/ws/rest/occurrence/'
	if (!is.null(args)) {
		args <- trim(as.character(args))
		args <- paste('&', paste(args, collapse='&'), sep='')
	}
	
	if (returnConcept) concept <- TRUE
	if (concept) {
		if (getkey) {
			key <- .GBIFKey(spec)
		}
		if (returnConcept) return(key)
		if (is.na(key)) {
			concept <- FALSE
			url <- paste(base, 'count?scientificname=', spec, cds, ex, args, sep='')
		} else {
			url <- paste(base, 'count?taxonconceptkey=', key, cds, ex, args, sep='')
		}
	} else {
		url <- paste(base, 'count?scientificname=', spec, cds, ex, args, sep='')
	}	
	
	tries <- 0
    while (TRUE)  {
	
		tries <- tries + 1
		if (tries > 5) { # if you cannot do this in 5 tries, you might as well stop
			stop('GBIF server does not return a valid answer after 5 tries')
		}
    	x <- try(readLines(url, warn = FALSE))
		if (class(x) != 'try-error') break
    }
    xn <- x[grep('totalMatched', x)]
	if (length(xn) == 0) {
		xe <- x[grep('gbif:exception', x)]
		if (length(xe)== 1) {
			xe <- unlist(strsplit(unlist(strsplit(xe, '>'))[2], '<'))[1]		
			cat(url, "\n")
			stop(xe)
		} else if (length(xe) > 1) {
			cat(url, "\n")
			stop(xe)
		} else {
			cat(url, "\n")
			stop("invalid request")
		}	
	}
    n <- as.integer(unlist(strsplit(xn, '\"'))[2])
    if (!download) {
        return(n)
    }
	
    if (n==0) {
		cat(gensp, ': no occurrences found\n')
        return(invisible(NULL))
    } else {
		if (feedback > 0) {
			cat(gensp, ':', n, 'occurrences found\n')
			flush.console()
		}
	}

	ntries <- min(max(ntries, 1), 100)
	if (! download) { return(n) }
	nrecs <- min(max(nrecs, 1), 1000)
	
	start <- max(1, start)
	if (start > n) {
		stop('"start" is larger than the number of records')
	}
	if (is.null(end)) {
		end <- n
	} else {
		stopifnot(end >= start)
	}
	
    iter <- n %/% nrecs
	breakout <- FALSE
	if (start > 1) {
		ss <- floor(start/nrecs)
	} else {
		ss <- 0
	}
	z <- NULL
	start <- start-1
    for (group in ss:iter) {
		if (group > 0) {
			start <- group * nrecs
			if (end < start) break
		}	
		if (group == iter) { 
			thisend <- min(end, n) - 1
			nrecs <- thisend-start+1
		} else { 
			thisend <- start+nrecs-1 
			thisend <- min(end, thisend)
		}
		
		if (feedback > 1) {
			if (group == ss) { 
				cat((start+1), '-', thisend+1, sep='')  
			} else { 
				cat('-', thisend+1, sep='')  
			}
			if ((group > ss & group %% 20 == 0)  |  group == iter ) { cat('\n') }
			flush.console()
		}

		if (concept) {
			aurl <- paste(base, 'list?taxonconceptkey=', key, '&mode=processed&format=darwin&startindex=', format(start, scientific=FALSE), '&maxresults=', format(nrecs, scientific=FALSE), cds, ex, args, sep='')
		} else {
			aurl <- paste(base, 'list?scientificname=', spec, '&mode=processed&format=darwin&startindex=', format(start, scientific=FALSE), '&maxresults=', format(nrecs, scientific=FALSE), cds, ex, args, sep='')
		}
		tries <- 0
        #======= if download fails due to server problems, keep trying  =======#
        while (TRUE) {
			tries <- tries + 1
			if (tries > ntries) {
				warning('GBIF did not return the data in ', ntries, '  tries for:')
				print(aurl)
				breakout <- TRUE
				break
			}
			test <- try (download.file(aurl, tmpfile, quiet=TRUE))
			if (class(test) == 'try-error') {
				print('download failure, trying again...')
			} else {
				xml <- scan(tmpfile, what='character', quiet=TRUE, sep='\n')
				xml <- chartr('\a\v', '  ', xml)
				zz <- try( gbifxmlToDataFrame(xml))
				if (class(zz) == 'try-error') {
					print('parsing failure, trying again...')
				}
				break
			}
	    }
		
		if (breakout) {
			break
		} else {
			z <- rbind(z, zz)
		}
	}

	d <- as.Date(Sys.time())
	z <- cbind(z, d)
	names(z) <- c("species", "continent", "country", "adm1", "adm2", "locality", "lat", "lon", "coordUncertaintyM", "maxElevationM", "minElevationM", "maxDepthM", "minDepthM", "institution", "collection", "catalogNumber",  "basisOfRecord", "collector", "earliestDateCollected", "latestDateCollected",  "gbifNotes", "downloadDate")
	z[,'lon'] <- gsub(',', '.', z[,'lon'])
	z[,'lat'] <- gsub(',', '.', z[,'lat'])
	z[,'lon'] <- as.numeric(z[,'lon'])
	z[,'lat'] <- as.numeric(z[,'lat'])

	k <- apply(z[ ,c('lon', 'lat')], 1, function(x) isTRUE(any(x==0)))

	if (removeZeros) {
		if (geo) {
			z <- z[!k, ]
		} else {
			z[k, c('lat', 'lon')] <- NA 
		}
	} else {
		z[k, c('lat', 'lon')] <- NA 
	}
		
	if (getAlt) {
		altfun <- function(x) {
					a <- mean(as.numeric(unlist(strsplit( gsub('-', ' ', gsub('m', '', ( gsub(",", "", gsub('\"', "", x))))),' ')), silent=TRUE), na.rm=TRUE)
					a[a==0] <- NA
					mean(a, na.rm=TRUE)
				}

		#elev <- apply(z[,c("maxElevationM", "minElevationM")], 1, FUN=altfun)
		#depth <- -1 * apply(z[,c("maxDepthM", "minDepthM")], 1, FUN=altfun)
		#alt <- apply(cbind(elev, depth), 1, FUN=function(x)mean(x, na.rm=TRUE))
		
		if (feedback<3) {
			w <- options('warn')
			options(warn=-1)
		}
		
		alt <- apply(z[,c("maxElevationM", "minElevationM", "maxDepthM", "minDepthM")], 1, FUN=altfun)
		
		if (feedback<3) options(warn=w)
	
		z <- cbind(z[,c("species", "continent", "country", "adm1", "adm2", "locality", "lat", "lon", "coordUncertaintyM")], 
		alt, 
		z[ ,c("institution", "collection", "catalogNumber",  "basisOfRecord", "collector", "earliestDateCollected", "latestDateCollected",  "gbifNotes", "downloadDate", "maxElevationM", "minElevationM", "maxDepthM", "minDepthM")])
	}

	if (dim(z)[1] > 0) {
	
		iso <- ccodes()
		z$ISO2 <- z$country
		i <- match(z$ISO2, iso[, 'ISO2'])
		z$country <- iso[i, 1]
		
		fullloc <- trim(as.matrix(z[, c('locality', 'adm1', 'adm2', 'country', 'continent')]))
		fullloc <- apply(fullloc, 1, function(x) paste(x, collapse=', '))
		fullloc <- gsub("NA, ", "", fullloc)
		fullloc <- gsub(", NA", "", fullloc)
		fullloc <- gsub('\"', "", fullloc)
		z$cloc <- fullloc

		if (sp) {
			coordinates(z) <- ~lon+lat
		}
	}	
#	if (inherits(ext, 'SpatialPolygons')) { overlay	}
	try(file.remove(tmpfile), silent=TRUE)
	
	return(z)
}

#sa <- gbif('solanum')
#sa <- gbif('solanum', '*')
#sa <- gbif('solanum', 'acaule*')
#sa <- gbif('solanum', 'acaule var acaule')

