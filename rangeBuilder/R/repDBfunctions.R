# functions to interface with Reptile Database data
# http://reptile-database.reptarium.cz/


## Version of the synonym matcher that will search within a span from the present to a specified year for the presence of a synonym

synonymMatchByYear_repDB <- function(x, year1=1950, year2=1900, fuzzy = TRUE, fuzzyDist=2, advancedSearch = TRUE, returnMultiple = FALSE, printReport = TRUE, nthreads = 1) {
	
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}

	res <- vector(mode='character',length=length(x))
	
	x <- gsub(' ', '_', x)
	uniqueSp <- unique(x)

	#return NA for genus only, or species only
	res[grepl("_NA$|^NA_", x)] <- NA
	uniqueSp <- uniqueSp[!grepl("_NA$|^NA_", uniqueSp)]
	
	#if year = NULL, perform search on accepted names only
	if (is.null(year1)) {
		
		uniqueRes <- lapply(uniqueSp, function(x) firstPass(x, synList = RepDBlistYr, synonyms = FALSE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))

	} else {
		
		#accepted + synonym searching
		
		#pull synonyms that have been described since specified year
		synList <- lapply(RepDBlistYr, function(x) unique(x$syn[which(x$year >= year1)]))

		if (nthreads > 1) {
			cl <- parallel::makePSOCKcluster(nthreads)
			parallel::clusterExport(cl = cl, varlist = c('firstPass', 'uniqueSp', 'synList', 'fuzzy', 'fuzzyDist', 'returnMultiple'), envir = environment())
		uniqueRes <- parallel::parLapply(cl, uniqueSp, function(x) firstPass(x, synList = synList, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
			parallel::stopCluster(cl)
		} else {
			uniqueRes <- lapply(uniqueSp, function(x) firstPass(x, synList = synList, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
		}
		
	}
	
	## Advanced Searching
	# if some were not matched, send to advanced search (if enabled)
	advancedSp <- uniqueSp[which(sapply(uniqueRes, function(x) x[[2]]) == 'not found')]
	
	if (length(advancedSp) > 0 & advancedSearch) {
		
		# create version of synonym list with all combinations of genus/species and alternate latin endings
		synList <- lapply(RepDBlistYr, function(x) unique(x$syn[which(x$year >= year2)]))
		
		for (i in 1:length(synList)) {
			allgenera <- unlist(unique(lapply(synList[[i]], function(x) strsplit(x, split='_')[[1]][1])))
			allspecies <- unlist(unique(lapply(synList[[i]], function(x) strsplit(x, split='_')[[1]][2])))
	
			#add masculine/feminine variations
			allspecies <- append(allspecies, gsub('a$', 'um', grep('a$', allspecies, value=TRUE)))
			allspecies <- append(allspecies, gsub('a$', 'is', grep('a$', allspecies, value=TRUE)))
			allspecies <- append(allspecies, gsub('a$', 'us', grep('a$', allspecies, value=TRUE)))
			allspecies <- append(allspecies, gsub('um$|is$|us$', 'a', grep('um$|is$|us$', allspecies, value = TRUE)))
			allspecies <- unique(allspecies)
	
			#create all combinations
			allcombos <- cbind(rep(allgenera, each=length(allspecies)), rep(allspecies, times = length(allgenera)))
			synList[[i]] <- unique(paste(allcombos[,1], allcombos[,2],sep='_'))
		}
		
		if (nthreads > 1) {
			
			cl <- parallel::makePSOCKcluster(nthreads)
			parallel::clusterExport(cl = cl, varlist = c('firstPass', 'advancedSp', 'synList', 'fuzzy', 'fuzzyDist', 'returnMultiple'), envir = environment())
			advancedRes <- parallel::parLapply(cl, advancedSp, function(x) firstPass(x, synList = synList, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
			parallel::stopCluster(cl)
			
		} else {
			advancedRes <- lapply(advancedSp, function(x) firstPass(x, synList = synList, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
		}
		
		#adjust status
		for (i in 1:length(advancedRes)) {
			advancedRes[[i]][[2]] <- gsub('accepted name|synonym', 'advanced search', advancedRes[[i]][[2]])
		}

		#combine advanced results with first pass
		for (i in 1:length(advancedSp)) {
			uniqueRes[[which(uniqueSp == advancedSp[i])]] <- advancedRes[[i]]
		}
	}
	
	# fill in results vector
	for (i in 1:length(uniqueSp)) {
		res[which(x == uniqueSp[i])] <- uniqueRes[[i]][[1]]
	}
	
	if (printReport) {
		report <- sapply(uniqueRes, function(x) x[[2]])
		report <- table(report)
		cat('\n')
		for (i in 1:length(report)) {
			cat('\t', names(report)[i], '\t\t', report[i], '\n', sep = "")
		}
		cat('\n')
	}
	
	return(res)
}
	
	


######## country based queries #########


#Function to return countries by species
getRepDBcountryList <- function(spname) {
	spname <- gsub(" ", "_", spname)
	res <- lapply(spname, function(x) repDBcountryList_bySp[[x]])
	names(res) <- spname
	return(res)
}

#Function to return species by country
getRepDBSpFromCountry <- function(country) {
	res <- lapply(toupper(country), function(x) repDBcountryList[[x]])
	names(res) <- country
	return(res)
}

#Function to see list of countries
getRepDBcountries <- function() {
	return(names(repDBcountryList))
}


#----------------------------------


firstPass <- function(sp, synList, synonyms, fuzzy, fuzzyDist, returnMultiple) {
	
	found <- FALSE
	status <- 'not found'
	
	# check for exact match with accepted names
	if (sp %in% names(synList)) {
		res <- sp
		found <- TRUE
		status <- 'exact match | accepted name'
	}
	
	# check for fuzzy match with accepted names
	if (fuzzy & !found) {
		fmatch <- as.vector(adist(sp, names(synList), partial = FALSE, ignore.case = TRUE))
		if (length(which(fmatch <= fuzzyDist)) > 0) {
			found <- TRUE
			if (length(which(fmatch <= fuzzyDist)) == 1) {
				res <- names(synList)[which(fmatch <= fuzzyDist)]
				status <- 'fuzzy match | accepted name'
			} else {
				status <- 'multiple fuzzy hits | accepted name'
				if (returnMultiple) {
					res <- paste(names(synList)[which(fmatch <= fuzzyDist)], collapse = ' | ')
				} else {
					res <- NA
				}
			}
		}
	}			
		
	# check for strict match in synonyms from year1 to present
	if (!found & synonyms) {
		matches <- names(which(unlist(lapply(synList, function(x) sp %in% x)) == TRUE))
		if (length(matches) > 0) {
			found <- TRUE
			if (length(matches) == 1) {
				res <- matches
				status <- 'strict match | synonym'
			}
			if (length(matches) > 1) {
				if (returnMultiple) {
					res <- paste(matches, collapse = ' | ')
				} else {
					res <- NA
				}
				status <- 'multiple hits | synonym'
			}
		}
	}
	
	# check for fuzzy match in synonyms from year1 to present
	if (fuzzy & !found & synonyms) {
		matches <- lapply(synList, function(x) as.vector(adist(sp, x, partial = FALSE, ignore.case = TRUE)))
		matches <- names(which(unlist(lapply(matches, function(x) (any(x <= fuzzyDist))))) == TRUE)
		if (length(matches) > 0) {
			found <- TRUE
			if (length(matches) == 1) {
				res <- matches
				status <- 'fuzzy match | synonym'
			}
			if (length(matches) > 1) {
				res <- paste(matches, collapse = ' | ')
				status <- 'multiple fuzzy hits | synonym'
			}
		}
	}
	
	if (!found) {
		res <- NA
	}
	
	return(list(match = res, status = status))
}











