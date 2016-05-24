## Function, that will take a taxon name and search for (in this order) strict matches, fuzzy matches, strict synonym matches, and fuzzy synonym matches. 

# Function for birds, mammals, amphibians


synonymMatch_other <- function(x, db, fuzzy = TRUE, fuzzyDist = 2, advancedSearch = TRUE, searchSynonyms = TRUE, returnMultiple = FALSE, printReport = TRUE, nthreads = 1) {

	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}

	if (!is.vector(x)) {
		stop('taxon must be a vector of taxon names.')
	}

	if (!db %in% c('birds', 'mammals', 'amphibians')) {
		stop('db must be either birds, mammals or amphibians.')
	}

	if (db == 'birds') {
		synonymDB <- birdList
	}
	if (db == 'mammals') {
		synonymDB <- mammalList
	}
	if (db == 'amphibians') {
		synonymDB <- amphibList
	}
	
	# convert NULL to NA
	for (i in 1:length(synonymDB)) {
		if (is.null(synonymDB[[i]])) {
			synonymDB[[i]] <- NA
		}
	}

	res <- vector(mode='character',length=length(x))
	
	x <- gsub(' ', '_', x)
	uniqueSp <- unique(x)

	#return NA for genus only, or species only
	res[grepl("_NA$|^NA_", x)] <- NA
	uniqueSp <- uniqueSp[!grepl("_NA$|^NA_", uniqueSp)]
	
	#if year = NULL, perform search on accepted names only
	if (!searchSynonyms) {
		
		uniqueRes <- lapply(uniqueSp, function(x) firstPass(x, synList = synonymDB, synonyms = FALSE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))

	} else {
		
		#accepted + synonym searching
		
		#reduce list to only those entries with synonyms
		#synonymBD <- synonymBD[sapply(synonymBD, function(x) !is.null(x))]

		if (nthreads > 1) {
			cl <- parallel::makePSOCKcluster(nthreads)
			parallel::clusterExport(cl = cl, varlist = c('firstPass', 'uniqueSp', 'synonymDB', 'fuzzy', 'fuzzyDist', 'returnMultiple'), envir = environment())
		uniqueRes <- parallel::parLapply(cl, uniqueSp, function(x) firstPass(x, synList = synonymDB, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
			parallel::stopCluster(cl)
		} else {
			uniqueRes <- lapply(uniqueSp, function(x) firstPass(x, synList = synonymDB, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
		}
		
	}
	
	## Advanced Searching
	# if some were not matched, send to advanced search (if enabled)
	advancedSp <- uniqueSp[which(sapply(uniqueRes, function(x) x[[2]]) == 'not found')]
	
	if (length(advancedSp) > 0 & advancedSearch) {
		
		# create version of synonym list with all combinations of genus/species and alternate latin endings
		for (i in 1:length(synonymDB)) {
			allgenera <- unlist(unique(lapply(synonymDB[[i]], function(x) strsplit(x, split='_')[[1]][1])))
			allspecies <- unlist(unique(lapply(synonymDB[[i]], function(x) strsplit(x, split='_')[[1]][2])))
	
			#add masculine/feminine variations
			allspecies <- append(allspecies, gsub('a$', 'um', grep('a$', allspecies, value=TRUE)))
			allspecies <- append(allspecies, gsub('a$', 'is', grep('a$', allspecies, value=TRUE)))
			allspecies <- append(allspecies, gsub('a$', 'us', grep('a$', allspecies, value=TRUE)))
			allspecies <- append(allspecies, gsub('um$|is$|us$', 'a', grep('um$|is$|us$', allspecies, value = TRUE)))
			allspecies <- unique(allspecies)
	
			#create all combinations
			allcombos <- cbind(rep(allgenera, each=length(allspecies)), rep(allspecies, times = length(allgenera)))
			synonymDB[[i]] <- unique(paste(allcombos[,1], allcombos[,2],sep='_'))
		}
		
		if (nthreads > 1) {
			
			cl <- parallel::makePSOCKcluster(nthreads)
			parallel::clusterExport(cl = cl, varlist = c('firstPass', 'advancedSp', 'synonymDB', 'fuzzy', 'fuzzyDist', 'returnMultiple'), envir = environment())
			advancedRes <- parallel::parLapply(cl, advancedSp, function(x) firstPass(x, synList = synonymDB, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
			parallel::stopCluster(cl)
			
		} else {
			advancedRes <- lapply(advancedSp, function(x) firstPass(x, synList = synonymDB, synonyms = TRUE, fuzzy = fuzzy, fuzzyDist = fuzzyDist, returnMultiple = returnMultiple))
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

