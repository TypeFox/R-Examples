# standardize country name to the set used in the rangeBuilder package

standardizeCountry <- function(country, nthreads = 1) {
	
	if (nthreads > 1) {
		if (!"package:parallel" %in% search()) {
			stop("Please load package 'parallel' for using the multi-thread option\n");
		}
	}

	country <- toupper(country)
	country <- gsub('_|-', ' ', country)
	country <- gsub('\\.', '', country)
	country <- stringi::stri_trans_general(country, "Latin-ASCII")
	country <- gsub('(^|\\s)ST\\.?\\s', 'SAINT ', country)
	country <- gsub('\\?|\\[|\\]', '', country)
	country <- gsub('\\s+', ' ', country)
	country <- trim(country)
	
	if (nthreads > 1) {
		cl <- parallel::makePSOCKcluster(nthreads)
		parallel::clusterExport(cl = cl, varlist = c('country', 'countryList'), envir = environment())
		res <- parallel::parSapply(cl, country, function(x) {
			if (x %in% unlist(countryList)) {
				ind <- which(sapply(countryList, function(y) x %in% y) == TRUE)
				return(names(countryList)[ind])
			} else {
				return('')
			}
		}, simplify = TRUE, USE.NAMES = FALSE)
		parallel::stopCluster(cl)
	} else {	
		res <- sapply(country, function(x) {
			if (x %in% unlist(countryList)) {
				ind <- which(sapply(countryList, function(y) x %in% y) == TRUE)
				return(names(countryList)[ind])
			} else {
				return('')
			}
		}, simplify = TRUE, USE.NAMES = FALSE)
	
	}
	
	return(res)
}
	


