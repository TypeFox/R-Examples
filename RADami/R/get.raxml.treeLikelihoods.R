get.raxml.treeLikelihoods <-
function(x, logfile = NA) {
## gets likelihoods from the RAxML_info file
## updated 2014-01-23 to deliver a fail if no trees written
	# cat(paste('working on file', x, '\n'), file = logfile)
	fileIn <- readLines(x)
	if(length(grep("Tree ", fileIn)) == 0) {
	  if(!is.na(logfile)) cat('... file', x, 'had no trees in it.', '\n')
	  return('FAIL')
	  }
	out <- as.double(sapply(grep("Tree ", fileIn, value=T), function(x) strsplit(x, ": ")[[1]][2]))
	names(out) <- as.character(1:length(out))
	out
	}
