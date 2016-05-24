lengths.report <-
function(dat, numtodo = 10, reportInterval = 2000, high.mem = TRUE) {
## set numtodo to 0 if you want to do all loci
  if(class(dat) != 'pyRAD.loci') stop("This function runs on a pyRAD data object")
  last.lines <- dat$cons - 1
  num.loci <- length(last.lines)
  datSeqs <- as.character(dat$seqs)
  last.locus <- ifelse(numtodo < 1, num.loci, numtodo)
  if(high.mem) block.lengths <- sapply(datSeqs[last.lines][1:last.locus], function(x) nchar(as.character(x)))
  else {
    start.time <- Sys.time()
	block.lengths = integer(0)
	for(i in 1:last.locus) block.lengths = c(block.lengths, nchar(as.character(datSeqs[last.lines[i]])))
	if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', last.locus, 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (last.locus - i), attr(Sys.time() - start.time, 'units')
  	   ))
	   }

	}
  names(block.lengths) <- dat$locus.index[last.lines][1:last.locus]
  return(block.lengths)
  }
