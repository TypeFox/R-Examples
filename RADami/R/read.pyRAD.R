read.pyRAD <-
function(filename, reportInterval = 20000, breakLinesSeparate = FALSE, doSummary = TRUE, ...) {
## reads the all.aligned file out of pyRAD, parses into names, loci, sequences
## updated with breakLinesSeparate in Oct 2012 because pyRAD switched to single-line summaries at the end of each aligned file
## updated 2012-11-16 to keep breaklines intact
## updated 2012-12-19 to get rid of all conversions to factors... apparently no longer needed for space considerations in R
  message("Reading data...")
  dat <- readLines(filename, ...)
  dat.breakLines <- dat.consensusLines <- grep("//", dat, fixed = TRUE) # this is slow, but only ca. 1 sec for data runs of 10s of thousands
  dat.breakLines.vector <- dat[dat.breakLines] #added 2012-11-16; ignores possibility of separate breakLines
  message("Splitting data...")
  dat.split <- strsplit(dat, " {1,100}") #uses whitespace to separate taxon names from sequences
  dat.names <- sapply(dat.split, function(x) x[1])
  dat.seqs <- sapply(dat.split, function(x) x[2])
  # dat.seqs[dat.breakLines] <- dat.breakLines.vector # shoves the consensus seqs back into the sequence vector, assuming breakLinesSeparate = F
  dat.firstLocusLines <- c(1, (dat.breakLines[1:(length(dat.breakLines)-1)] + 1))
  if(breakLinesSeparate) {
    dat.lastLocusLines <- dat.breakLines - 2
    dat.consensusLines <- dat.breakLines - 1
	}
  else dat.lastLocusLines <- dat.breakLines - 1
  dat.locus.index <- character(length(dat))
  locusCounter <- 1
  message("Assigning locus number...")
  start.time <- Sys.time()
  ## crazy slow! before vectorization:
  #for (i in 1:length(dat)) {
  #  if(i-1 %in% dat.breakLines) locusCounter <- locusCounter + 1
  #	dat.locus.index[i] <- paste("locus.", locusCounter, sep = "")
  #	if(i / reportInterval - i %/% reportInterval == 0) {
  #	   message(paste('...', i, 'of', length(dat.locus.index), 
  #	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (length(dat.locus.index) - i), attr(Sys.time() - start.time, 'units')
  #	   ))
  #	   }
  #	}

  ## after some vectorization:
  for (i in 1:length(dat.firstLocusLines)) {
    dat.locus.index[dat.firstLocusLines[i]:dat.lastLocusLines[i]] <- names(dat.breakLines.vector)[i] <- paste("locus.", i, sep = "")
	if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', length(dat.firstLocusLines), 
 	   '-- Estimated time remaining =', round(((Sys.time() - start.time) / i) * (length(dat.firstLocusLines) - i), 1), attr(Sys.time() - start.time, 'units')
  	   ))
	}
  }
  # dat.locus.index <- as.factor(dat.locus.index) # done only for memory considerations... slows things down in analysis unless you transform back to character
  out = list(tips = dat.names, seqs = dat.seqs, breaks = dat.breakLines, break.vectors = dat.breakLines.vector, cons = dat.consensusLines, locus.index = dat.locus.index, file.read = filename, timestamp = date())
  class(out) <- 'pyRAD.loci'
  if(doSummary) out$radSummary <- summary(out)
  return(out)
  }
