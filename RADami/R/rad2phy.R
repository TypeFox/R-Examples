rad2phy <-
function(pyDat, inds = row.names(pyDat), loci = dimnames(pyDat)[[2]], outfile = 'pyMat.out.phy', padding = 50, verbose = FALSE) {
## makes a phylip-style data matrix from rad.mat output, limiting by individuals and loci
  if(class(pyDat) != "rad.mat") warning("I'm expecting output from rad.mat")
  outfile = file(outfile, "wt")
  open(outfile)
  cat(paste(length(inds), sum(sapply(pyDat[inds[1], loci], nchar)), "\n"), file = outfile) #header: number of individuals, number of bases
  for(i in inds) {
    if(!verbose) message(paste("Writing DNA line for individual", i))
	cat(i, file = outfile)
	cat(paste(rep(" ", padding - nchar(i)), collapse = ""), file = outfile)
	cat(paste(pyDat[i, loci], collapse = ""), file = outfile)
	cat("\n", file = outfile) # endline
	}
  close(outfile)
  }
