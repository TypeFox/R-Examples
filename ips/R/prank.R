prank <-
function(x, outfile = "PRANK", guidetree = NULL, gaprate = 0.025, gapext = 0.75,
         path){
	
	rwd <- getwd()
	
	# check and handle guidetree
	# ----------------------------
	if (class(guidetree) == "phylo"){
		
		missingseqs <- which(!guidetree$tip.label %in% x$nam)
		if (length(missingseqs) > 0)
			guidetree <- drop.tip(guidetree, missingseqs)
		missingtips <- which(!x$nam %in% guidetree$tip.label)
		if (length(missingtips) > 0){
			txt <- paste(x$nam[missingtips], 
                   "not contained in guide tree")
			stop(txt)
		}
		if (!all(guidetree$tip.label %in% x$nam))
				stop("Guidetree does not match sequences")
	
		write.tree(guidetree, "prank_guidetree.tre")
		gtfile <- paste(rwd, "/prank_guidetree.tre", sep = "")
		
		
		# check and handle sequences
		# ----------------------------
		id <- match(guidetree$tip.label, x$nam)
		x$nam <- x$nam[id]
		x$seq <- x$seq[id]
	}
	
	# write FASTA file
	# ----------------------------
	write.fas(x, "prankinput.fas", interleave = FALSE)
	infile <- paste("'", rwd, "/prankinput.fas'", sep = "")
	
	# call PRANK
	# ----------------------------
	if (class(guidetree) == "phylo") {
	  call.prank <- paste(path, " -F -d=", infile, " -t=", gtfile, 
                        " -o=", outfile, " -gaprate=",   		
	                      gaprate, " -gapext=", gapext, sep = "")
	}	
  else {
    call.prank <- paste(path, " -F -d=", infile, " -o=", outfile, 
                        " -gaprate=", gaprate, " -gapext=", gapext, sep = "")	
  }
	system(call.prank)
	
	# read alignment back in
	# ----------------------------
	if (class(guidetree) == "phylo")
		fn <- paste(outfile, "1.fas", sep = ".")	
  else 												
    fn <- paste(outfile, "2.fas", sep = ".")
	read.fas(fn)
}