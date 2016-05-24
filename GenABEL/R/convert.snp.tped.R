"convert.snp.tped" <-
  function(tpedfile,tfamfile,outfile,strand="u",bcast=10000) {

	alcodes <- alleleID.codes()

	poss <- c("u","+","-","file")
	intstrand <- which(poss == strand)-1
	if (length(intstrand)!=1) {cat("strand argument must be one of",poss,"\n");stop();}

# strand == 0 -> unknown ('u')
# strand == 1 -> plus ('+')
# strand == 2 -> minus ('-')
# strand == 3 -> four columns (name,chr,pos,strand) expected at the beginning

    .C("convert_snp_tped",as.character(tpedfile),as.character(tfamfile),as.character(outfile),as.integer(intstrand),as.integer(bcast),as.character(alcodes),as.integer(length(alcodes)),PACKAGE="GenABEL")
    return(invisible(0))
  }
