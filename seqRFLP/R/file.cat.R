file.cat <-
function (dir = NULL, appendix = ".fas", file = NULL) {
	if(is.null(file))
	{stop("You have to specify a file name.")}
	fname = list.files(dir)
	fname = fname[regexpr(appendix, fname) > 0]
    for (k in 1:length(fname)) {
        dna = readLines(file.path(dir, fname[k]))
        if (k == 1) 
            DNA = dna
        if (k > 1) 
            DNA = c(DNA, dna)
    }
	# class(DNA) <- "fasta"
    writeLines(DNA, con = file)
}

