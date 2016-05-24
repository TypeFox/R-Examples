raw2Fas <-
function (dir = NULL, appendix = ".fasta") 
{
    if(is.null(dir))
	{stop("You have to choose a directory where sequence files exist.")}
	fname = list.files(dir)
	fname = fname[grepl(appendix, fname)]
    for (k in 1:length(fname)){
        rawdat = readLines(file.path(dir, fname[k]))
		dna = rawdat[!grepl(">",rawdat)]
		stri = paste(dna, collapse = "")
		if(!any(grepl(">", rawdat))){
		   warning(paste("The input file ",fname[k],"seems not in FASTA format.\n",
		   "its file name has been added as the sequence's name. \n"))
		   fnam = gsub(appendix, "", fname[k])
		   fnam <- paste(">seq_", fnam, sep = "")
		}
		else{
		fnam = rawdat[grepl(">",rawdat)]
		}
		if(grepl("[0-9]", stri)){
		  warning(paste("Numbers appeared the dna sequence in file", fname[k], ".\n Is this file contains only dna sequence data? "))
		}
        if (k == 1){
            DNA = c(fnam, stri)
		}
        if (k > 1){ 
            DNA = c(DNA, fnam, stri)
        }
	}
	class(DNA) <- "fasta" 
    return(DNA)
}

