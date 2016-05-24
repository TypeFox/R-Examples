get.data.dims <-
function(genome.file)
{
	# Obtain the size of the genome file using UNIX function call.
        # Then parse the result, nonempty entry 1 is num of rows,
        #  and entry 2 is num of words (= nrows * ncols for a matrix).
        #  Thus num of cols = entry 2 / entry 1.
        nrows <- 0
        nrowscols <- 0

        wcounts <- try(system(paste("wc ", genome.file, sep=""), intern=TRUE))
        wlist <- unlist(strsplit(wcounts, " ", fixed=TRUE))

        if(length(wlist) < 3) {
                print(paste("error occurred, failed to compute number of rows and columns for the input file ", genome.file, sep=""))
                return(list(nrows=0, ncols=0))
        }
		
        # Since beginning can contain arbitrary number of spaces,
        # skip emptys. First non-empty is number of rows,
        # second non-empty is numrow*numcol
        for (v in (1:length(wlist)))
        {
                if (wlist[v] != "" && nrows == 0)
                        nrows <- as.numeric(wlist[v])

                else if (wlist[v] != "" && nrowscols == 0)
                        nrowscols <- as.numeric(wlist[v])
        }

        ncols <- nrowscols / nrows
	return(list(nrows=nrows, ncols=ncols))
}

