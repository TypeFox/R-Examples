combineToSingleMatrix<-function(outputPrefix, segFile, header)
{
	m = c()
	seg = read.table(segFile, sep = "\t", header = header)
	n = nrow(seg)
	error.message = "Not enough memory to allocate a single matrix from pairwise matrices!"
	tryCatch({	
		m = matrix(0, nrow = n, ncol = n)
    }, error = function(c) error.message , warning = function(c) error.message )
	if (length(m) == 0)
	{
		rm(m)
		return()
	}
    seg = seg[,1:3]
	colnames(seg) = c("chr", "start", "end")
	chrs = as.vector(unique(seg$chr))
	nchrs = length(chrs)
	for (i in 1:nchrs)
	{
		name1 = chrs[i]
        
        for (j in i:nchrs)
		{
			name2 = chrs[i]
            filename = paste(outputPrefix, "_", name1, "_", name2, ".txt", sep = "")
            mLocal = as.matrix(read.table(filename, sep = '\t'))
            indicesChr1 = which(seg$chr == name1)
            indicesChr2 = which(seg$chr == name2)
            i1 = indicesChr1[1]
            iN = indicesChr1[length(indicesChr1)]
            j1 = indicesChr2[1]
            jN = indicesChr2[length(indicesChr2)]
            m[i1:iN, j1:jN] = mLocal[1:nrow(mLocal), 1:ncol(mLocal)]
		}
	}
    ind <- lower.tri(m)
    m[ind] = t(m)[ind]
	write.table(m, paste(outputPrefix, "_entire_contat_matrix.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote =  F)
    cat("Completed writing single matrix\n")

}

buildCIM<-function(HiCFile, segFile, format, outputPrefix, resolution, header = FALSE, inclusive = FALSE, verbose = TRUE, combineToSingle = TRUE)
{
	 supportedFormats = c("sam", "nodup", "maq")
	 if (!file.exists(HiCFile) | !file.exists(segFile))
	 {
		cat("Error: cannot open HiC file or segmentation file.\n")
		return()
	 }
	 if (!(format %in% supportedFormats))
	 {
		cat(paste("Error: format:", format, "is not supported.\n"))
	 }
	 
	pathToScript = paste(find.package("chromoR") ,"/exec/parse_hic_to_matrix.py", sep = "")
    cmd = paste("python", pathToScript, segFile, HiCFile, format, outputPrefix, as.integer(header), as.integer(inclusive), as.integer(verbose), as.integer(resolution))
	 system(cmd)
	 if (verbose)
	 {
		cat("Completed writing contact pairwise matrices\n")
	 }
	 if (combineToSingle)
	 {
		cat("Trying to combine matrices to one single matrix...\n")
		combineToSingleMatrix(outputPrefix, segFile, header)
	 }

}





