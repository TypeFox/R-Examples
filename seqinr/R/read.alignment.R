#
# Read files of aligned sequences in various formats
#
read.alignment <- function(file, format, forceToLower = TRUE)
{
  #
  # Check that we have read permission on the file:
  #
  if(file.access(file, mode = 4) != 0) stop(paste("File", file, "is not readable"))
  
  ali <- switch( format,
	fasta = .Call("read_fasta_align", file, PACKAGE = "seqinr"), 
	FASTA = .Call("read_fasta_align", file, PACKAGE = "seqinr"), 
	mase = .Call("read_mase", file, PACKAGE = "seqinr"),
	MASE = .Call("read_mase", file, PACKAGE = "seqinr"),
	phylip = .Call("read_phylip_align", file, PACKAGE = "seqinr"),
	PHYLIP = .Call("read_phylip_align", file, PACKAGE = "seqinr"),
	msf = .Call("read_msf_align", file, PACKAGE = "seqinr"),
	MSF = .Call("read_msf_align", file, PACKAGE = "seqinr"),
	CLUSTAL = .Call("read_clustal_align", file, PACKAGE = "seqinr"),
	clustal = .Call("read_clustal_align", file, PACKAGE = "seqinr"),
	stop("Wrong format name: Format available are fasta,mase,phylip,msf,clustal")
  )

  ali <- lapply(ali, as.character)
  if(forceToLower) ali[[3]] <- lapply(ali[[3]], tolower)
  if(format == "mase"){
    ali <- list(nb = as.numeric(ali[[1]]), nam = ali[[2]], seq = ali[[3]], com = ali[[4]]) 
  } else {
    ali <- list(nb = as.numeric(ali[[1]]), nam = ali[[2]], seq = ali[[3]], com = NA)
  }
  class(ali) <- "alignment"
  return(ali)
}
