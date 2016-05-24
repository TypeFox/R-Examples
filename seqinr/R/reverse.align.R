reverse.align <- function(nucl.file,
                          protaln.file,
                          input.format = "fasta",
                          out.file,
                          output.format = "fasta",
                          align.prot = FALSE,
                          numcode = 1,
                          clustal.path = NULL,
                          forceDNAtolower = TRUE,
                          forceAAtolower = FALSE){
  #
  # Sequence import section
  #
  seq.nucl <- read.fasta(nucl.file, forceDNAtolower = forceDNAtolower)
  nseqs <- length(seq.nucl) # number of sequences
  
  if(!isTRUE(align.prot)){  ## the protein alignment file is provided 
    protaln <- read.alignment(protaln.file, format = input.format,
                              forceToLower = forceAAtolower)
  } else { ## protein alignment file not provided, have to align with clustal 
    tmp <- tempfile(pattern = "clustal")
    protseq.file <- tempfile(pattern = "protein")
    write.fasta(sequences = lapply(seq.nucl, function(x)
      translate(x, numcode = numcode)), names = names(seq.nucl), file.out = protseq.file) 
    system(paste(clustal.path, " -outfile=", tmp ," -infile=", protseq.file, sep = ""))
    protaln <- read.alignment(tmp, format = "clustal", forceToLower = forceAAtolower)
    input.format <- "clustal"
  }
  #
  # Force sequences to be in the same order in the nucleic and protein
  # versions. Use the protein order.
  #
  ordername <- unlist(lapply(protaln$nam, function(x) which(names(seq.nucl) == x)))
  seq.nucl <- seq.nucl[ordername]
  #
  # The character used to represent gaps is function of the alignment file format.
  #
  gapchar <- NULL
  if(input.format %in% c("fasta", "clustal", "phylip", "mase")){
    gapchar <- "-"
  }
  if(input.format == "msf"){
    gapchar <- "."
  }
  if(is.null(gapchar)) stop("no known gap character for this alignment format")
  #
  # Memory allocation, cds.aln is the result to write in file.out
  #
  cds.aln <- vector(mode = "list", length = nseqs)
  names(cds.aln) <- protaln$nam
  #
  # index[[j]] is for current position in nucleic acid sequences number j
  # expressed in codon units.
  #
  index <- as.list(rep(0, nseqs))
  names(index) <- protaln$nam
  #
  # Main loop to build the reverse alignment.
  # "allaln" is a flag still TRUE if no gap was found in a column.
  #
  ncharprot <- nchar(protaln$seq[1]) # number of AA + gaps
  for(k in seq_len(ncharprot)){
    allaln <- TRUE
    for(j in seq_len(nseqs)){
      if(substr(protaln$seq[j], k, k) != gapchar){
        index[[j]] <- index[[j]] + 1
      } else {
        allaln <- FALSE
      }
    }
    if(allaln){
      #
      # There was no gap in this column, we include the corresponding codon
      # in the nucleic acid alignment:
      #
      for(j in seq_len(nseqs)){
        cds.aln[[j]] <- c(cds.aln[[j]], seq.nucl[[j]][(3*index[[j]]-2):(3*index[[j]])])
      }
    }
  }
  #
  # Write into output file the reverse alignment.
  #
  write.fasta(sequences = cds.aln, names = names(seq.nucl), file.out = out.file, open = "w")
  return(NULL)
}
