CNOGpro <-
function(hitsfile, gbkfile=stop("You must provide a GenBank reference file"), windowlength=100, name="Default organism"){
  copynumber <- list()
  class(copynumber) <- "CNOGpro"
  copynumber$Name <- name
  copynumber$windowlength <- windowlength
  
  # Read Genbankfile and create table of genes and intergenic regions:
  cat ("Reading contents of GenBank file...\n")
  table <- readGenBank(gbkfile)
  table <- updateGeneTable(table)
  
  cat ("Finished creating table of genes!\n")
  
  # Retrieve the FASTA sequence:
  cat ("Storing fasta sequence...\n")

  gbkfilename <- basename(gbkfile)  
  fullsequence <- getsequence(paste(tempdir(),"/",gsub(".gbk|.gb|.gff|.genbank|.txt|.text","",gbkfilename),".fasta",sep=""))
  
  # Load hits
  cat("Attempting to read best-hit read-location files. This might take a while...\n")
  if(missing(hitsfile)) {cat("No hits-file provided. You will not be able to perform copy number analysis, but can still inspect genetic elements.\n")}
  else{hits <- gethits(hitsfile)}
  
  #Everything that is per-chromosome...
  copynumber$accession <- attr(fullsequence, which="name")
  copynumber$genes <- table
  copynumber$chrlength <- table$Right[nrow(table)]

  # Calculate GC:
  cat ("Calculating GC percentage in sliding windows.\n")
  copynumber$GCperwindow <- getGCinwindows(fullsequence, windowlength)
    
  # Then run sliding window counting on it.
  cat ("Read-location file successfully read. Counting reads in sliding windows...\n")
  copynumber$GC_weights <- NA
  copynumber$is_GC_normalized <- F
  if (!missing(hitsfile)){
    copynumber$ReadsprWindow <- getreadcounts(hits, windowlength, copynumber$chrlength)
    cat ("Successfully counted the number of reads per sliding window!\n")
      
    # Calculating approximate mean and variance counts
    estimates <- sampleChromosome(copynumber)
    copynumber$mean <- estimates[[1]]
    copynumber$variance <- estimates[[2]]
    if (copynumber$mean > copynumber$variance){
      stop("Your data seem to be suffering from underdispersion. Shutting down...\n", call.=FALSE)
    }
  }

  # Remove large objects from memory
  rm(table)
  if (!missing(hitsfile)) {rm(hits)}
  cat("Raw data loaded.\n")
  
  return(copynumber)
}
