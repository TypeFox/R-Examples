"seqaln" <-
function(aln, id=NULL, profile=NULL,
                   exefile = "muscle",
                   outfile = "aln.fa",
                   protein = TRUE,
                   seqgroup = FALSE,
                   refine = FALSE,
                   extra.args = "",
                   verbose = FALSE) {

  ## Log the call
  cl <- match.call()

  ## alignment to fasta object
  aln <- as.fasta(aln, id=id)

  ## nothing to align?
  if(!nrow(aln$ali) > 1 && is.null(profile)) {
    warning("nothing to align")
    aln$ali <- aln$ali[ , !is.gap(aln$ali), drop=FALSE]
    colnames(aln$ali) <- NULL
    return(aln)
  }
  
  if(!is.null(profile) & !inherits(profile, "fasta"))
    stop("profile must be of class 'fasta'")

  if(grepl("clustalo", tolower(exefile))) {
    prg <- "clustalo"
    ver <- "--version"
    
    if(!is.null(profile))
      args <- c("", "--profile1", "--in", "--out")
    else
      args <- c("--in", "--out")

    extra.args <- paste(extra.args,"--force")
    if(seqgroup)
      extra.args <- paste(extra.args, "--output-order=tree-order")
    else
      extra.args <- paste(extra.args, "--output-order=input-order")
    
    if(verbose)
      extra.args <- paste(extra.args,"--verbose")

    if(!is.null(profile) && length(grep("dealign", extra.args))==0)
      warning("profile alignment with clustalo: consider using extra.args='--dealign'")
      
    #if(protein)
    #  extra.args <- paste(extra.args,"--seqtype Protein")
    #else
    #  extra.args <- paste(extra.args,"--seqtype DNA")
  }
  else {
    prg <- "muscle"
    ver <- "-version"
    
    if(!is.null(profile))
      args <- c("-profile", "-in1", "-in2", "-out")
    else
      args <- c("-in", "-out")

    if(refine)
      extra.args <- paste(extra.args,"-refine")
    if(protein)
      extra.args <- paste(extra.args,"-seqtype protein")
    else
      extra.args <- paste(extra.args,"-seqtype dna")
  }
  
  ## Check if the program is executable
  os1 <- .Platform$OS.type
  status <- system(paste(exefile, ver),
                   ignore.stderr = TRUE, ignore.stdout = TRUE)

  if(!(status %in% c(0,1)))
    stop(paste("Launching external program failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))
  

  ## Generate temporary files
  toaln <- tempfile()  
  write.fasta(aln, file=toaln)

  profilealn <- NULL
  if(!is.null(profile)) {
    profilealn <- tempfile()  
    write.fasta(profile, file=profilealn)
  }
  
  if(is.null(outfile))
    fa <- tempfile() 
  else
    fa <- outfile

  ## Build command to external program
  if(is.null(profile)) {
    cmd <- paste(exefile, args[1], toaln, args[2],
                 fa, extra.args, sep=" ")
  }
  else {
    cmd <- paste(exefile, args[1], args[2], profilealn, args[3], toaln, args[4],
                 fa, extra.args, sep=" ")
  }
  
  if(verbose)
    cat(paste("Running command:\n ", cmd , "\n"))
  
  ## Run command
  if (os1 == "windows")
    success <- shell(shQuote(cmd), ignore.stderr = !verbose, ignore.stdout = !verbose)
  else
    success <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)
  
  if(success!=0)
    stop(paste("An error occurred while running command\n '",
               exefile, "'", sep=""))

  ## Re-group sequences to initial alignment order
  ## (muscle groups similar sequences by default)
  naln <- read.fasta(fa, rm.dup=FALSE)
  if(!seqgroup) {
    if(is.null(profile)) {
      ord <- match(aln$id, naln$id)
      naln$id <- naln$id[ord]
      naln$ali <- naln$ali[ord,]
    }
  }

  ## Delete temporary files
  if(!is.null(profile))
    unlink(profilealn)
  unlink(toaln)
  if(is.null(outfile)) unlink(fa)
  naln$call=cl
  return(naln)
}

