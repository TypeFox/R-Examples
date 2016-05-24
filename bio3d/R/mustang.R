"mustang" <- function(files, exefile="mustang", outfile="aln.mustang.fa",
                      cleanpdb=FALSE, cleandir="mustangpdbs", verbose=TRUE) {
  ## Check if the program is executable
  os1 <- .Platform$OS.type
  status <- system(paste(exefile, "--version"),
                   ignore.stderr = TRUE, ignore.stdout = TRUE)
  
  if(!(status %in% c(0,1)))
    stop(paste("Launching external program failed\n",
               "  make sure '", exefile, "' is in your search path", sep=""))
  
  if(!all(file.exists(files)))
    stop(paste("Missing files:", paste(files[ !file.exists(files) ], collapse=", ")))

  ## produce cleaned CA pdb files for mustang
  if(cleanpdb) {
    if(!file.exists(cleandir))
      dir.create(cleandir)
    newfiles <- c()
    for(i in 1:length(files)) {
      tmpout <- paste(cleandir, basename(files[i]), sep="/")
      pdb     <- read.pdb(files[i])
      sele    <- atom.select(pdb, "calpha", verbose=verbose)
      new     <- trim.pdb(pdb, sele)
      seq1    <- aa321(new$atom$resid)
      seq3     <- aa123(seq1)
      new$atom$type  <- "ATOM"
      new$atom$resid <- seq3
      write.pdb(new, file=tmpout)
      newfiles <- c(newfiles, tmpout)
    }
    files <- newfiles
  }
  
  infile <- tempfile()
  tmpout <- tempfile()
  dirn <- unique(dirname(files))

  if(length(dirn)>1)
    stop("All files must be in one directory")
  
  files <- basename(files)

  rawlines <- NULL
  rawlines <- c(rawlines, paste(">", dirn))

  for ( i in 1:length(files) )
    rawlines <- c(rawlines, paste("+", files[i], sep=""))
  
  write.table(rawlines, file=infile, quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  cmd <- paste(exefile, "-f", infile, "-o", tmpout, "-F fasta")

  if(verbose)
    cat("Running command\n", cmd, "\n")
  
  if (os1 == "windows")
    success <- shell(shQuote(cmd), ignore.stderr = !verbose, ignore.stdout = !verbose)
  else
    success <- system(cmd, ignore.stderr = !verbose, ignore.stdout = !verbose)

  if(success!=0)
    stop(paste("An error occurred while running command\n '",
               exefile, "'", sep=""))
  
  aln <- read.fasta(paste(tmpout, ".afasta", sep=""))
  rownames(aln$ali) <- paste(dirn, rownames(aln$ali), sep="/")
  aln$id <- rownames(aln$ali)

  unlink(infile); unlink(tmpout);

  if(!is.null(outfile))
    write.fasta(aln, file=outfile)
  
  return(aln)
}
