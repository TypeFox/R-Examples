"dssp.pdbs" <- function(pdbs, ...) {
  if(!is.pdbs(pdbs))
    stop("provide a pdbs object as obtained from pdbaln()")
  
  dots <- list(...)
  if(any(c("resno", "full") %in% names(dots)))
    stop("arguments resno and full not allowed in dssp.pdbs()")
  
  gaps.res <- gap.inspect(pdbs$ali)
  sse <- matrix(NA, ncol=ncol(pdbs$resno), nrow=nrow(pdbs$resno))
  
  for ( i in 1:length(pdbs$id) ) {
    ##- Check for local and/or online PDB file to run dssp on
    file <- pdbs$id[i]
    toread <- file.exists(file)
    if ((substr(file, 1, 4) == "http") | (nchar(file) == 4)) {
      toread <- TRUE
    }
    if (!toread) {
      stop(paste("Corresponding PDB file could not be found for entry:\n\t-", pdbs$id[i]))
    }
 
    tmp.pdb = read.pdb(pdbs$id[i])
    tmp.sse = dssp.pdb(tmp.pdb, resno=FALSE, full=FALSE, ...)

    ##sse[i,  which(gaps.res$bin[i,]==0)] = tmp.sse$sse

    ##- this old way (line 30) will have problems if alignment contains 
    ##  only a portion of the full PDB structure read on line 27 above. 
    ##  Thus the edit below uses residue number and chain id from the 
    ##  alignment (ali.names) as a reference to populate the sse matrix 
    ##  from the PDB dssp result (pdb.names) 

    ali.names <- paste(pdbs$resno[i,], pdbs$chain[i,], sep="_")

    pdb.names <- paste(tmp.pdb$atom$resno[tmp.pdb$calpha],
                       tmp.pdb$atom$chain[tmp.pdb$calpha], sep="_")
    names(tmp.sse$sse) <- pdb.names
        
    sse[i,] <- tmp.sse$sse[ ali.names ]
  }

  return(sse)
}
