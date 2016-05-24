"pdbs2pdb" <- function(pdbs, inds=NULL, rm.gaps=FALSE) {
  if(!inherits(pdbs, "pdbs")) {
    stop("Input 'pdbs' should be of class 'pdbs', e.g. from pdbaln() or read.fasta.pdb()")
  }
  
  if(is.null(inds))
    inds <- seq(1, length(pdbs$id))

  ## Temporaray file
  fname <- tempfile(fileext = "pdb")

  ## Set indicies
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  
  all.pdbs <- list()
  for ( i in 1:length(inds) ) {
    j <- inds[i]
    
    ## Set indices for this structure only
    f.inds <- NULL
    if(rm.gaps) {
      f.inds$res <- gaps.res$f.inds
      f.inds$pos <- gaps.pos$f.inds
    }
    else {
      f.inds$res <- which(gaps.res$bin[j,]==0)
      f.inds$pos <- atom2xyz(f.inds$res)
    }

    ## Make a temporary PDB object
    write.pdb(pdb=NULL,
              xyz  =pdbs$xyz[j,f.inds$pos],   resno=pdbs$resno[j,f.inds$res],
              resid=pdbs$resid[j,f.inds$res], chain=pdbs$chain[j,f.inds$res],
              file=fname)

    all.pdbs[[i]] <- read.pdb(fname)
  }
  names(all.pdbs) <- sub(".pdb$", "", basename(pdbs$id[inds]))
  return(all.pdbs)
}
